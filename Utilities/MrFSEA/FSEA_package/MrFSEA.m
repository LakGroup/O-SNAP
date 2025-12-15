function [res_pos,res_neg,res_descr,p_feature] = MrFSEA(work_dir,data,group,data_entrez,file_name,opts)
% Computing Feature Set Enrichment Analysis (FSEA) with various ranking metrics.
% Major adjustments:
% - 12 feature ranking metrics implemented (plus its ranks and/or absolute values)
% - possibility to use external ranking metric method
% - possibility to use external Feature Set Database
% - phenotype permutation available for all ranking methods (even external)
% - increased efficiency by parallel computing
%
% Input:
% data - feature expression matrix [mxn]
% group - sample labels [1xn]
% data_entrez - feature entrez IDs [mx1]
% file_name - save file name 
% opts - structure of parameters
% Output:
% res_pos - table with positive FSEA results
% res_neg - table with negative FSEA results
% res_desc - description of tables with results
% p_feature - features p-values from permutation test
%
% Author:
% Michal.Marczyk@polsl.pl

%check inputs
if nargin < 5
    opts = default_FSEA_opts();
end

if nargin < 3
    error('Too few input arguments.')
end

if length(unique(data_entrez)) ~= length(data_entrez)
    error('Non-unique feature identifiers found.')
end

%load FS database
if opts.show; disp(['Load ' char(opts.FS_name) ' database.']); end
tmp = load(opts.FS_name+".mat");
FS = tmp.FS;
[~,sign_sort] = sort(abs([FS.entrez{:}]));
opts.direction = [FS.entrez_dir{:}];
opts.direction = opts.direction(sign_sort);

%calculate ranks
[data_ranks,opts.tail] = calc_ranks_FSEA(data,group,opts.rank_type,opts.abs,opts.tied_rank,opts.direction);

%sort data by ranks
[data_ranks,ind] = sort(data_ranks,opts.sort_type);
data_entrez = data_entrez(ind);
data = data(ind,:);
n = length(data_ranks);

%filter database by removing features in FS that are not present in data
ind_d = cell(1,FS.nb);  %indices of all FS features in dataset
RI = ind_d;
for a=1:FS.nb
    [common,ind_d{a}] = intersect_fast(data_entrez,FS.entrez{a});  %find features that are in given FS
    if length(common) ~= FS.entrez_nb(a)
        FS.entrez{a} = common;
        FS.entrez_dir{a} = FS.entrez_dir{a}(ind_d{a});
        FS.entrez_nb(a) = length(common);    %find features with ranks within FS
    end
    RI{a} = data_ranks(ind_d{a});
end
    
%filter database by min. number of features in FS
if opts.FS_filt(1) > 0 || opts.FS_filt(2) < Inf
    del = FS.entrez_nb < opts.FS_filt(1) | FS.entrez_nb > opts.FS_filt(2);
    FS.ID(del) = [];
    FS.descr(del) = [];
    FS.nb = length(FS.ID);
    FS.entrez(del) = [];
    FS.entrez_dir(del) = [];
    FS.entrez_nb(del) = [];
    ind_d(del) = [];
    RI(del) = [];
    if opts.show; disp([num2str(sum(del)) ' FS filtered, ' num2str(FS.nb) ' FS left in database.']); end
else
    if opts.show; disp([num2str(FS.nb) ' FS in database found.']); end
end

%pre-calculate ranks for permutation test
if opts.show
    disp('Calculate ranking metrics for all permutations.')
end
ranks_perm = pre_calc_ranks_FSEA(data_ranks,data,group,opts);
switch opts.tail
    case 'both'
        p_feature = 1-(sum(repmat(abs(data_ranks),1,opts.perm_nb) > ranks_perm,2))/opts.perm_nb;
    case 'left'
        p_feature = 1-(sum(repmat(data_ranks,1,opts.perm_nb) < ranks_perm,2))/opts.perm_nb;
    case 'right'
        p_feature = 1-sum(repmat(data_ranks,1,opts.perm_nb) > ranks_perm,2)/opts.perm_nb;
end
p_feature(ind) = max(p_feature,1/opts.perm_nb);

% separate calculations for each FeatureSet
out = zeros(FS.nb,9); NES_perm = zeros(FS.nb,opts.perm_nb);
LEF = cell(FS.nb,1);
if opts.show
    fprintf('Progress:\n');
    fprintf(['\n' repmat('.',1,FS.nb) '\n\n']);
end
for a=1:FS.nb 
    opts_tmp = opts;
    out_tmp = zeros(1,9);
    FS_tmp = FS;
        
    %sort ranks in FS to find LEF
    RI_tmp = sort(RI{a},opts_tmp.sort_type); 
    
    %Enrichment Score calculation
    Phit = zeros(n,1); 
    Pmiss = ones(n,1) * 1/(n-FS_tmp.entrez_nb(a));
    Phit(ind_d{a}) = (abs(RI{a}).^opts_tmp.p)/sum(abs(RI{a}));
    Pmiss(ind_d{a}) = 0;
    ES = cumsum(Phit - Pmiss);    
    [~,ind] = max(abs(ES));     %find maximum ES 
    ESmax = ES(ind);
    
    if opts.show
        plot_classic_FSEA(replace(file_name,'FSEA_','FSEA_plots_'),ES,ESmax,data_ranks,ind_d{a},FS.ID{a},Phit,Pmiss);
    end
    
    % Find p-value by permutation test
    ES_perm = ES_perm_test(ranks_perm,data_entrez,FS_tmp.entrez{a},opts_tmp.perm_nb,opts_tmp.sort_type);
        
    out_tmp(1) = a;        % FS name  
    out_tmp(2) = FS_tmp.entrez_nb(a);     % FS Size  
    out_tmp(3) = ESmax;           % Enrichment Score 
    
    if ESmax>0    %positive ES
        [LEF_tmp,ind_LEF] = intersect_fast(FS_tmp.entrez{a},data_entrez(1:ind));   % features from LEF 
        LEF_tmp(:,2) = RI_tmp(ind_LEF);      % ES from features of LEF   
        ind = ES_perm > 0;
        tmp_NES = zeros(1,opts.perm_nb);
        tmp_NES(ind) = ES_perm(ind)/mean(ES_perm(ind));
        NES_perm(a,:) = tmp_NES;
        if ~isempty(ES_perm(ind))
            out_tmp(4) = max(1,sum(ES_perm >= ESmax))/length(ES_perm(ind));  % p-value
            out_tmp(5) = ESmax/mean(ES_perm(ind)); % Normalized Enrichment Score 
        else
            out_tmp(4) = NaN;  % p-value
            out_tmp(5) = NaN; % Normalized Enrichment Score 
        end
        
    else    %negative ES
        [LEF_tmp,ind_LEF] = intersect_fast(FS_tmp.entrez{a},data_entrez(ind:end));   % features from LEF 
        LEF_tmp(:,2) = RI_tmp(ind_LEF);      % ES from features of LEF   
        ind = ES_perm < 0;
        tmp_NES = zeros(1,opts.perm_nb);
        tmp_NES(ind) = ES_perm(ind)/abs(mean(ES_perm(ind)));
        NES_perm(a,:) = tmp_NES;
        if ~isempty(ES_perm(ind))
            out_tmp(4) = max(1,sum(ES_perm <= ESmax))/length(ES_perm(ind));  % p-value
            out_tmp(5) = ESmax/abs(mean(ES_perm(ind))); % Normalized Enrichment Score 
        else
            out_tmp(4) = NaN;  % p-value
            out_tmp(5) = NaN; % Normalized Enrichment Score 
        end
    end
    out_tmp(9) = size(LEF_tmp,1);          % observed counts in LEF
    
    out(a,:) = out_tmp;
    LEF{a} = LEF_tmp;
    if opts.show
        fprintf('\b|\n');
    end
end

%calculate NES p-values using NES null distribution for pos and neg
%separately
NES_perm_pos = NES_perm(NES_perm>0);
NES_perm_neg = NES_perm(NES_perm<0);
ind_pos = out(:,3) >= 0;
for a=1:FS.nb
    if ind_pos(a)
        out(a,6) = max(sum(NES_perm_pos >= out(a,5))/length(NES_perm_pos),1/length(NES_perm_pos));          %NES p-value
        out(a,7) = min(out(a,6)/(sum(out(ind_pos,5) >= out(a,5))/sum(ind_pos)),1);%NES q-value
    else
        out(a,6) = max(sum(NES_perm_neg <= out(a,5))/length(NES_perm_neg),1/length(NES_perm_neg));
        out(a,7) = min(out(a,6)/(sum(out(~ind_pos,5) <= out(a,5))/sum(~ind_pos)),1);
    end
end

out(ind_pos,8) = min(1,sum(ind_pos)*out(ind_pos,6)); % NES FWER
out(~ind_pos,8) = min(1,sum(~ind_pos)*out(~ind_pos,6));

%divide results into pos and neg
res_pos = out(ind_pos,:);
LEF_pos = LEF(ind_pos);
res_neg = out(~ind_pos,:);
LEF_neg = LEF(~ind_pos);

%sort result by FS significance
switch opts.FS_sort_type
    case 'ES_pval'
        [~,ind1] = sort(res_pos(:,4));
        [~,ind2] = sort(res_neg(:,4));
    case 'NES_qval'
        [~,ind1] = sort(res_pos(:,6));
        [~,ind2] = sort(res_neg(:,6));
    case 'NES_FWER'
        [~,ind1] = sort(res_pos(:,7));
        [~,ind2] = sort(res_neg(:,7));
    case 'none'
        ind1 = 1:size(res_pos,1);
        ind2 = 1:size(res_neg,1);
    otherwise
        error('Wrong FS sorting method selected.')
end
res_pos = res_pos(ind1,:);
ind_tmp = res_pos(:,1);
res_pos = num2cell(res_pos);
res_pos(:,1) = FS.ID(ind_tmp);

res_neg = res_neg(ind2,:);
ind_tmp = res_neg(:,1);
res_neg = num2cell(res_neg);
res_neg(:,1) = FS.ID(ind_tmp);
res_descr={'FS.NAME','FS size','ES','ES_p-val','NES','NES_p-val','NES_q-val','NES_FWER','# features in LEF','LEF features'};

if opts.save
    warning off
    
    % save results in .mat file
    save(fullfile(work_dir,[file_name '.mat']),'res_neg','res_pos','res_descr','LEF_neg','LEF_pos')
    
    if ispc
        % save results in .xls file
        xlswrite(fullfile(work_dir,file_name),res_descr,'pos','A1')
        if ~isempty(res_pos)
            xlswrite(fullfile(work_dir,file_name),res_pos,'pos','A2');
            xlswrite(fullfile(work_dir,file_name),LEF_pos,'pos','J2')
        end
        xlswrite(fullfile(work_dir,file_name),res_descr,'neg','A1')
        if ~isempty(res_neg)
            xlswrite(fullfile(work_dir,file_name),res_neg,'neg','A2');
            xlswrite(fullfile(work_dir,file_name),LEF_neg,'neg','J2')
        end
    else
        disp('Saving to .xls on that platform is not supported.')
    end
end