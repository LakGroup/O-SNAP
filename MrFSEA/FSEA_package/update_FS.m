function update_FS(data_name)
% Function for reading FeatureSet databaset from .xls file and storing as .mat
% file.
% Input:
% data_name - filename of FS database with filename extension (.xls,.xlsx, etc.)

dot_ind = find(data_name=='.');
if isempty(dot_ind)
    error('No file extension included.')
end

%load data
tmp1 = readtable(data_name,"NumHeaderLines",0);
disp([data_name(1:dot_ind-1) ' dataset loaded.'])

FS.ID = tmp1{:,1};        %FeatureSet ID
FS.descr = tmp1{:,2};       %FeatureSet descriptions
FS.nb = length(FS.ID);    %Number of FS in databaset

FS.entrez = cell(FS.nb,1);      %Feature entrez ID within each FS
FS.entrez_nb = zeros(FS.nb,1);  %Number of features within each FS
for a=1:FS.nb
    % disp(['FS no. ' num2str(a) '/' num2str(FS.nb)])
    tmp_a = tmp1{a,3:end};
    FS.entrez{a} = unique(abs(tmp_a(~isnan(tmp_a))));
    FS.entrez_dir{a} = sign(unique(tmp_a(~isnan(tmp_a))));
    FS.entrez_nb(a) = length(FS.entrez{a});
    
    if isnumeric(FS.ID{a})
        FS.ID{a} = num2str(FS.ID{a});
    end        
end

save([data_name(1:dot_ind-1) '.mat'],'FS')