function ranks_perm = pre_calc_ranks(data_ranks,data,group,opts)
% Calculate ranks for permutation test.

[m,n] = size(data);
rank_type = opts.rank_type;
if_abs = opts.abs;
if_trank = opts.tied_rank;
direction = opts.direction;

data_ranks_p = cell(1,maxNumCompThreads);
data_ranks_p(:) = {data_ranks};
group_p = cell(1,maxNumCompThreads);
group_p(:) = {group};
perm_nb_p = num2cell(get_n_samples_per_split(opts.perm_nb,maxNumCompThreads));
ranks_perm_p = cellfun(@(a) zeros(m,a), perm_nb_p,'uni',0);

if strcmp(opts.perm_type,'entrez')
    parfor p=1:maxNumCompThreads
        for a=1:perm_nb_p{p}
            ind = RandPermFast(m);
            ranks_perm_p{p}(:,a) = data_ranks_p{p}(ind);
        end
    end
elseif strcmp(opts.perm_type,'pheno')
    parfor p=1:maxNumCompThreads
        for a=1:perm_nb_p{p}
            ind = RandPermFast(n);
            ranks_perm_p{p}(:,a) = calc_ranks(data,group_p{p}(ind),rank_type,if_abs,if_trank,direction);
        end
    end
end
ranks_perm = [ranks_perm_p{:}];