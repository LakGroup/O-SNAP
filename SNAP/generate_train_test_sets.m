function [train_idx, test_idx, split_method] = generate_train_test_sets(feature_data,options)
arguments
    feature_data table
    options.split_method string = "replicate"
    options.k double = 5
    options.proportion double = 0.2 % ratio of test:train
end
    switch options.split_method
        case "replicate"
            % check for proper number of replicates
            bio_reps_by_group=arrayfun(@(x) unique(feature_data{strcmp(feature_data.group,x),"biological_replicate"}), unique(feature_data.group),'uni',0);
            bio_reps_all = bio_reps_by_group{1};
            for i=2:numel(bio_reps_by_group)
                bio_reps_all = intersect(bio_reps_all,bio_reps_by_group{i});
            end
            n_reps = numel(bio_reps_all);
            if n_reps >= 4
                [train_idx, test_idx] = split_for_replicates(feature_data);
                split_method = "replicate";
            else 
                fprintf("  Warning: Number of replicates (=%.0f) < 4, swapping to 'k-fold' option w/ values: 'k'=%.0f\n",n_reps,options.k)
                [train_idx, test_idx] = split_for_k_fold(feature_data,options.k);
                split_method = "k-fold";
            end
        case "bootstrap"
            [train_idx, test_idx] = split_for_bootstrap(feature_data,options.k,options.proportion);
            split_method = "bootstrap";
        case "k-fold"
            [train_idx, test_idx] = split_for_k_fold(feature_data,options.k);
            split_method = "k-fold";
        case "none"
            train_idx = {1:size(feature_data,1)};
            test_idx = [];
        otherwise
            ME = MException('SNAP:invalid_argument_for_feature_selection', ...
                sprintf("Argument 'by' is invalid: %s\nChoose from: 'replicate' (default),'bootstrap','k-fold','none'",options.split_method));
            throw(ME)
    end
end

function [train_idx, test_idx] = split_for_replicates(data)
    [bio_reps,~,idx] = unique(data.biological_replicate);
    n_reps = numel(bio_reps);
    test_idx = cell(1,n_reps);
    train_idx = cell(1,n_reps);
    for i=1:n_reps
        test_idx{i} = find(i==idx);
        train_idx{i} = find(i~=idx);
    end
end

function [train_idx, test_idx] = split_for_k_fold(data,k)
    n_samps = size(data,1);
    n_groups = numel(unique(data.group));
    all_groups_present = false;
    % split data until it is ensured that there are multiple groups in each division
    tries = 1;
    while ~all_groups_present && tries < 30
        split_idx_components = split_data_to_n(1:n_samps,k);
        all_groups_present = all(cellfun(@(x) numel(unique(data(x,"group"))),split_idx_components,'uni',1) == n_groups);
        tries = tries + 1;
    end
    test_idx = cell(1,k);
    train_idx = cell(1,k);
    for i=1:k
        test_idx{i} = split_idx_components{i};
        train_idx{i} = [split_idx_components{(1:k)~=i}];
    end
end

function [train_idx, test_idx] = split_for_bootstrap(data,k,proportion)
    n_samps = size(data,1);
    n_groups = numel(unique(data.group));
    n_test = floor(proportion*n_samps);
    if n_test < n_groups
        ME = MException('SNAP:bootstrap_proportion_low', ...
                sprintf("W/ option 'proportion'=%.2f', bootstrap test set n=%.0f < number of groups n=%.0f",proportion,n_test,n_groups));
        throw(ME)
    end
    test_idx = cell(1,k);
    train_idx = cell(1,k);
    for i=1:k
        all_groups_present = false;
        % split data until it is ensured that there are multiple groups in each division
        tries = 1;
        while ~all_groups_present && tries < 30
            idx = randperm(n_samps);
            idx_cutoff = n_test;
            idx_test = idx(1:idx_cutoff);
            idx_train = idx(idx_cutoff+1:end);
            all_groups_present = all([numel(unique(data(idx_test,"group"))) numel(unique(data(idx_train,"group")))] == n_groups);
            tries = tries + 1;
        end
        test_idx{i} = idx_test;
        train_idx{i} = idx_train;
    end
end