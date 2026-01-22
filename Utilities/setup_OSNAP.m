% -------------------------------------------------------------------------
% setup_OSNAP.m
% -------------------------------------------------------------------------
% Sets up environment for utilities that OSNAP depends on. Ideally, only
% needs to be run once.
%
% Example on how to use it:
%   setup_OSNAP
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%%
function flag = setup_OSNAP()
    %% Check that OSNAP is in the path
    disp("Checking path...")
    check_for_OSNAP_path()
    %% Check for all required Toolboxes
    disp("Checking for required toolboxes...")
    check_for_toolboxes();
    %% Ensure correct functions if shadowed another installation (ex. PLS toolbox)
    disp("Ensuring functions are not shadowed...")
    top_in_search_path_OSNAP('skewness','stats\stats\skewness.m');
    top_in_search_path_OSNAP('pca','stats\stats\pca.m');
    savepath
    flag = 1;
    %% start parallel pool
    p_pool = gcp('nocreate');
    if isempty(p_pool)
        parpool("Processes");
    elseif class(p_pool) == "parallel.ThreadPool"
        delete(p_pool)
        parpool("Processes");
    end
end

function check_for_OSNAP_path()
    if isempty(which('run_OSNAP_pipeline'))
        path_to_OSNAP = uigetdir(getenv('USERPROFILE'),"Select directory and subdirectories containing O-SNAP package");
        addpath(genpath(fullfile(path_to_OSNAP,'Utilities')));
    end
end

function check_for_toolboxes()
    toolbox_req = {...
        'Bioinformatics Toolbox',...
        'Deep Learning Toolbox',...
        'Parallel Computing Toolbox',...
        'Signal Processing Toolbox',...
        'Statistics and Machine Learning Toolbox',...
        };
    toolbox_list = ver;
    has_toolboxes = zeros(1,numel(toolbox_req));
    for i=1:numel(toolbox_list)
        [flag, idx] = ismember(toolbox_list(i).Name,toolbox_req);
        if flag
            has_toolboxes(idx)=1;
        end
    end
    if any(~has_toolboxes)
        message = "Missing the following MATLAB Toolboxes:";
        toolbox_missing = toolbox_req(~has_toolboxes);
        for i=1:sum(~has_toolboxes)
            message_append = sprintf("\n  - %s",toolbox_missing{i});
            message = message + message_append;
        end
        message = message + sprintf("\nPlease install Toolboxes (Home > Add-Ons > Manage Add-Ons)");
        error(message,'OSNAP:missingToolboxes');
    end
end