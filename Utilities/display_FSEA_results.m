% -------------------------------------------------------------------------
% display_FSEA_results.m
% -------------------------------------------------------------------------
% Prints the results from FSEA analysis for a given pair of phenotypes
%
% Example on how to use it:
%   display_FSEA_results(T_FSEA, group_pair, options)
% -------------------------------------------------------------------------
% Input:
%   T_FSEA: The output from FSEA with information on feature family ES,
%           NES, and statistics
%   group_pair: [1x2] array containing group identifiers
% Options:
%   alpha: Significance testing threshold
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   H. H. Kim, J. A. Martinez-Sarmiento, F. R. Palma, A. Kant, E. Y. Zhang,
%   Z. Guo, R. L. Mauck, S. C. Heo, V. Shenoy, M. G. Bonini, M. Lakadamyali,
%   O-SNAP: A comprehensive pipeline for spatial profiling of chromatin
%   architecture. bioRxiv, doi: 10.1101/2025.07.18.665612 (2025).
% -------------------------------------------------------------------------
%%
function display_FSEA_results(T_FSEA, group_pair, options)
arguments
    T_FSEA table
    group_pair string
    options.alpha double = 0.05
end
    % Filter for values satisfying alpha
    families_significant = T_FSEA.("NES_q-val") < options.alpha;
    % Print to screen results
    if any(families_significant)
        disp(group_pair(1) + " VS " + group_pair(2))
        disp("- - - - - - - - - - - - - - - - -")
        T_FSEA_filt = sortrows(T_FSEA(families_significant,:),"NES_q-val");
        for j=1:sum(families_significant)
            fprintf("%s:\t%4.2f\t%4.2e\n\n",...
                string(T_FSEA_filt{j,"FS.NAME"}),...
                T_FSEA_filt{j,"NES"},...
                T_FSEA_filt{j,"NES_q-val"});
        end
        disp("- - - - - - - - - - - - - - - - -")
    end
end


