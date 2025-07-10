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
function setup_OSNAP()
try
    % TODO: add instructions directing user to all relevant installations
    path_to_insidepoly = uigetdir(getenv('USERPROFILE'),"Please select directory containing inside_poly");
    old_folder = cd(path_to_insidepoly);
    insidepoly_install
    cd(old_folder)
    % Ensure correct functions if shadowed by the PLS toolbox installation
    top_in_search_path_OSNAP('skewness','stats\stats\skewness.m');
    top_in_search_path_OSNAP('pca','stats\stats\pca.m');
catch ME
    if(strcmp(ME.identifier,'MATLAB:mex:No_compiler_found_link_Win64'))
        msgbox({'Install C compiler for Matlab:',...
            '   1. On the MATLAB Home tab, in the Environment section, click Add-Ons > Get Add-Ons.', ...
            '   2. Search for Min_gW or select from Features.',...
            '   3. Follow instructions on page to install Min_gW-w64 C/C++/Fortran Compiler.',...
            '      Note: Different MATLAB versions may have different installation steps.',...
            '   4. Restart Matlab.'})
    end
    rethrow(ME)
end
end
