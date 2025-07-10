% -------------------------------------------------------------------------
% .m
% -------------------------------------------------------------------------
% 
%
% Example on how to use it:
%   
% -------------------------------------------------------------------------
% Input:
%   : 
%   : 
%   : 
%   : 
% Output:
%   : 
% Options:
%   : 
%   : 
%   : 
%   : 
%   : 
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%%%% generate volcano plot
function plot_OSNAP_volcano(feature_compare_data,group_1,group_2,options)
    arguments
        feature_compare_data table;
        group_1 string;
        group_2 string;
        options.alpha double = 0.05;
        options.fold_change_threshold double = 2;
        % options.plot_labels logical = false;
        options.save_path string = "";
    end
    log2_fold_change_threshold = log2(options.fold_change_threshold);
    idx_sig = all([abs(feature_compare_data.log2_fold_change)>log2_fold_change_threshold feature_compare_data.adj_p_value < options.alpha],2);
    idx_sig_up = all([feature_compare_data.log2_fold_change>0 idx_sig],2);
    idx_sig_down = all([feature_compare_data.log2_fold_change<0 idx_sig],2);
    v_fig_size = 20*ones(size(feature_compare_data,1),1);
    v_fig_size(idx_sig) = 40*ones(sum(idx_sig),1);
    v_fig_color = repmat([0.5 0.5 .5],size(feature_compare_data,1),1);
    v_fig_color(idx_sig_up,:) = repmat([1 0 0],sum(idx_sig_up),1);
    v_fig_color(idx_sig_down,:) = repmat([0 0 1],sum(idx_sig_down),1);
    v_fig = figure('visible','off');
    scatter(feature_compare_data.log2_fold_change,-log10(feature_compare_data.adj_p_value),v_fig_size,v_fig_color,'filled','MarkerEdgeColor','k');
    hold on
    plot([log2_fold_change_threshold log2_fold_change_threshold],ylim(),'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
    hold on
    plot(-[log2_fold_change_threshold log2_fold_change_threshold],ylim(),'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
    hold on
    plot(xlim(), [-log10(options.alpha) -log10(options.alpha)],'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
    ylabel("-log_{10}(adj. p-value)","FontSize",18)
    xlabel("log_{2}(Fold Change)","FontSize",18)
    title(replace("\color{blue}" + group_1 + " \color{black}vs " + "\color{red}" + group_2,"_","-"),...
        'FontSize',20);
    set(gca,'FontSize',14)
    if options.save_path ~= ""
        save_path = sprintf("%s_%s_%s",options.save_path,group_1,group_2);
        savefig(save_path+".fig");
        saveas(v_fig,save_path+".png");
    end
end

