function copy_labels_limits(ax1, ax2)
    ax1.XLabel.String = ax2.XLabel.String;
    ax1.YLabel.String = ax2.YLabel.String;
    ax1.ZLabel.String = ax2.ZLabel.String;
    ax1.XLim = ax2.XLim;
    ax1.YLim = ax2.YLim;
    ax1.ZLim = ax2.ZLim;
end

