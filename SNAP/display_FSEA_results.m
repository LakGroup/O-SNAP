function display_FSEA_results(T_FSEA, group_pair, options)
arguments
    T_FSEA table
    group_pair string
    options.alpha double = 0.05
end
    families_significant = T_FSEA.("NES_q-val") < options.alpha;
    if any(families_significant)
        disp(group_pair(1) + " VS " + group_pair(2))
        disp("- - - - - - - - - - - - - - - - -")
        T_FSEA_filt = sortrows(T_FSEA(families_significant,:),"NES_q-val");
        for j=1:sum(families_significant)
            fprintf("%s:\t%4.2f\t%4.2e\n\n",string(T_FSEA_filt{j,"FS.NAME"}),T_FSEA_filt{j,"NES"},T_FSEA_filt{j,"NES_q-val"});
        end
        disp("- - - - - - - - - - - - - - - - -")
    end
end

