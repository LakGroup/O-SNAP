function validIdx = idxVoronoiAreaDataContainingVar(dataInfoTable,varsToLoad)
    filepaths = dataInfoTable{:,'filepath'};
    validIdx = ones(size(dataInfoTable,1),1);
    for i=1:length(filepaths)
        % disp(dataInfoTable{i,'filepath'});
        try
            if ~hasVars(dataInfoTable{i,'filepath'},varsToLoad)
                validIdx(i) = 0;
            end
        catch
            validIdx(i) = 0;
        end    
    end
    if ~all(validIdx)
        disp('Problem files:')
        f = dataInfoTable{~validIdx,'filepath'};
        for i=1:length(f) 
            disp("   " + num2str(i) + "   " + f(i));
        end
    end
    validIdx = find(validIdx);
end