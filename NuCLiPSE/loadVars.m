% load specified variables from file
function data = loadVars(filePath, varsToLoad)
    idx = cellfun(@(x) ~isempty(who('-file', filePath,x)),varsToLoad,'uni',1);
    varsToLoadInData = varsToLoad(idx);
    data = load(filePath,varsToLoadInData{:});
    % varsNotThere = varsToLoad(~idx);
end