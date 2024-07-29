% Fast way to check file for variables
function result = hasVars(filepath, varsToCheck)
    idx = cellfun(@(x) ~isempty(who('-file', filepath,x)),varsToCheck,'uni',1);
    result = all(idx);
end