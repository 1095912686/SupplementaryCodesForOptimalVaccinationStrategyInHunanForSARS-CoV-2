function caseCountInGroups = findCaseCountInGroups(data, agePartition)

groupCount = numel(agePartition) - 1;
[~, uniqueCaseNameIndex] = unique(data.caseName);

groupOfCase = whichGroup(data.caseAge(uniqueCaseNameIndex), agePartition);

caseCountInGroups = zeros(groupCount, 1);
for i = 1:groupCount    % inspect cases in each age group
    idx = groupOfCase == i; % index of cases in group i
    caseCountInGroups(i) = sum(idx);  % case count in group i
end

end