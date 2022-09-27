function nC = caseCountViaData()
ImportAgePartition;
caseContactData = readtable('contactDataProcessed.xlsx');
[uniqueCases, uniqueID] = unique(caseContactData.caseNameInCaseData);

caseAge = caseContactData.caseAge(uniqueID);
vaccines = caseContactData(uniqueID,30:end);

nC = zeros(numel(agePartition)-1, 1);
for i = 1:size(uniqueCases,1)
    
    % skip the missing contactAges
    if isnan(caseAge(i))
        continue;
    end
        
%     % how many doses this contact has taken
%     dosei = ~isempty(vaccines.time1{i}) + ~isempty(vaccines.time2{i}) + ~isempty(vaccines.time3{i});
    
    % which group this contact belongs to 
    groupi = whichGroup(caseAge(i), agePartition);
    
    % add this contact to corresponding age-stratified convergence table
    nC(groupi, 1) = nC(groupi, 1) + 1;
end

writematrix(nC,'nC.xlsx');
end