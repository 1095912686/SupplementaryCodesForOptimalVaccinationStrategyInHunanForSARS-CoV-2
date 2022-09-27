function [nVC,nVCC] = vaccineCoverageViaDataByIsComplete()

% read data
ImportAgePartition;

opts = detectImportOptions('contactDataProcessed.xlsx');
[opts.VariableTypes{[7:8, 28, 31, 34]}] = deal('datetime', 'datetime', 'datetime', 'datetime', 'datetime'); % illnessOnset, diagnosedDate, time1,2,3
caseContactData = readtable('contactDataProcessed.xlsx', opts);

caseNames = unique(caseContactData.caseNameInCaseData); 
caseNames(1) = []; % delete the empty
contactNames = caseContactData.contactName;
contactAge = caseContactData.contactAge;
vaccines = caseContactData(:,27:end);

% main loop
nVC = zeros(numel(agePartition)-1, 4);
nVCC = zeros(numel(agePartition)-1, 4);
for i = 1:size(caseContactData,1)
    
    if i == 7882
        1;
    end

    % skip the missing contactAges
    if isnan(contactAge(i))
        continue;
    end

    % how many doses this contact has taken
    doseNum = ~isnat(vaccines.time1(i)) + ~isnat(vaccines.time2(i)) + ~isnat(vaccines.time3(i));
    doseDemand = max([vaccines.doseDemand1(i), vaccines.doseDemand2(i), vaccines.doseDemand3(i)]);

    if doseNum == 0
        dosei = 0;    % none
    elseif doseNum < doseDemand 
        dosei = 1;    % unfully vaccinated
    elseif doseNum == doseDemand
        dosei = 2;    % fully vaccinated
    elseif doseNum > doseDemand
        dosei = 3;    % booster vaccinated
    end 

    % which group this contact belongs to
    groupi = whichGroup(contactAge(i), agePartition);

    % add this contact to corresponding age-stratified coverage table
    nVC(groupi, dosei + 1) = nVC(groupi, dosei + 1) + 1;

    % if the contact is a case, save to nVCC
    if sum(strcmp(contactNames(i), caseNames)) >= 1
        nVCC(groupi, dosei + 1) = nVCC(groupi, dosei + 1) + 1 / sum(strcmp(contactNames,contactNames(i)));
    end
end

writematrix(nVC,'nVC.xlsx');
writematrix(nVCC,'nVCC.xlsx');

end