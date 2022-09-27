clear; close all; clc;

opts = detectImportOptions('contactDataProcessed.xlsx');
opts.VariableTypes(1,[7,8]) = repmat({'datetime'},1,2);
data = readtable('contactDataProcessed.xlsx',opts);

% caseNames
caseNames = data.caseNameInCaseData;
emptyID = cellfun(@isempty,caseNames);
data(emptyID,:) = [];

caseNames = data.caseNameInCaseData;
[uniqueCaseNames,uniqueCaseIndex] = unique(caseNames);

% diagnoseDate
diagnosedDate = data.diagnosedDate(uniqueCaseIndex,1);

minDate = dateshift(min(diagnosedDate),'start','day')
maxDate = dateshift(max(diagnosedDate),'end','day')
edges = minDate:maxDate;
histogram(diagnosedDate, edges);
%set(gca,'xticklabel',['', sortedNames(id)]);


% age partition
agePartition = [0,10,20,30,40,50,60,70,200]; %linspace(0,100,11); % [0,10,20,...,100]

% read agePopulation
agePopulation = readmatrix('agePopulationVector.xlsx');

% case-wise duration of making contacts set as 4-day constant according to data
duration = 4 * ones(size(data, 1));  

%% main loop, find the segment partition that maximizing the F norm of difference of two contact matrices
timeVector = dateshift(min(diagnosedDate),'end','day') : dateshift(max(diagnosedDate),'start','day');
d = numel(timeVector);
Fnorm = zeros(d,1);
for i = 1:d
    tau = timeVector(i);
    data1 = data(data.diagnosedDate <= tau,:);
    data2 = data(data.diagnosedDate > tau,:);

     caseCountInGroups1 = findCaseCountInGroups(data1, agePartition);
     caseCountInGroups2 = findCaseCountInGroups(data2, agePartition);

    contactDataMatrix_1 = computeContactDataMatrix(data1, duration, data1.caseAge, data1.contactAge, agePartition);
    contactDataMatrix_2 = computeContactDataMatrix(data2, duration, data2.caseAge, data2.contactAge, agePartition);
    C1_MLE = estimateContactMatrix_MLE(contactDataMatrix_1, caseCountInGroups1, agePopulation, 0);
    C2_MLE = estimateContactMatrix_MLE(contactDataMatrix_2, caseCountInGroups2, agePopulation, 0);
    
    Fnorm(i) = sum(abs(C1_MLE(:) - C2_MLE(:)), 'omitnan');
end

plot(Fnorm, 'rp');
xlabel('Days since July 29, 2021');
ylabel('norm(C_1 - C_2, ''fro'')');
title('Difference between C_1, C_2 for all possible date of partitioning');

