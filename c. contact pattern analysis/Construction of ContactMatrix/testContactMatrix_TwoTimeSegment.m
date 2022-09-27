clear; close all; clc;

ImportAgePartition;
ImportFigureLegends;
opts = detectImportOptions('contactDataProcessed.xlsx');
opts.VariableTypes(1,[7,8]) = repmat({'datetime'},1,2);
caseContactData = readtable('contactDataProcessed.xlsx',opts);

%  extract data with (caseAge, contactAge) available 
caseContactData = caseContactData(~isnan(caseContactData.caseAge) & ~isnan(caseContactData.contactAge), :); 
vaildCaseContactAgePairsCount = size(caseContactData, 1)  
[uniqueCaseName, uniqueCaseNameIndex] = unique(caseContactData.caseName);
uniqueCaseNameCount = numel(uniqueCaseName) 


% caseCount in age groups
caseCountInGroups = findCaseCountInGroups(caseContactData, agePartition)


% caseCountInGroups = zeros(8, 1);
% groupOfCase = whichGroup(caseContactData.caseAge(uniqueCaseNameIndex), agePartition);
% for i = 1:8    % inspect cases in each age group
%     idx = groupOfCase == i; % index of cases in group i
%     caseCountInGroups(i) = sum(idx);  % case count in group i
% end


% diagnoseDate
diagnosedDate = caseContactData.diagnosedDate(uniqueCaseNameIndex,1);

minDate = dateshift(min(diagnosedDate),'start','day')
maxDate = dateshift(max(diagnosedDate),'end','day')
edges = minDate:maxDate;
histogram(diagnosedDate, edges);
%set(gca,'xticklabel',['', sortedNames(id)]);

plot(diagnosedDate, caseContactData.caseAge(uniqueCaseNameIndex), 'rp');

% datetime partition threshold
 tau = datetime(2021,08,04);
%tau = datetime(2021,08,07);
% tau = datetime(2021,08,01);
% tau = datetime(2021,08,14);
data1 = caseContactData(caseContactData.diagnosedDate <= tau,:);
data2 = caseContactData(caseContactData.diagnosedDate > tau,:);
caseCountInGroups1 = findCaseCountInGroups(data1, agePartition)
caseCountInGroups2 = findCaseCountInGroups(data2, agePartition)

duration = 4 * ones(size(caseContactData, 1));  % case-wise duration of making contacts set as 4-day constant according to data

contactDataMatrix = computeContactDataMatrix(caseContactData, duration, caseContactData.caseAge, caseContactData.contactAge, agePartition, caseCountInGroups);
contactDataMatrix_2 = computeContactDataMatrix(data2, duration, data2.caseAge, data2.contactAge, agePartition, caseCountInGroups2);
contactDataMatrix_1 = computeContactDataMatrix(data1, duration, data1.caseAge, data1.contactAge, agePartition, caseCountInGroups1);


%% read agePopulation
agePopulation = readmatrix('agePopulationVector.xlsx');


%% maximum likelihood estimation of contact matrix C_1, C_2
C1_MLE = estimateContactMatrix_MLE(contactDataMatrix_1, caseCountInGroups1, agePopulation, 0);
C2_MLE = estimateContactMatrix_MLE(contactDataMatrix_2, caseCountInGroups2, agePopulation, 0);

figure;
heatmap(double(C1_MLE>C2_MLE));
sum(C1_MLE>C2_MLE,'all')

C1_MLE(isnan(C1_MLE)) = 0;
lambdaMax.C1 = max(eig(C1_MLE));
lambdaMax.C2 = max(eig(C2_MLE));
lambdaMax

%% plot and export
matrixplot(C1_MLE, 'contactMatrix of first segment, estimated by MLE', ageLegend);  % heatmap of contactDataMatrix
exportgraphics(gca,'contactMatrixSegment1_MLE.jpg','Resolution',300); % export figure

matrixplot(C2_MLE, 'contactMatrix of second segment, estimated by MLE', ageLegend);  % heatmap of contactDataMatrix
exportgraphics(gca,'contactMatrixSegment2_MLE.jpg','Resolution',300); % export figure


figure;
plot(C1_MLE(:),C2_MLE(:),'rp'); hold on;
t = linspace(0,12,1e3);

% linearEq = fit(C1_MLE(:), C2_MLE(:), 'linear');
% k = predict(linearEq,t);

plot(t,t,'b-');


function caseCountInGroups = findCaseCountInGroups(data, agePartition)
groupCount = numel(agePartition) - 1;
caseCountInGroups = zeros(groupCount, 1);
[~, uniqueCaseNameIndex] = unique(data.caseName);
groupOfCase = whichGroup(data.caseAge(uniqueCaseNameIndex), agePartition);
for i = 1:groupCount    % inspect cases in each age group
    idx = groupOfCase == i; % index of cases in group i
    caseCountInGroups(i) = sum(idx);  % case count in group i
end
end