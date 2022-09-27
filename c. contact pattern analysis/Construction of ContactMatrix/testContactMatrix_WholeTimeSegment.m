clear; close all; clc;

%% read the processed case-contact data table
opts = detectImportOptions('contactDataProcessed.xlsx');
opts.VariableTypes(1,[7,8]) = repmat({'datetime'},1,2);
caseContactData = readtable('contactDataProcessed.xlsx',opts);


%  extract data with (caseAge, contactAge) available 
caseContactData = caseContactData(~isnan(caseContactData.caseAge) & ~isnan(caseContactData.contactAge), :); 
vaildCaseContactAgePairsCount = size(caseContactData, 1)  
[uniqueCaseName, idx] = unique(caseContactData.caseName);
uniqueCaseNameCount = numel(uniqueCaseName) 



% ImportAgePartition;
agePartition = [0,10,20,30,40,50,60,70,200];
ageLabels = {'0 to 9', '10 to 19','20 to 29', '30 to 39', '40 to 49', '50 to 59', '60 to 69', '\geq 70'};

% for age groups in Xiamen
%agePartition = [0, 3, 10, 20, 30, 40, 50, 60, 70, 80, 200]; 
% ageLabels = {'0 to 2', '3 to 9', '10 to 19','20 to 29', '30 to 39', '40 to 49', '50 to 59', '60 to 69', '70 to 79', '\geq 80'};

% caseCount in age groups
caseCountInGroups = findCaseCountInGroups(caseContactData(idx,:), agePartition)
% caseCountInGroups = zeros(8, 1);
% groupOfCase = whichGroup(caseContactData.caseAge(idx), agePartition);
% for i = 1:8    % inspect cases in each age group
%     idx = groupOfCase == i; % index of cases in group i
%     caseCountInGroups(i) = sum(idx);  % case count in group i
% end



%% Construct A, the Contact Data Matrix
% contactDataMatrix(i,j) denotes the daily average number of close contact in group
% j, will produced by a individual in group i. 


duration = 4 * ones(size(caseContactData, 1), 1);  % case-wise duration of making contacts set as 4-day constant according to data

contactDataMatrix = computeContactDataMatrix(caseContactData, duration, caseContactData.caseAge, caseContactData.contactAge, agePartition, caseCountInGroups);
contactDataMatrix(isnan(contactDataMatrix)) = 0;

matrixplot(contactDataMatrix, 'contactDataMatrix', ageLabels);  % heatmap of contactDataMatrix
exportgraphics(gca,'contactDataMatrix.jpg','Resolution',300); % export figure



%% read population Data
populationData = readtable('populationData.xlsx');
% head(populationData)
ageGroup = populationData.ageGroup;
population = populationData.all;

agePopulation = zeros(numel(agePartition)-1,1);
for i = 1:numel(agePartition)-1
    agePopulation(i) = sum(population(ageGroup>=agePartition(i) & ageGroup < agePartition(i+1)),'omitnan');
end
writematrix(agePopulation,'agePopulationVector.xlsx');

plotFlag = 1;
%% least square estimation of contactMatrix C
C_LSE = estimateContactMatrix_LSE(contactDataMatrix, caseCountInGroups, agePopulation, plotFlag, ageLabels);
exportgraphics(gca,'contactMatrix_LSE.jpg','Resolution',300);

%% weighted least square estimation of contactMatrix C
C_WLSE = estimateContactMatrix_WLSE(contactDataMatrix, caseCountInGroups, agePopulation, plotFlag, ageLabels);
exportgraphics(gca,'contactMatrix_WLSE.jpg','Resolution',300);

%% maximum likelihood estimation of contact matrix C
C_MLE = estimateContactMatrix_MLE(contactDataMatrix, caseCountInGroups, agePopulation, plotFlag, ageLabels);
exportgraphics(gca,'contactMatrix_MLE.jpg','Resolution',300);

%% entrywise 95% confidential intervals (based MLE)
% lower bound
C_lower = icdf('Poisson',0.025,C_MLE);
matrixplot(C_lower, 'the lower bounds of the entrywise 95% CI of the contactMatrix', ageLabels); % plot
exportgraphics(gca, 'contactMatrix_MLE_lowerCI.jpg','Resolution',300); % export figure

% upper bound 
C_upper = icdf('Poisson',0.975,C_MLE);
matrixplot(C_upper, 'the upper bounds of the entrywise 95% CI of the contactMatrix', ageLabels); % plot
exportgraphics(gca, 'contactMatrix_MLE_upperCI.jpg','Resolution',300);  % export figure


%% write MLE C, C_lower, C_upper
%cd('H:\xmuph\0 HunanCovid-19 paper\Codes\step3_SEIAR_and_VEFIAQR_Simulation');
writematrix(contactDataMatrix, 'contactDataMatrix.xlsx');
writematrix(C_LSE, 'C_LSE.xlsx');
writematrix(C_WLSE, 'C_WLSE.xlsx');
writematrix(C_MLE,'C_MLE.xlsx');
writematrix(C_lower,'contactMatrixLowerBound.xlsx');
writematrix(C_upper,'contactMatrixUpperBound.xlsx');
%cd('H:\xmuph\0 HunanCovid-19 paper\Codes\step2_Construction of ContactMatrix');

lambdaMax.A = max(eig(contactDataMatrix));
lambdaMax.C_LSE = max(eig(C_LSE));
lambdaMax.C_WLE = max(eig(C_WLSE));
lambdaMax.C_MLE = max(eig(C_MLE));
lambdaMax.C_lower = max(eig(C_lower));
lambdaMax.C_upper = max(eig(C_upper));
lambdaMax
