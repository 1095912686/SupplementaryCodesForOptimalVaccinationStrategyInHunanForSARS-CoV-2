clear; close all; clc;

ImportFigureLegends; % import the following cell array for figure legend: groupLegend; ageLegend;

% read data
%% read data of CaseAndContacts
opts2 = detectImportOptions('contactDataProcessed.xlsx');
[opts2.VariableTypes{[7:8, 28, 31, 34]}] = deal('datetime', 'datetime', 'datetime', 'datetime', 'datetime'); % illnessOnset, diagnosedDate, time1,2,3
opts2.VariableTypes{15} = 'double'; % contactAge
caseContactData = readtable('contactDataProcessed.xlsx', opts2);

% denote exposure ranking groups by digits
%exposureRanking = cellfun(@(A)str2double(A(isstrprop(A,'digit'))), caseContactData.exposureRanking, 'UniformOutput', false);
exposureRanking = caseContactData.exposureRanking;
%caseContactData.exposureRanking = cell2mat(exposureRanking);

opts3 = detectImportOptions('caseData129.xlsx');
caseData129 = readtable('caseData129.xlsx',opts3);

%% extract info of a specific case
allCaseName = unique(caseData129.CaseName);

%% isCase as the outcome
isCase = cellfun(@(A)sum(strcmp(A, allCaseName)), caseContactData.contactName);
sum(isCase)     % 385 cases in close contacts due to repeat contacts

%% variables that may impact the outcome of isCase
% vaccineStatus
vaccines = caseContactData(:, 27:end);
vaccinationStatus = zeros(size(caseContactData, 1), 1);
doseNum = ~isnat(vaccines.time1) + ~isnat(vaccines.time2) + ~isnat(vaccines.time3);
doseDemand = max([vaccines.doseDemand1, vaccines.doseDemand2, vaccines.doseDemand3], [], 2);
for i = 1:size(caseContactData,1)

    if doseNum(i) == 0
        vaccinationStatus(i) = 0;    % none
    elseif doseNum(i) < doseDemand(i) 
        vaccinationStatus(i) = 1;    % unfully vaccinated
    elseif doseNum(i) == doseDemand(i)
        vaccinationStatus(i) = 2;    % fully vaccinated
    elseif doseNum(i) > doseDemand(i)
        vaccinationStatus(i) = 3;    % booster vaccinated
    elseif isnan(doseDemand(i))
        vaccinationStatus(i) = doseNum(i);  % if doseDemand missing, set vaccinationStatus to be doses finished as defalut
    end 

end

% contact gender
contactGender = zeros(size(caseContactData,1),1);
maleIndex = strcmp(caseContactData.contactGenderUpdate, '男');
contactGender(maleIndex) = 1;
contactGender(~maleIndex) = 0;

% other possible impact factors (contactAge, exposureRanking)
dataForTraining = caseContactData(:, [15, 19]); % contactAge, exposureRanking
dataForTraining.contactGender = contactGender;
dataForTraining.vaccinationStatus = vaccinationStatus;
dataForTraining.isCase = isCase;
dataForTraining(isnan(dataForTraining.contactAge), :) = [];

% resampling for sample size balance
positiveIndex = find(dataForTraining.isCase);
resampledIndex = repmat(positiveIndex, [round(size(dataForTraining,1) / sum(dataForTraining.isCase))-1, 1]);
dataForTrainingBalanced = [dataForTraining; dataForTraining(resampledIndex,:)];

%%
% based on the resampled data set: dataForTrainingBalanced, we use the
% Regression Learner toolbox in MATLAB to train the logistic model.
%%
% Let Y denotes a 0-1 random variable with Y=1 represent the close contact
% is infected, and Y=0 for not infected.
modelspec = 'logit(isCase) ~ vaccinationStatus + contactAge + contactGender + exposureRanking';
model_logistic = fitglm(dataForTrainingBalanced, modelspec, 'Distribution','binomial')

% modelspec = 'logit(isCase) ~ vaccinationStatus + contactAge + contactGender + exposureRanking';
% model_logistic = stepwiseglm(dataForTrainingBalanced, modelspec, 'Distribution','binomial','upper','linear')

% modelspec = 'logit(isCase) ~ vaccinationStatus + contactAge + contactGender';
% model_logistic = stepwiseglm(dataForTrainingBalanced, modelspec, 'Distribution','binomial', 'upper', 'linear')

% extract coefficients of vaccinationStatus and contactAge
coeff1 = model_logistic.Coefficients{'vaccinationStatus', 'Estimate'};  % coefficient of vaccinatation status
coeff2 = model_logistic.Coefficients{'contactAge', 'Estimate'}; % coefficient of age

% group-specific OR and VE
N = readmatrix('agePopulationVector.xlsx'); % the population vector (stratified by age)
nVC = readmatrix('nVC.xlsx'); % vaccineConvergenceViaData();

averageAge = (5:10:75) * N / sum(N)
averageVaccination = (0:3) * sum(nVC, 1)' / sum(nVC, 'all')

OR1 = exp(((0:3) * coeff1 + (5:10:75).' * coeff2) - (0*coeff1 + 75 * coeff2))  % versus (un-vaccinated, > 70 years old)
OR2 = exp(((0:3) * coeff1 + (5:10:75).' * coeff2) - (averageVaccination*coeff1 + averageAge * coeff2))  % versus (average vaccination, average age)
OR3 = exp(((0:3) * coeff1 - (0*coeff1)))  % Odds ratio
OR4 = exp(((5:10:75)' * coeff2 - (averageAge*coeff2)))
VE.pointEstimation = 1 - OR3
sigma.pointEstimation = 1 - OR4

writematrix(OR1, 'OR1.xlsx');
writematrix(OR2, 'OR2.xlsx');
%writematrix(VE, 'VE.xlsx');


%% pair bootstrap, resampling on contacts
%plotSlice(model_logistic)
%plotDiagnostics(model_logistic)
%plotResiduals(model_logistic)
[ynew,ynewci] = predict(model_logistic, dataForTrainingBalanced(:,1:end-1));
err = (ynew - dataForTrainingBalanced.isCase);

batchSize = 1e4; % size(dataForTrainingBalanced, 1); 
repeat = 1e3;

coefficientsRecord = zeros(repeat, 4); % []
for i = 1:repeat
    idx = randsample(size(dataForTrainingBalanced,1), batchSize, 'true'); % with replacement
    bootstrapDataSet = dataForTrainingBalanced(idx,:); % pair bootstrap
    
%     % error bootstrap
%     bootstrapDataSet.isCase = bootstrapDataSet.isCase + err(idx);
%     bootstrapDataSet.isCase = min(max(bootstrapDataSet.isCase + err(idx), 0), 1);


    model_logistic = fitglm(bootstrapDataSet, modelspec, 'Distribution','binomial');
    coefficientsRecord(i,1) = model_logistic.Coefficients{'vaccinationStatus', 'Estimate'};  % coefficient of vaccinatationStatus
    coefficientsRecord(i,2) = model_logistic.Coefficients{'contactAge', 'Estimate'}; % coefficient of contactAge
    coefficientsRecord(i,3) = model_logistic.Coefficients{'exposureRanking', 'Estimate'}; % coefficient of contactAge
    coefficientsRecord(i,4) = model_logistic.Coefficients{'contactGender', 'Estimate'}; % coefficient of contactAge
end

%% OR and VE versus (un-vaccinated &&  > 70 years old)
% resampled OR
resampledOR = arrayfun(@(coeff1, coeff2) exp(((0:3) * coeff1 + (5:10:75).' * coeff2) - (0*coeff1 + 75 * coeff2)).',...
                       coefficientsRecord(:,1), coefficientsRecord(:, 2), 'UniformOutput', 0);  % versus (un-vaccinated, > 70 years old)
resampledOR = cell2mat(resampledOR);   
resampledOR = reshape(resampledOR', 8, 4, []);

meanOR = mean(resampledOR, 3)
stdOR = std(resampledOR, [], 3)
writematrix(stdOR, 'stdOR.xlsx');
writematrix(meanOR, 'meanOR.xlsx');

% resampled VE
resampledVE = 1 - cell2mat(arrayfun(@(coeff1) exp(((0:3) * coeff1) - (0*coeff1)),...
                       coefficientsRecord(:,1), 'UniformOutput', 0));  % versus (un-vaccinated)
VE.resampledMean = mean(resampledVE);
VE.resampledStd = std(resampledVE);
VE


%% visualize
% boxchart for coefficients distribution
fig = figure;
tiledlayout(1,4,"TileSpacing","compact");
coefficientsNames = {'vaccinationStatus', 'contactAge', 'exposureRanking', 'contactGender'};
writetable(array2table(coefficientsRecord, 'VariableNames', coefficientsNames), 'coefficientsResampledRecord.xlsx');
for i = 1:4
    nexttile;
    bi = boxchart(coefficientsRecord(:,i), 'Notch', 'on');
    xticklabels([]);
    ylabel(coefficientsNames{i});

    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
    set(get(gca, 'xlabel'), 'FontSize', 14, 'FontName', 'Times New Roman');
    set(get(gca, 'ylabel'), 'FontSize', 14, 'FontName', 'Times New Roman');
    set(get(gca, 'title'), 'FontSize', 14, 'FontName', 'Times New Roman');

end
exportgraphics(fig, 'coefficientDistribution.jpg','Resolution',300);  % export figure
savefig(fig, 'coefficientDistribution.fig');

% boxchart for resampled VE
fig2 = figure;
boxplot(resampledVE, 'Notch', 'on');
ylabel('Vaccine Efficacy');
xticklabels(doseLegend);
fig2 = setupCurrentFigure(fig2);
exportgraphics(fig2, 'resampledVE.jpg','Resolution',300);  % export figure
savefig(fig, 'resampledVE.fig');


% % for boxchart of OR distribution 
% tempNames = cell(32, 1);
% for k = 1:32
%     j = ceil(k/8);
%     i = k - 8*(j-1);
%     tempNames{k} = [num2str(i), num2str(j-1)];
% end
% 
% resampledOR = arrayfun(@(coeff1, coeff2) exp(((0:3) * coeff1 + (5:10:75).' * coeff2) - (0*coeff1 + 75 * coeff2)),...
%                        coefficientsRecord(:,1), coefficientsRecord(:, 2), 'UniformOutput', 0);  % versus (un-vaccinated && > 70 years old)
% 
% 
% tableOR = cellfun(@(A)A(:)', resampledOR, 'UniformOutput', 0);
% tableOR = cell2mat(tableOR);
% tableOR = array2table(tableOR, 'VariableNames', tempNames);
% writetable(tableOR, 'tableOR.xlsx');
