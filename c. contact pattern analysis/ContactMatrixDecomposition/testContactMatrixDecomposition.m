clear all; close all; clc;

% read data
ImportAgePartition;  
ImportFigureLegends;
groupCount = numel(agePartition)-1;
caseContactData = readtable('contactDataProcessed.xlsx');
caseContactData = caseContactData(~isnan(caseContactData.caseAge) & ~isnan(caseContactData.contactAge), :); %  extract data with (caseAge, contactAge) available
vaildCaseContactAgePairsCount = size(caseContactData, 1)  
[uniqueCaseName, idx] = unique(caseContactData.caseName);
uniqueCaseNameCount = numel(uniqueCaseName) 

% key variables
relation = caseContactData.relation;
setting = caseContactData.setting;

uniqueSettingCN = unique(setting); % uniqueSetting in Chinese
uniqueSettingEN = {'Empty', 'Public Transports', 'Amusement', 'School', 'Household', 'Hospital',...
    'Working Area', 'Service Industry', 'Uncategorized', 'Community', 'Catering industry'};  % uniqueSetting in English

% contacts count in each settings
contactCountInSettings = zeros(numel(uniqueSettingCN),1);
for i = 1:numel(uniqueSettingCN)
    contactCountInSettings(i) = sum( strcmp(setting, uniqueSettingCN{i}) );
end
contactCountInSettings
contactCountInNonEmptySettings = sum(contactCountInSettings(2:end)) % 


% caseCount in age groups
caseCountInGroupsAll = findCaseCountInGroups(caseContactData(idx,:), agePartition)
caseCountInGroups = zeros(size(caseCountInGroupsAll, 1), numel(uniqueSettingEN)); % save caseCount in each setting

figureCount = 1;
lambdaMax = zeros(numel(uniqueSettingEN), 2);
C = zeros(11, groupCount, groupCount);
A = zeros(11, groupCount, groupCount);
for i = 1:numel(uniqueSettingEN)
    
    
    % extract rows corresponding to uniqueSetting{i}
    indices = strcmp(setting, uniqueSettingCN{i});
    caseContactDatai = caseContactData(indices,:);
    %N2(i) = numel(unique(caseContactDatai.caseName));
    
    % for caseContactData corresponding to uniqueSetting{i}, estimate the contact matrix via MLE
    [Ci, contactDataMatrix, caseCountInGroupsi] = estimateCi(caseContactDatai, agePartition);
    caseCountInGroups(:,i) = caseCountInGroupsi; % save caseCount in this setting

    % replace nan by 0
    Ci(isnan(Ci)) = 0;
    contactDataMatrix(isnan(contactDataMatrix)) = 0;
   
    
    % heatmap of contactDataMatrix
    plotTitle = ['(', num2str(figureCount), '-', '1) ', uniqueSettingEN{i}, ', Contact Data Matrix'];
    matrixplot(contactDataMatrix, plotTitle, ageLegend);  % heatmap of contactDataMatrix
    exportgraphics(gca, ['contactDataMatrix', num2str(i), '.jpg'],'Resolution',300);  % export figure

    
    % heatmap of contactMatrix
    plotTitle = ['(', num2str(figureCount), '-', '2) ', uniqueSettingEN{i}, ', Contact Matrix Estimated via MLE'];
    matrixplot(Ci, plotTitle, ageLegend);  % heatmap of contactDataMatrix
    exportgraphics(gca, ['contactMatrixByMLE', num2str(i), '.jpg'],'Resolution',300);  % export figure

    figureCount = figureCount + 1;
    
    C(i,:,:) = Ci;
    A(i,:,:) = contactDataMatrix;
    
    writematrix(Ci, ['contactMatrix_', uniqueSettingEN{i}, '.xlsx']);
    writematrix(contactDataMatrix, ['contactDataMatrix_', uniqueSettingEN{i}, '.xlsx']);

    lambdaMax(i,1) = max(eig(contactDataMatrix));
    lambdaMax(i,2) = max(eig(Ci));
    
end

lambdaTable = table;
lambdaTable.SettingNames = uniqueSettingEN';
lambdaTable.lambdaMaxForContactDataMatrices = lambdaMax(:,1);
lambdaTable.lambdaMaxForContactMatrices = lambdaMax(:,2);
writetable(lambdaTable, 'lambdaMaxForContactMatricesInDifferentSettings.xlsx');

% %% for recovering
% recoveredC = zeros(8);
% recoveredA = zeros(8);
% for i = 1:11
%     recoveredC = recoveredC + squeeze(C(i,:,:)) .* caseCountInGroups(:,i);
%     recoveredA = recoveredA + squeeze(A(i,:,:)) .* caseCountInGroups(:,i);
% end
% 
% recoveredC = recoveredC ./ sum(caseCountInGroups, 2) .* caseCountInGroupsAll;
% recoveredA = recoveredA ./ sum(caseCountInGroups, 2) .* caseCountInGroupsAll;
% 
% 
% % temp = repmat(caseCountInGroups.', [1, 1, 8]);
% % recoveredC = squeeze(sum(C , 1)) ;
% figure; heatmap(recoveredC); title('Recovered Contact Matrix');
% exportgraphics(gca, 'recoveredContactMatrix.jpg','Resolution',300);  % export figure
% % 
% % recoveredA = squeeze(sum(A , 1));
% figure; heatmap(recoveredA); title('Recovered Contact Data Matrix');
% exportgraphics(gca, 'recoveredContactDataMatrix.jpg','Resolution',300);  % export figure


% %% try to represent Empty and Uncategorized as linear combination of other contact matrices of 9 settings, initialize for solve linear equation
% vecC = zeros(groupCount^2, 9); 
% vecA = zeros(groupCount^2, 9); 
% count = 1;
% for i = 1:numel(uniqueSettingEN)    % total 11 categories
%     if i ~= 1 && i ~= 9  % non-Empty and non-Uncategorized
%         vecC(:, count) = squeeze(C(i,:));
%         vecA(:, count) = squeeze(A(i,:));
%         count = count + 1;
%     end
% end
% 
% 
% C1 = squeeze(C(1,:,:));  vecC1 = C1(:)
% C9 = squeeze(C(9,:,:));  vecC9 = C9(:)
% A1 = squeeze(A(1,:,:));  vecA1 = A1(:)
% A9 = squeeze(A(9,:,:));  vecA9 = A9(:)
% 
% 
% %%% optimizing the non-negative coefficients
% 
% % non-constrained least square, which lead to negative coefficient
% coeffEmptyC = vecC \ vecC1
% coeffEmptyA = vecA \ vecA1
% coeffUncategorizedC = vecC \ vecC9
% coeffUncategorizedA = vecA \ vecA9
% 
% % constrained 
% objfunC1 = @(p) norm(vecC9 - vecC * p);
% p0 = ones(9,1); % iniital value
% lb = zeros(9,1);
% [p,fval,exitflag,output,lambda,grad,hessian] = fmincon(objfunC1, p0, [], [], [], [], lb, []);
% res = reshape(vecC1 - vecC * p, [8,8])
% fittedCoefficients = p
% 
% figure; 
% tiledlayout('flow');
% nexttile;
% heatmap(C1); 
% nexttile;
% heatmap(res)
% 
