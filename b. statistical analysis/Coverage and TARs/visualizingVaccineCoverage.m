clear; close all; clc;

n = 8;
ImportFigureLegends;

caseData129 = readtable('caseData129.xlsx');
%[nVC,nVCC] = vaccineCoverageViaDataByIsComplete()
%[nVC,nVCC] = vaccineCoverageViaDataByDoseCount()
nVC = readmatrix('nVC.xlsx'); % contactCount in (ageGroups, dosesGroups)
nVCC = readmatrix('nVCC.xlsx'); % caseCount in (ageGroups, dosesGroups)
nC = readmatrix('nC.xlsx'); % caseCount in ageGroups
caseAge = readmatrix('caseAge.xlsx'); % caseAge
N = readmatrix('agePopulationVector.xlsx');

%% Figure 0: visualize VC, the vaccine Coverage rate
VC = nVC ./ sum(nVC,2) * 100; % persentage normalized by age groups
VC = [VC; sum(nVC,1) / sum(nVC,"all") * 100]; % add one row of summation of age groups
VC = VC(:,[3,2]); % extract fully, un-fully vaccinated columns for visualize

figure(1);
ba = bar(VC,'stacked','FaceColor','flat');

c = colormap("summer");
% c = [102,194,165;...
%     141, 160, 203;...
%     252,141,98] /256;
for i = 1:size(VC,2)
    ba(i).CData = c(60*i+60,:);
    %ba(i).CData = c(i,:);
end


set(gca,'XTickLabel',[ageLegend,'Total'], 'FontSize', 16);
%legend('3-doses only','2-doses only','1-dose only');
%legend('booster vaccinated', 'fully vaccinated', 'vaccinated but un-fully')
legend('fully vaccinated', 'vaccinated but un-fully', 'Location', 'northwest')
ylabel('Coverage Ratio');


%%% add number for each bins of VC
%[~, ~, ytips1] = ba.YEndPoints; ytips1 = ytips1; %  + 2;
ytips1 = ba.YEndPoints;
ytips1 = ytips1 + 2;
xtips1 = ba.XEndPoints;
xtips1 = xtips1;

for i=1:size(VC,1)
    % if this bin is too short
    if sum(VC(i,:)) < 5
        labels1 = num2str(sum(VC(i,:),2), '%.1f %%');
        text(xtips1(i)-0.3,ytips1(i),labels1,'VerticalAlignment','middle');
        continue;
    end
    
    % otherwise, there are places to label the stacked bins
    for j=1:size(VC,2)
        if VC(i,j)>0
            labels_stacked = num2str(sum(VC(i,1:j)),'%.1f %%');
            hText = text(xtips1(i), sum(VC(i,1:j),2)+4, labels_stacked);
            set(hText, 'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize',10, 'Color','black');
        end
    end

end
fig1 = setupCurrentFigure(figure(1));
exportgraphics(fig1, 'vaccineCoverage_barplot.jpg','Resolution',300);  % export figure
savefig(fig1, 'vaccineCoverage_barplot.fig');

%% Figure 1: Total Attack Rate with different doses of vaccination
P = nVCC ./ nVC;
contactInDoseGroup = sum(nVC, 1);
caseInDoseGroup = sum(nVCC, 1);
TAR = caseInDoseGroup ./ contactInDoseGroup;


figure(2)
b2 = bar(TAR);
set(b2, 'facecolor',[1,1,1]);

xlabel('Condition of Vaccination');
ylabel('Total attack rate (TAR)');
set(gca,'XTickLabel',{'none','unfully vaccinated','fully vaccinated','boster vaccinated'});

for i = 1:4
    labels_stacked=[num2str(sum(nVCC(:,i))),' cases'];
    hText = text(i, TAR(i)+3e-7, labels_stacked);
    set(hText, 'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize',10, 'Color','black');
end



%% Figure 2: Distribution of Number of Close Contact of those 129-th cases
numberOfCloseContact = caseData129.NumberOfCloseContact;
numberOfCloseContact(isnan(numberOfCloseContact)) = [];

XLABEL = 'Number of Close Contact';
YLABEL = 'density (frequency / binWidth)';
TITLE = '';
binWidth = 20;

fitDist(numberOfCloseContact, binWidth, XLABEL, YLABEL, TITLE);



%% Figure 3: Age distribution of cases
caseAge(isnan(caseAge)) = [];

XLABEL = 'Age';
YLABEL = 'density (frequency / binWidth)';
TITLE = '';
binWidth = 5;

fitDist(caseAge, binWidth, XLABEL, YLABEL, TITLE);




%% Figure 4: time interval from illness onset to positively tested
D1 = caseData129.IllnessOnset;
D2 = caseData129.PositiveTestDate;
dD21 = days(D2 - D1);
dD21(isnan(dD21)) = [];
dD21(dD21 == 0) = [];

XLABEL = 'Time Interval from Illness Onset to Diagnosed (in days)';
YLABEL = 'density (frequency / binWidth)';
TITLE = '';
binWidth = 0.5;

fitDist(dD21, binWidth, XLABEL, YLABEL, TITLE);

%% Figure 5: time interval from possible date of infection to illness onset
D0 = caseData129.ExposedDate;
idx = ~isnat(D0);
D0 = D0(idx);
%D0 = datetime(D0,'ConvertFrom','excel');

XLABEL = 'Possible Date of Infection to Illness Onset';
YLABEL = 'density (frequency / binWidth)';
TITLE = '';
binWidth = 1;
dD10 = days(D1(idx) - D0);
dD10(dD10==0) = [];

fitDist(dD10, binWidth, XLABEL, YLABEL, TITLE);

%% Figure 6: time interval from possible date of infection to illness onset
%% GT distribution
data = readtable('GTdata.xlsx');
gt = data.generation_time;
maxGT = max(gt);
minGT = min(gt);
len = maxGT-minGT;
freqGT = zeros(len,1);
for i = 1:len+1
    freqGT(i) = sum(gt == minGT + i - 1);
end
freqGT(freqGT==0) = [];

XLABEL = 'Generation Time (in days)';
YLABEL = 'frequency';
TITLE = '';
binWidth = 1;
fitDist(gt,binWidth, XLABEL, YLABEL, TITLE);
