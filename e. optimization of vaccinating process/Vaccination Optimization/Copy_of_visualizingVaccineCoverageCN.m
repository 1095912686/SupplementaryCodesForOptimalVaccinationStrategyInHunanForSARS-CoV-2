clear; close all; clc;

n = 8;
ImportFigureLegends;

caseData129 = readtable('caseData129.xlsx');
nVC = readmatrix('nVC.xlsx'); % contactCount in (ageGroups, dosesGroups)
nVCC = readmatrix('nVCC.xlsx'); % caseCount in (ageGroups, dosesGroups)
nC = readmatrix('nC.xlsx'); % caseCount in ageGroups
caseAge = readmatrix('caseAge.xlsx'); % caseAge
N = readmatrix('agePopulationVector.xlsx');

%% Figure 0: visualize VC, the vaccine Coverage rate
VC = nVC ./ sum(nVC,2) * 100;
VC = [VC; sum(nVC,1) / sum(nVC,"all") * 100];
VC=VC(:,[4,3,2]); % 3, 2, 1-doses

figure;
ba = bar(VC,'stacked','FaceColor','flat');

c = colormap("summer");
% c = [102,194,165;...
%     141, 160, 203;...
%     252,141,98] /256;
for i = 1:size(VC,2)
    ba(i).CData = c(60*i+60,:);
    %ba(i).CData = c(i,:);
end

% set(gca,'XTickLabel',[ageLegendCN,'总计']);
% legend('仅接种了3剂次的','仅接种了2剂次的','仅接种了1剂次的');
% ylabel('覆盖率');
set(gca,'XTickLabel',[ageLegend,'Total']);
%legend('3-doses only','2-doses only','1-dose only');
legend('booster vaccinated', 'fully vaccinated', 'vaccinated but un-fully')
ylabel('Coverage Rate');


%%% add number for each bins of VC
[~, ~, ytips1] = ba.YEndPoints; ytips1 = ytips1 + 2;
[~, ~, xtips1] = ba.XEndPoints; xtips1 = xtips1 - 0.3;

for i=1:size(VC,1)
    % if this bin is too short
    if sum(VC(i,:)) < 5
        labels1 = num2str(sum(VC(i,:),2), '%.1f %%');
        text(xtips1(i),ytips1(i),labels1,'VerticalAlignment','middle');
        continue;
    end
    
    % otherwise, there are places to label the stacked bins
    for j=1:size(VC,2)
        if VC(i,j)>0
            labels_stacked=num2str(sum(VC(i,1:j)),'%.1f %%');
            hText = text(i, sum(VC(i,1:j),2), labels_stacked);
            set(hText, 'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize',10, 'Color','black');
        end
    end

end

%% Figure 1: Total Attack Rate with different doses of vaccination
P = nVCC ./ nVC;
contactInDoseGroup = sum(nVC, 1);
caseInDoseGroup = sum(nVCC, 1);
TAR = caseInDoseGroup ./ contactInDoseGroup;


figure;
b2 = bar(TAR);
set(b2, 'facecolor',[1,1,1]);
% xlabel('疫苗接种情况');
% ylabel('不同群体的患病风险');
% title('不同疫苗接种条件下的罹患率');
% set(gca,'XTickLabel',{'未接种','接种1剂','接种2剂','接种3剂'});
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

XLABEL = '密接人数';
YLABEL = '频率/组距';
TITLE = '';
binWidth = 20;

fitDist(numberOfCloseContact, binWidth, XLABEL, YLABEL, TITLE);



%% Figure 3: Age distribution of cases
caseAge(isnan(caseAge)) = [];

XLABEL = '年龄';
YLABEL = '频率/组距';
TITLE = '';
binWidth = 5;

fitDist(caseAge, binWidth, XLABEL, YLABEL, TITLE);




%% Figure 4: time interval from illness onset to positively tested
D1 = caseData129.DateOfIllnessOnset;
D2 = caseData129.DateOfPositiveTest;
dD21 = days(D2 - D1);
dD21(isnan(dD21)) = [];
dD21(dD21 == 0) = [];

XLABEL = '发病到阳性诊断的时间间隔（天）';
YLABEL = '频率/组距';
TITLE = '';
binWidth = 0.5;

fitDist(dD21, binWidth, XLABEL, YLABEL, TITLE);

%% Figure 5: time interval from possible date of infection to illness onset
D0 = caseData129.PossibleDateOfInfection;
idx = ~isnan(D0);
D0 = D0(idx);
D0 = datetime(D0,'ConvertFrom','excel');

XLABEL = '可能感染日期到发病的时间间隔（天）';
YLABEL = '频率/组距';
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
