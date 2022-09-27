function C = estimateContactMatrix_LSE(contactDataMatrix, caseCountInGroups, agePopulation, varargin)
%%% least square estimation of contactMatrix C

switch nargin
    case 3
        plotFlag = 1; % plot the heatmap as a default
    case 4
        plotFlag = varargin{1};
    case 5
        plotFlag = varargin{1};
        ageLabels = varargin{2};
end


A = contactDataMatrix;
N = agePopulation;
groupCount = size(A,1);
C = zeros(groupCount);

for i = 1 : groupCount-1
    for j = i+1 : groupCount
       C(i,j) = ( A(i,j) + N(i)/N(j) * A(j,i) ) / (1 + N(i)^2/N(j)^2);
       C(j,i) = N(i)/N(j) * C(i,j);
    end
end

for i = 1 : groupCount
    C(i,i) = A(i,i);
end


matrixplot(C, 'contactMatrix for age groups, estimated by LS', ageLabels);


end