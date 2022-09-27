function C = estimateContactMatrix_MLE(contactDataMatrix, caseCountInGroups, agePopulation, varargin)
%%% maximum likelihood estimation of contact matrix C

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
K = caseCountInGroups;
groupCount = size(A,1);

C = ((K).*(A) + (K.') .* (A.')) ./ (K + N ./ (N.') .* (K.'));



if plotFlag == 1
    matrixplot(C, 'contactMatrix for age groups, estimated by MLE', ageLabels);
end

end