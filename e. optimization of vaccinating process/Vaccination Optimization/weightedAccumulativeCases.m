function weightedOutput = weightedAccumulativeCases(fun, params, tInit, xInit, tFinal, stepSize, varargin)
%%% Compute all accumulative number of cases 
%%% INPUT:  
%   fun         : the derivative of the systems of ODEs
%   params      : parameters in ODEs, a struct
%   tInit       : time instance for simulation start
%   xInit       : initial value for the systems of ODEs
%   tFinal      : time instance for simulation end
%   stepSize    : step size for ode solvers
%   weight      : (optional) 8*4 weight matrix for accumulative cases (to obtain the hospital count and mobility count)
%
%%% OUTPUT:
% accumulativeNumberOfCases

% optional weight (ones as default)
if nargin > 6
    weight = varargin{1};
else 
    weight = ones(8,4);
end


% solve
[x,t] = odeSolveRK4(fun, tInit, xInit, tFinal, stepSize);
[V, E, F, I, A, Q, R] = extractVariableFromMatrix_forVEFIAQR(x);

% daily New Cases
p = params.p;
omegap = params.omegap;
omega = params.omega;
dailyNewCases = p(:)'.*omega(:)'.*E + (1-p(:)').*omegap(:)'.*E;

% accumulative via integral
accumulativeNumberOfCases = trapz(t, dailyNewCases);

% output
weight = weight(:);
weightedOutput = accumulativeNumberOfCases * weight;

end