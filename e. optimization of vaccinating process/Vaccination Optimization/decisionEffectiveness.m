function [df, fInit, dm] = decisionEffectiveness(i, j, dm, fun, params, tInit, xInit, fInit, tFinal, stepSize, varargin)
%%% evaluate the decision effectiveness of assign this dose to the ij-th
%%% group (those finished (j-1) doses in age group i)
%
%%% INPUT:
%   i,j         : the ij-th group
%   dm          : number of doses for one-time distributing
%   fun         : objective loss funtion
%   params      : model parameters
%   tInit       : time instance for simulation start
%   xInit       : current vaccine covorage at tInit
%   fInit       : the value of the objective function at current xInit
%   tFinal      : time instance for simulation end
%   stepSize    : stepSize for RK-4 solver
%
%%% OUTPUT:
%   df          : change of the objective function (the effectness of this dose)
%   fInit       : value of the objective function after update
%   dm          : output dm unchanged

% optional weight (ones as default)
if nargin > 10
    weight = varargin{1};
end

% if population size of group ij equals to 0
xInit = reshape(xInit, 8, []);
if xInit(i,j) <= 0
    df = 0;
    dm = 0;
    return;
end

% if population size is greater 0 but lesser than xInit(i,j)
if xInit(i,j) - dm < 0
    dm = xInit(i,j);
end

xInit2 = xInit;
xInit2(:,1:4) = xInit2(:,1:4) + dm * feasibleDirection(i,j);
f1 = fInit; %f1 = AccumulativeCases(fun, params, tInit, xInit1, tFinal, stepSize);
f2 = weightedAccumulativeCases(fun, params, tInit, xInit2(:), tFinal, stepSize, weight); 

df = (f2 - f1) / dm; % accumulative cases prevention of per doses
fInit = f2;

end