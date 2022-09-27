function dxdt = dxdt_VEFIAQR(x,t,params)
%%% Consider the differencial equation: dx/dt = f(x,t,params), where x denotes a vector of
%%% m dimension. (Let n denotes the number of subgroups, and m = 5*n)
%%% This function evaluates the derivative f(x,t) at a given point (x,t)
% INPUT:
%  x = [V0, V1, V2, V3,...
%       E0, E1, E2, E3,...
%       F0, F1, F2, F3,...
%       I0, I1, I2, I3,...
%       A0, A1, A2, A3,...
%       Q0, Q1, Q2, Q3,...
%       R0, R1, R2, R3];   a n*28 matrix

% % %  x = [[V0, E0, F0, I0, A0, Q0, R0],...
% % %       [V1, E1, F1, I1, A1, Q1, R1],...
% % %       [V2, E2, F2, I2, A2, Q2, R2], ...
% % %       [V3, E3, F3, I3, A3, Q3, R3]];      
%
% 
%    n   : number of age groups
%    t   : time variable
%  params: a struct of model parameters (assume that all parameters are group-distincted)
%
%       params.Beta   : n*n matrix, whose entries are beta_ij (beta_ij denotes the transimission capacity from group i to group j)
%       params.VE     : vaccine effecacy
%       params.kappa  : describe the difference of Transmission capacity between Iij and Aij
%       params.omega  : transition rate from Eij to Fij
%       params.omegap : (omegap is the short of omega_prime) transition rate from Fij to Iij
%       params.omegapp: transition rate from Fij to Aij
%       params.p      : proportion of asymptomatic cases
%       params.gamma  : the inverse of infectious period of Iij
%       params.gammap : (gammap is the short of gamma_prime) the inverse of average infectious period of A
%       params.gammapp: remove rate of Qij
%       params.mu     : rate of quarantine
%       params.f      : fatality rate of Iij
%       params.br     : birth rate of group ij (INTERSECTION of age group i and vaccinated group j)
%       params.dr     : natrual death rate of group ij 
%
%  NOTE:
%      n*n matrices: Beta
%      n*4 matrices: kappa, omega, omegap, p, gamma, gammap, f, br, dr
%          vectors : Vi, Ei, Fi, Ii, Ai, Qi, Ri, i = 1, 2,..., n.
%          scalars : mu
%               
%
% OUTPUT:   
%       f = f(x,t,params), a vector of derivatives
%
tempID = 1:4;
V = x(:, tempID + 4*0);
E = x(:, tempID + 4*1);
F = x(:, tempID + 4*2);
I = x(:, tempID + 4*3); 
A = x(:, tempID + 4*4);
Q = x(:, tempID + 4*5);
R = x(:, tempID + 4*6);
% tempID = [1, 8, 15, 22];
% V = x(:, tempID);
% E = x(:, tempID + 1);
% F = x(:, tempID + 2);
% I = x(:, tempID + 3); 
% A = x(:, tempID + 4);
% Q = x(:, tempID + 5);
% R = x(:, tempID + 6);
N = sum(x,2); 

Beta = params.Beta;
VE = params.VE;
kappa = params.kappa;
omega = params.omega;
omegap = params.omegap;
omegapp = params.omegapp;
p = params.p;
gamma = params.gamma;
gammap = params.gammap;
gammapp = params.gammapp;
mu = params.mu;
f = params.f;
br = params.br;
dr = params.dr;


infectiveVector = sum(I + kappa .* (E + A), 2); % \sum_k  (I_ji + kappa(E_jk + Ajk))
infectionForce = Beta.' * infectiveVector;     % lambda, a n*1 vector


dVdt = -infectionForce.*(1-VE).*V;
dEdt =  infectionForce.*(1-VE).*V - omega.*E;
dFdt = omega.*E - (1-p).*omegap.*F - p.*omegapp.*F;
dIdt = (1-p).*omegap.*F - mu.*I - gamma.*I;
dAdt = p.*omegapp.*F - mu.*A - gammap.*A;
dQdt = mu.*A + mu.*I - gammapp.*Q;
dRdt = gamma.*I + gammap.*A + gammapp.*Q;


dxdt = [dVdt, dEdt, dFdt, dIdt, dAdt, dQdt, dRdt];
end