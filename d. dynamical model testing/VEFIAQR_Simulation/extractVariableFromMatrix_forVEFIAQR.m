function [V, E, F, I, A, Q, R] = extractVariableFromMatrix_forVEFIAQR(x)
% the dimension of temp told everything
tempID = 1:32;
V = x(:,tempID + 0*32);
E = x(:,tempID + 1*32);
F = x(:,tempID + 2*32);
I = x(:,tempID + 3*32);
A = x(:,tempID + 4*32);
Q = x(:,tempID + 5*32);
R = x(:,tempID + 6*32);
end