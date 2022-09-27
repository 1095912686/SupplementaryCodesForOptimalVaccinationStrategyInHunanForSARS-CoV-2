function [S,E,I,A,R] = extractVariableFromCells_forSEIAR(x)



m = numel(x);
[n,~] = size(x{1});

% save S
if nargout >= 1
    S = zeros(m,n);
    for i = 1:m
        S(i,:) = (x{i}(:,1)).';
    end
    varargout{1} = S;
end

% save E
if nargout >= 2
    E = zeros(m,n);
    for i = 1:m
        E(i,:) = (x{i}(:,2)).';
    end
    varargout{2} = E;
end

% save I
if nargout >= 3
    I = zeros(m,n);
    for i = 1:m
        I(i,:) = (x{i}(:,3)).';
    end
    varargout{3} = I;
end


% save A
if nargout >= 4
    A = zeros(m,n);
    for i = 1:m
        A(i,:) = (x{i}(:,4)).';
    end
    varargout{4} = A;
end


% save R
if nargout >= 5
    R = zeros(m,n);
    for i = 1:m
        R(i,:) = (x{i}(:,5)).';
    end
    varargout{5} = R;
end