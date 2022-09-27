function dVC = feasibleDirection(i,j)

dVC = zeros(8,4);
dVC(i,j) = -1;
dVC(i,j+1) = 1;
end