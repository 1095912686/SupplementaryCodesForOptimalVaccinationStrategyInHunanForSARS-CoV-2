function group = whichGroup(age, agePartition)
% given a partition like: agePartition = [0,10,20,30,40,50,60,70,200];
% input a scalar age
% ouput the group that such age belongs to.

groupCount = numel(agePartition) - 1;

group = zeros(size(age));

for k = 1:groupCount
    idx = age >= agePartition(k) & age < agePartition(k+1);
    group(idx) = k;
end

end

