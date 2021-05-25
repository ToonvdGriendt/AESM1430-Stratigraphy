%% 2P

 clear all; clc; close all;
% tic


variables

for iter=1:50
runName_iter = ['model_' num2str(iter)]

real = createConditionedReal(type,tdim,param,y_cond, well_pos,tnReal);


for r=1:tnReal;
    for i=1:tdim(1)+1
        for j=1:tdim(1)+1
            realz(i,j,r) = real((i-1)*(tdim(1)+1)+j,r);
        end
    end
end
R = realz(1:50,1:50,:);
R = 0.5+(R - min(min(min(R)))) / (max(max(max(R))) - min(min(min(R))));


if tplotV == 1 && iter == 50
   figure(iter)
   imagesc(realz)
   axis xy equal tight   
end

save(runName_iter)
end

