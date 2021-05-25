%% 2P

clear all; clc; close all;
tic


variables

real = createConditionedReal(type,dim,param,y_cond, well_pos,nReal);

for r=1:nReal;
    for i=1:dim(1)+1
        for j=1:dim(1)+1
            realz(i,j,r) = real((i-1)*(dim(1)+1)+j,r);
        end
    end
end

