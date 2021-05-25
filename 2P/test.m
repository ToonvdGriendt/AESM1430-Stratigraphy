clear all; clc; close all;

type = 'Spherical';         % variogram type (Gaussian, Exponential, Spherical)
dim = [40 60 1];            % reservoir dimensions 
param = [10 10 0 1 0.1];    %[SigmaX,SigmaY,Teta,Variance,RNugget]
well_pos = [10 10 ; 40 40]; % well locations
y_cond = [1 ; 1 ];          % well data
nReal = 1;                  


% [H,h,dim,XYCovariogram] = GeoStat_Covariance2Filter(type,dim,param);

real = createConditionedReal(type,dim,param,y_cond, well_pos,nReal);

for r=1:nReal;
    for i=1:dim(1)+1
        for j=1:dim(1)+1
            realz(i,j,r) = real((i-1)*(dim(1)+1)+j,r);
        end
    end
end
