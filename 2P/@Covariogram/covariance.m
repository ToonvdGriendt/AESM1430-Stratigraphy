function value = covariance(covariogram,xDist,yDist)
% Covariogram/covariance returns the value of the covariance
% value = covariance(covariogram,xDist,yDist)

weightedDistanceSquared = ...
          covariogram.CDistInv(1,1) .* xDist.^2 ...
    + 2 * covariogram.CDistInv(1,2) .* xDist .* yDist ...
        + covariogram.CDistInv(2,2) .* yDist.^2;
weightedDistance = sqrt(weightedDistanceSquared);

VarZ = covariogram.RNugget * covariogram.Variance;
VarC = (1 - covariogram.RNugget) * covariogram.Variance;
VarioType = covariogram.VarioType;

if(strcmp(VarioType,'Gaussian') == 1) 
   % disp('Gaussian VarioType has been defined')
   value = VarC .* exp(- 3 .* weightedDistanceSquared);
elseif(strcmp(VarioType,'Spherical') == 1) 
   % disp('Spherical VarioType has been defined')
   weightedDistance = min(1.,weightedDistance);
   value = VarC .* (1 - 1.5 .* weightedDistance + 0.5 .* weightedDistance.^3);
elseif(strcmp(VarioType,'Exponential') == 1) 
   % disp('Exponential VarioType has been defined')
   value = VarC .* exp(- 3 .* weightedDistance);
else
   error('Invalid VarioType has been defined');
end

nuggetList = find(weightedDistance < 0.001);

value(nuggetList) = value(nuggetList) + VarZ;