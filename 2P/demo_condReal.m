type    = 'Spherical';      % variogrma type: Spherical, Exponential or Gausian
gridDim = [101 71 1];       %dimensions of grid
param   = [50,20,30,.10,.011];   %variogram parameters: [SigmaX,SigmaY,Teta,Variance,RNugget] 
y_cond  = .5*ones(4,1);    %values in wells
wellPos = [15 15;30 30;20 50;90,2]; %position of wells
nReal   = 10;               %number of conditioned realizations

realizations = createConditionedReal(type,gridDim,param,y_cond,wellPos,nReal);
max(max(realizations))