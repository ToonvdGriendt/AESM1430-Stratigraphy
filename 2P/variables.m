%% variables


runName = 'testRun';        % name that data will be saved


type = 'Spherical';         % variogram type (Gaussian, Exponential, Spherical)
tdim = [50 50 1];           % reservoir dimensions 
param = [10 20 0 1 0.1];    %[SigmaX,SigmaY,Theta,Variance,RNugget]
well_pos = [10 10 ; 40 40]; % well locations
y_cond = [1 ; 1 ];          % well data
tnReal = 1;                 % number of realizations to make

tplotV = 1;                 % 0: do not plot, 1: plot




