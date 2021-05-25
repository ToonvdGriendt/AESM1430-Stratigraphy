function [xSize, ySize] = size(covariogram)
% Covariogram/returns typical size of covariogram
% [xSize, ySize] = size(covariogram)


xSize = sqrt(covariogram.CDist(1,1));
ySize = sqrt(covariogram.CDist(2,2));