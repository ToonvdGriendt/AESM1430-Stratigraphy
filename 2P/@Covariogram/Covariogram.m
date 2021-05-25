function covariogram = Covariogram(Sigma1,Sigma2,Teta,Variance,RNugget,VarioType)
% Covariogram/Covariogram Covariogram class constructor
% covariogram = Covariogram(SigmaX,SigmaY,Teta,Variance,RNugget,VarioType)
%
% members of Covariogram object are
% covariogram.Sigma1
% covariogram.Sigma2
% covariogram.Teta
% covariogram.CDist
% covariogram.CDistInv
% covariogram.Variance
% covariogram.RNugget
% covariogram.VarioType

covariogram.Sigma1 = Sigma1;
covariogram.Sigma2 = Sigma2;
covariogram.Teta   = Teta;

CosTeta = cos(Teta * pi / 180);
SinTeta = sin(Teta * pi / 180);


C = [ [ (CosTeta^2 * Sigma1^2 + SinTeta^2 * Sigma2^2) ...
        ((CosTeta * SinTeta) * (Sigma1^2 - Sigma2^2)) ]; ...
      [ ((CosTeta * SinTeta) * (Sigma1^2 - Sigma2^2)) ...
        (SinTeta^2 * Sigma1^2 + CosTeta^2 * Sigma2^2) ] ];

covariogram.CDist    = C;    
covariogram.CDistInv = inv(C);
covariogram.Variance = Variance;
covariogram.RNugget = RNugget;
covariogram.VarioType = VarioType;

covariogram = class(covariogram,'Covariogram');
