function display(covariogram)
% Covariogram/display displays Covariogram object

fprintf('Mebers of Covariogram object are:\n');
fprintf('Sigma1     is %12.6g\n',covariogram.Sigma1);
fprintf('Sigma2     is %12.6g\n',covariogram.Sigma2);
fprintf('Teta       is %12.6g\n',covariogram.Teta);
fprintf('Variance   is %12.6g\n',covariogram.Variance);
fprintf('RNugget    is %12.6g\n',covariogram.RNugget);
fprintf('VarioType  is %s\n',covariogram.VarioType);
fprintf('\n');
fprintf('Typical size of Covariogram object is:\n');
fprintf('LX         is %12.6g\n',sqrt(covariogram.CDist(1,1)));
fprintf('LY         is %12.6g\n',sqrt(covariogram.CDist(2,2)));