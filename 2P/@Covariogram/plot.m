function plot(covariogram,XMax,YMax,DeltaX,DeltaY)
% Covariogram/plot displays Covariogram object
% plot(covariogram,XMax,YMax,DeltaX,DeltaY)

if nargin == 1
    [xSize, ySize] = size(covariogram)
    XMax = 1.5 * xSize
    YMax = 1.5 * ySize
    DeltaX = 2.0 * XMax / 50
    DeltaY = 2.0 * YMax / 50
end

X = (-XMax:DeltaX:XMax);
X = ones(length(X),1) * X;

Y = (-YMax:DeltaY:YMax);
Y = Y' * ones(1,length(Y));

fprintf('size(X) is (%d,%d)\n',size(X))
fprintf('size(Y) is (%d,%d)\n',size(Y))

XYCovariogram = covariance(covariogram,X,Y);

fprintf('size(XYCovariogram) is (%d,%d)\n',size(XYCovariogram))


% mesh((-XMax:DeltaX:XMax),(-YMax:DeltaY:YMax),XYCovariogram)
surf((-XMax:DeltaX:XMax),(-YMax:DeltaY:YMax),XYCovariogram)
% shading interp
grid on
box on

