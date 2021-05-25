function Xpad = GeoStat_padzero(X,dim)

% function X = GeoStat_padzero(X,dim)
%
% Function to pad zeros around a matrix

[nx,ny]     = size(X);
Xpad        = zeros(max(nx,dim(1)),max(ny,dim(2)));
[nxp,nyp]   = size(Xpad);

if mod(nx,2) == 0
    Ix          = [ceil(nxp/2)-nx/2+1:ceil(nxp/2)+nx/2];
else
    Ix          = [1+floor(nxp/2)-floor(nx/2):floor(nxp/2)+1+floor(nx/2)];
end
if mod(ny,2) == 0
    Iy          = [ceil(nyp/2)-ny/2+1:ceil(nyp/2)+ny/2];
else
    Iy          = [1+floor(nyp/2)-floor(ny/2):floor(nyp/2)+1+floor(ny/2)];
end

Xpad(Ix,Iy) = X;




% test
return


X = [1:25];X = reshape(X,5,5);
X = X(1:5,1:4)
dim = [5,5];
Xpad =  GeoStat_padzero(X,dim)