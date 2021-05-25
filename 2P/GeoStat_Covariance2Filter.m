function [H,h,dim,XYCovariogram] = GeoStat_Covariance2Filter(type,dim,param)

% function [H,H_AR,XYCovariogram] = GeoStat_Covariance2Filter(type,dim,param)
%
% This function provides the frequency response H of a linear filter,
% corresponding to a covariance function C(tau) = c*(1-M)/(2-M)
% The filter will be centred around (ceil(dim(1)/2),ceil(dim(2)/2)).
%
% Input:
% type              'Gaussian'      : The Gaussian filter corresponds to C(tau) with M = 1-exp(-(t/a).^2)
%                   'Exponential'   : The exponential filter corresponds to C(tau) with M = 1-exp(-abs(t/a))
%                   'Spherical'     : The spherical filter corresponds to C(tau) with M = 3*t/2*a + t^3/(2*a^3) for t < a and C(tau)= 0        , t => a
%                   
% dim               Dimensions of the parameter field
% param             [SigmaX,SigmaY,Teta,Variance,RNugget]
%
% Output:
% H                 Frequency response of the filter
% h                 Pulse response of the filter with origin in the centre
% XYCovariogram     Corresponding covariogram
%
% Version 1, June 2006
% S.G. Douma, SHELL-EPT-RXX



% Initialize, make sure the dimensions are uneven to have a unique centre.
Nz = floor(dim(3)/2)*2+1; 
Ny = floor(dim(2)/2)*2+1;
Nx = floor(dim(1)/2)*2+1;
dim = [Nx,Ny,Nz];

% Set the type
covariogram{1} = Covariogram(param(1),param(2),param(3),param(4),param(5),type);


% Compute the covariance
gridLayout.DeltaX = 1;
gridLayout.DeltaY = 1;
if mod(Nx,2) == 0
    xArray = gridLayout.DeltaX .* (mod((0:Nx-1)-1+Nx/2,Nx)+1-Nx/2);
else
    xArray = gridLayout.DeltaX .* (mod((0:Nx-1)+(Nx-1)/2,Nx)-(Nx-1)/2);
end
X = xArray' * ones(1,Ny);
if mod(Ny,2) == 0
    yArray = gridLayout.DeltaY .* (mod((0:Ny-1)-1+Ny/2,Ny)+1-Ny/2);
else
    yArray = gridLayout.DeltaY .* (mod((0:Ny-1)+(Ny-1)/2,Ny)-(Ny-1)/2);
end
Y = ones(Nx,1) * yArray;
XYCovariogram = covariance(covariogram{1},X,Y);


% Compute the corresponding filter
H     = fft2(ifft2(sqrt(abs(fft2(XYCovariogram)))));
h     = fftshift(ifft2(H));





% Test the function
return
type  = 'Gaussian'
dim   = [100 70 1];
param = [15,25,30,1,0];
[H,h,dim] = GeoStat_Covariance2Filter(type,dim,param);

A = zeros(size(H));
A(10,10) = 1; A(end-10,end-10) = 1;A(end-10,10) = 1; A(10,end-10) = 1; A(30,25)=1;

% map
figure;
subplot(221);imagesc(ifft2(H));
subplot(222);imagesc(A);
subplot(223);imagesc(ifft2(H.*fft2(A)));

% inverse
figure;
D3 = ifft2(H.*fft2(A));
subplot(321);imagesc(ifft2(H))
subplot(322);imagesc(A)
subplot(323);imagesc(ifft2(H.*fft2(A)));
subplot(324);imagesc(ifft2(fft2(D3)./H))
subplot(325);imagesc(ifft2(fft2(D3)./fft2(A)));


% Pad with zeros
hpad = zeros(dim(1)*2,dim(2)*2);
hpad(floor(end/4)+1:floor(end/4)+dim(1),floor(end/4)+1:floor(end/4)+dim(2)) = h;
apad = zeros(dim(1)*2,dim(2)*2);
apad(ceil(end/4)+1:ceil(end/4)+dim(1),ceil(end/4)+1:ceil(end/4)+dim(2)) = A;
figure;
subplot(221);imagesc(hpad)
subplot(222);imagesc(apad)
dpad = circshift(ifft2(fft2(hpad).*fft2(apad)),[ceil(dim(1)/2) ceil(dim(2)/2)]);
dpad = dpad(1:dim(1),1:dim(2))
subplot(223);imagesc(dpad);
subplot(224);imagesc(ifft2(H.*fft2(A)));



Y = GeoStat_filter(A,h,dim_padzeros)
figure;imagesc(Y);


% new realizations
figure;
for k = 1:50

    E = randn(size(H));
    subplot(221);imagesc(ifft2(H));colorbar
    subplot(222);imagesc(E);colorbar
    subplot(224);imagesc(ifft2(H.*fft2(E)));;colorbar
    pause
    
end


% Condition the realization
Nx = dim(1);Ny = dim(2);Nz = dim(3);
y_cond      = 4*ones(4,1);
I_y = zeros(Nx,Ny); 
I_y(15,15)=1;I_y(35,15)=1;I_y(15,65)=1;I_y(75,65)=1;
I_y = I_y(:);I_y = find(I_y==1);
R_sqrtm     = [];
Voidblocks  = [];
J = [];
J = zeros(Nx*Ny,length(I_y)); 
J(I_y,:) = eye(length(I_y));
dim_padzero = [];

figure;
for k = 1:15
    
    E = randn(size(H));
    X_uncond = ifft2(H.*fft2(E));
    [X_cond,Ct,Cinv,y_error_w,y_condt,x_correctiont,y_error_wt]= GeoStat_condition(h,y_cond,I_y,X_uncond,dim,R_sqrtm,Voidblocks,J,dim_padzero);
    subplot(221);imagesc(ifft2(H));colorbar
    subplot(222);imagesc(X_uncond);colorbar
    subplot(223);imagesc(X_cond);;colorbar
    subplot(224);imagesc(x_correctiont);colorbar
    
    X_condt = X_cond(:);
    X_condt(I_y)
  
    pause
    
end



figure;imagesc(h);colorbar
figure
subplot(222);surf(X_cond);hold on;plot3(repmat([1:Ny]',1,Nx),repmat([1:Nx]',1,Ny)',y_condt(:,:,1)','x','Linewidth',3)
subplot(221);imagesc(X_cond);colorbar;
subplot(223);imagesc(d);colorbar%imagesc(filter2((h_ar),y_error));;colorbar
subplot(224);imagesc(X_uncond);colorbar%imagesc(filter2((h_ar),y_error));;colorbar

subplot(224);imagesc(y_error)


% autoregressive
E    = randn(size(H));
H = H;
Y    = ifft2(H.*fft2(E));
h_ar = ifft2(1-1./H);

h_ar2 = ifft2(1./H);
h_ar2(1,1) = 0;
figure;imagesc(h_ar2)
Y2   = ifft2(fft2(Y).*fft2(h_ar2));

figure;
subplot(221);imagesc(h_ar2);
subplot(222);imagesc(Y);
subplot(223);imagesc(Y2);
subplot(224);imagesc(Y2+E);
norm(Y2+E-Y)

figure;
subplot(221);surf(h_ar2);
subplot(222);surf(Y);
subplot(223);surf(Y2);
subplot(224);surf(Y2+E);
norm(Y2+E-Y)

figure;
subplot(221);surf(Y+ifft2(H.*fft2(y_error)));
subplot(222);surf(Y);
subplot(223);surf(Y2);
subplot(224);surf(Y2+E);
norm(Y2+E-Y)
figure;plot([Y(:) Y2(:)])




