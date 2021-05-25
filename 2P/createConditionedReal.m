function real = createConditionedReal(type,dim,param,y_cond, well_pos,nReal)

% == Select a variogram/covariance and compute associated filter == %
dim_padzeros = [];
[H,h,dim] = GeoStat_Covariance2Filter(type,dim,param);
Nx = dim(1); Ny = dim(2); Nz = dim(3);
figure;imagesc(h);title([type ' filter']);
close



% == Generate conditioned realizations == %
I_y    = zeros(Nx,Ny); 
for i = 1:size(y_cond,1)
    I_y(well_pos(i,1),well_pos(i,2))=1; 
end
I_y = I_y(:); I_y = find(I_y==1);

R_sqrtm  = .1*eye(4); R_sqrtm  = [];
Voidblocks  = [];
J = [];
J = zeros(Nx*Ny,length(I_y)); 
J(I_y,:) = eye(length(I_y));


% figure; 
for k = 1:nReal
    E = randn(size(H));
    X_uncond = ifft2(H.*fft2(E));
    [X_cond,Ct,Cinv,y_error_w,y_condt,x_correctiont,y_error_wt]= GeoStat_condition(h,y_cond,I_y,X_uncond,dim,R_sqrtm,Voidblocks,J,dim_padzeros);
%     subplot(221);imagesc(ifft2(H));colorbar;title('Filter')
%     subplot(222);imagesc(X_uncond);colorbar;title('Unconditioned')
%     subplot(224);imagesc(X_cond);;colorbar;title('Conditioned')
%     subplot(223);imagesc(x_correctiont);colorbar;title('Adaptation')
    X_condt = X_cond(:); X_condt(I_y);  
%     pause(.1)
    real(:,k) = X_condt;
end

