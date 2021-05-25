function [X_cond,Ct,Cinv,y_error_w,y_condt,x_correctiont,y_error_wt] = GeoStat_condition(h,y_cond,I_y,X_uncond,dim,R_sqrtm,Voidblocks,J,dim_padzero,Ct,Cinv)


%function [B,E,C,theta,theta_long,Vpred,Vresid,orders] = GeoStat_condition(method,h,dim,orders,n_ensemble,method,constraints,Voidblocks)
%
% Condition parameter field to direct measurements of that parameter field
% in a number of locations (e.g. well log data). The covariance of the parameter field
% determines the solution. This is based on the Kriging solution using a
% filter implementation of the covariance properties.
%
% Version 2, June 2006
% S.G. Douma, SHELL-EPT-RXX



% ================= %
% == Check input == %
% ================= %
if nargin < 11, Cinv = []; end
if nargin < 10, Cinv = [];Ct = [];end
if nargin < 9, Cinv = [];Ct = [];J = []; end
if nargin < 8, Cinv = [];Ct = [];J = [];dim_padzero = []; end
if nargin < 7, Cinv = [];Ct = [];J = [];Voidblocks = [];dim_padzero = [];end
if nargin < 6, Cinv = [];Ct = [];J = [];R_sqrtm = sparse(zeros(length(y_cond)));Voidblocks = [];dim_padzero = [];end
if nargin < 5, Cinv = [];Ct = [];J = [];R_sqrtm = sparse(zeros(length(y_cond)));dim = size(X_uncond); dim(3) = size(X_uncond,3);Voidblocks = [];dim_padzero = [];end
  


% ==================== %
% == Initialization == %
% ==================== %
if isempty(dim), dim = size(X_uncond); dim(3) = size(X_uncond,3); end
Nx = dim(1); Ny = dim(2);, Nz = dim(3); shape_array = 0;
if isempty(dim_padzero), dim_padzero = 2*dim; end
if isempty(Voidblocks), Voidblocks = zeros(Nx*Ny*Nz,1); end
if size(X_uncond,1) == Nx,
    shape_array = 1; X_uncond = X_uncond(:);
end
if isempty(R_sqrtm)
    R_sqrtm = sparse(zeros(length(y_cond)));
end
n_ensemble = size(X_uncond,2);
Voidblocks = Voidblocks(:);
I_nonvoid  = find(Voidblocks==0);
I_void     = find(Voidblocks~=0);



% ========================================================== %
% == Compute the covariance function and relevant inverse == %
% ========================================================== %
if isempty(Ct)
    Ct  = GeoStat_filter(h,h,dim_padzero);
    Ctt = [Ct(:);0];
end
if isempty(Cinv)
    if isempty(J),
        I_C = abs(repmat(I_y(:)',length(I_y),1) - repmat(I_y(:),1,length(I_y))) + find(Ctt==max(Ctt));
        I_C(I_C>length(Ctt)) = length(Ctt);
        Cinv = Ctt(I_C)\eye(length(I_y));
    else
        for jj = 1:size(J,2)
            CJt =  GeoStat_filter(reshape(J(:,jj),Nx,Ny),Ct,dim_padzero);
            CJ(:,jj) = CJt(:);
        end
        Cinv = (J'*CJ)\eye(length(I_y));
    end
end



% ================= %
% == Computation == %
% ================= %
   

for k = 1:n_ensemble

    
    % == Extract the present values at measurement locations == %
    y_uncond = X_uncond(I_y,k); 
    
    
    % == Optionally add noise to the measurement values == %
    y_condt = y_cond + R_sqrtm*randn(size(y_cond));
    
 
    % == Compute the weighted measurement error inv(J'*Cyy*J)*y_error
    y_error_w = Cinv*(y_condt - y_uncond);
    
    
    % == Create the measured field of pulses at measurement locations == %
    y_error_wt = zeros(length(Voidblocks),1);
    y_error_wt(I_nonvoid(I_y)) = y_error_w;
    y_error_wt = reshape(y_error_wt,Nx,Ny,Nz);

    
    % == Compute the conditioned realization == % 
    x_correctiont = GeoStat_filter(y_error_wt,Ct,dim_padzero); % Check: shouldn't that be CJ instead of Ct
    x_correction  = x_correctiont(:);
    X_cond(:,k)   = X_uncond(:,k) + x_correction(I_nonvoid);
    
    
end



% ============ %
% == Output == %
% ============ %

if shape_array == 1
    X_cond = reshape(X_cond,Nx,Ny,Nz);
end

