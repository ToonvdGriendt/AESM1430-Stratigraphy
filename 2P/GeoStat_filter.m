function Y = GeoStat_filter(X1,X2,dim_padzeros,inv);

% function Y = GeoStat_filter(X1,X2,dim_padzeros)
%
% This function performs a (multi-dimensional) filter operation
% of X1 and X2. The arrays can be padded with zeros to deal with boundary
% errors using dim_padzeros.
% The centre (ceil(dim(1)/2) , ceil(dim(2)/2)) of X2 is taken as the origin.
% The left upper corner of X1 is its origin.
%
% Illustrated with GeoStat_conv('demo')

if nargin < 4 , inv = []; end

if ~isstr(X1)
    

    % Dimensions
    dim = size(X1); dim(3) = size(X1,3);
    dim2 = size(X2); dim2(3) = size(X2,3);
    if dim ~= dim2,
        error('Dimensions X1 and X2 should be equal')
    end

    % Pad with zeros
    if nargin < 3, dim_padzeros = []; end
    
    
    if isempty(inv)

        if ~isempty(dim_padzeros)
            dim_padzeros = ceil(dim_padzeros/2)*2+1;
            X1pad = X1;
            X1pad(dim_padzeros(1),dim_padzeros(2))= 0 ;
            X2pad = GeoStat_padzero(X2,dim_padzeros);
        else
            X1pad = X1;
            X2pad = X2;
            dim_padzeros = dim*2;
        end

        % Compute the filter operation
        Y = circshift(ifftn(fftn(X1pad).*fftn(X2pad)),[ceil(dim_padzeros(1)/2) ceil(dim_padzeros(2)/2)]);
        Y = Y(1:dim(1),1:dim(2),1:dim(3));

    else


        % Compute the inverse filter
        
        if ~isempty(dim_padzeros)
            dim_padzeros = ceil(dim_padzeros/2)*2+1;
            X1pad = X1;
            X1pad(dim_padzeros(1),dim_padzeros(2))= 0 ;
            X2pad = GeoStat_padzero(X2,dim_padzeros);
        else
            X1pad = X1;
            X2pad = X2;
            dim_padzeros = dim*2;
        end

        % Compute the filter operation
        Y = circshift(ifftn(fftn(X1pad)./fftn(X2pad)),[ceil(dim_padzeros(1)/2) ceil(dim_padzeros(2)/2)]);
        Y = Y(1:dim(1),1:dim(2),1:dim(3));

        
    end

 
else


    % Illustrate
    type  = 'Gaussian';
    dim   = [100 71 1];
    param = [10,20,30,1,0];
    dim_padzeros = [152 130 0];
    [H,h] = GeoStat_Covariance2Filter(type,dim,param);
    h      = circshift(ifft2(H),[ceil(dim(1)/2),ceil(dim(2)/2)]);
    
    A = zeros(size(H));
    A(4,4) = 1; A(end-10,end-10) = 1;A(end,1) = 1; A(10,end) = 1; A(25,35)=1;
    
    
    % Pad with zeros
    if ~isempty(dim_padzeros)
        dim_padzeros = ceil(dim_padzeros/2)*2;
        Apad = A;
        Apad(dim_padzeros(1),dim_padzeros(2))= 0 ;
        hpad = GeoStat_padzero(h,dim_padzeros);
    end
   
    % Compute the filter operation
    Y = circshift(ifftn(fftn(Apad).*fftn(hpad)),[ceil(dim_padzeros/2) ceil(dim_padzeros/2)]);
    D = Y(1:dim(1),2:dim(2)+1);
    
    figure;
    subplot(221);imagesc(h);title('X2')
    subplot(222);imagesc(A);title('X1')
    subplot(223);imagesc(D);title('X1 filtered by X2')
    
    Y = [];
    
end