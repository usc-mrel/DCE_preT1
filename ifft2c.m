function res = ifft2c(x)
%
%
% res = ifft2c(x)
% 
% orthonormal centered 2D ifft
%
% (c) Michael Lustig 2005

sz = size(x);
res = sqrt(prod(sz(1:2)))*fftshift(ifft2(ifftshift(x)));

