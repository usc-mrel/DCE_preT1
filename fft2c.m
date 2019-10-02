function res = fft2c(x)

% res = fft2c(x)
% 
% orthonormal forward 2D FFT
%
% (c) Michael Lustig 2005

sz = size(x);
res = 1/sqrt(prod(sz(1:2)))*fftshift(fft2(ifftshift(x)));

