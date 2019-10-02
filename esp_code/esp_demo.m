clear;
close all;
clc;

%% Create gold standard and corrupted images
N = 256;
testIm = phantom(256); % gold standard

win = gausswin(N,10)*gausswin(N,10)';
version1 = fftshift(ifft2(ifftshift(win.*fftshift(fft2(ifftshift(testIm)))))); % First image
nrmse = norm(version1(:)-testIm(:))/norm(testIm(:));

noise = complex(randn(N),randn(N));
noise = noise*nrmse/norm(noise(:))*norm(testIm(:));
version2 = phantom+noise; % Second image
nrmse2 = norm(version2(:)-testIm(:))/norm(testIm(:));

%% Display images
figure;imagesc([abs(testIm) abs(version1) abs(version2)]);axis equal;axis tight;colormap(gray(256));axis off;caxis([0,1]);
title('Gold standard (left) together with two corrupted images');
figure;imagesc([abs(testIm-testIm) abs(testIm-version1) abs(testIm-version2)]);axis equal;axis tight;colormap(gray(256));axis off;caxis([0,0.5]);
title('Error images corresponding to the previous figure');

%% Display ESPs (Plotting in Nyquist units)
[esp_plot, grid, lambda] = esp(testIm,version1); % Compute ESP for the first image, using cross-validation to estimate lambda
figure;
plot(grid,esp_plot,'linewidth',3);
hold on;
esp(testIm,version2,[],lambda); % Compute ESP for the second image, reusing the regularization parameter from the first case
legend('First image','Second image');

%% Display ESPs (Plotting in mm)
[esp_plot, grid, lambda] = esp(testIm,version1,[2 2]); % Compute ESP for the first image, using cross-validation to estimate lambda (assuming 2mm isotropic resolution)
figure;
plot(grid,esp_plot,'linewidth',3);
hold on;
esp(testIm,version2,[2 2],lambda); % Compute ESP for the second image, reusing the regularization parameter from the first case
legend('First image','Second image');