clear; clc; close all;

%% Load ROI and set matname
load('WM_ROI.mat');
load('GM_ROI.mat');
load('supp.mat');
load('ref.mat');
T1_ref                  = T1;
Mo_ref                  = Mo;
R                       = [32 28 24 20 18 16 14 12 10 9 8 7 6 5 4 3 1]; % Acceleration factor

matnames                = {'Conventional_Full.mat', 'Conventional_rect_R3.mat', ...
    'Conventional_rect_R4.mat', 'Conventional_rect_R5.mat', 'Conventional_rect_R6.mat', ...
    'Conventional_rect_R7.mat', 'Conventional_rect_R8.mat', 'Conventional_rect_R9.mat', ...
    'Conventional_rect_R10.mat', 'Conventional_rect_R12.mat', 'Conventional_rect_R14.mat', ...
    'Conventional_rect_R16.mat', 'Conventional_rect_R18.mat', 'Conventional_rect_R20.mat', ...
    'Conventional_rect_R24.mat', 'Conventional_rect_R28.mat', 'Conventional_rect_R32.mat', ...
    'Direct_ellip_R1.mat', 'Direct_rect_R3.mat', ...
    'Direct_rect_R4.mat', 'Direct_rect_R5.mat', 'Direct_rect_R6.mat', ...
    'Direct_rect_R7.mat', 'Direct_rect_R8.mat', 'Direct_rect_R9.mat'...
    'Direct_rect_R10.mat', 'Direct_rect_R12.mat', 'Direct_rect_R14.mat', ...
    'Direct_rect_R16.mat', 'Direct_rect_R18.mat', 'Direct_rect_R20.mat', ...
    'Direct_rect_R24.mat', 'Direct_rect_R28.mat', 'Direct_rect_R32.mat'};

%%
n                       = size(matnames, 2);
WM_ref                  = WM_ROI.*T1_ref;
GM_ref                  = GM_ROI.*T1_ref;
res                     = zeros(n, 8);
res_esp                 = zeros(100, 2, n);
[nx, ny]                = size(T1_ref);
toShow                  = zeros(ny*2, nx*n/2, 2, 2);
fprintf('%10s %10s %10s %10s %10s %10s %10s %10s\n','WM mean','WM std','WM rmse','WM ssim','GM mean','GM std','GM rmse','GM ssim');

for ct = 1:n
    load(matnames{ct});
    
    WM                  = WM_ROI.*T1;
    GM                  = GM_ROI.*T1;
    
    % Error maps
    toShow((1:ny) + (floor((ct - 0.5)/n*2))*ny, (1:nx) + (ct -1  - floor((ct -0.5)/n*2)*n/2)*nx, 1, 1) = rot90(WM_ROI.*T1);
    toShow((1:ny) + (floor((ct - 0.5)/n*2))*ny, (1:nx) + (ct -1  - floor((ct -0.5)/n*2)*n/2)*nx, 2, 1) = rot90(WM_ROI.*abs(T1-T1_ref)./(abs(T1_ref)));
    toShow((1:ny) + (floor((ct - 0.5)/n*2))*ny, (1:nx) + (ct -1  - floor((ct -0.5)/n*2)*n/2)*nx, 1, 2) = rot90(WM_ROI.*Mo);
    toShow((1:ny) + (floor((ct - 0.5)/n*2))*ny, (1:nx) + (ct -1  - floor((ct -0.5)/n*2)*n/2)*nx, 2, 2) = rot90(WM_ROI.*abs(Mo-Mo_ref)./(abs(Mo_ref)));


    % MSE and SSIM
    [WM_mu, WM_sigma]   = mean_std(WM);
    [GM_mu, GM_sigma]   = mean_std(GM);
    WM_rmse             = rmse(WM_ref(WM_ref~=0), WM(WM~=0));
    GM_rmse             = rmse(GM_ref(GM_ref~=0), GM(GM~=0));
    WM_ssim             = ssim(WM, WM_ref);
    GM_ssim             = ssim(GM, GM_ref);
    
    % ESP
    [esp_plot, grid] = esp(supp.*T1_ref, supp.*T1, [1 2]);
    
    res(ct, :)          = [WM_mu WM_sigma WM_rmse WM_ssim GM_mu GM_sigma GM_rmse GM_ssim];
    res_esp(:, 1, ct)   = esp_plot;
    res_esp(:, 2, ct)   = grid;
    
    fprintf('%10f %10f %10f %10f %10f %10f %10f %10f\n',1000*WM_mu, 1000*WM_sigma, WM_rmse, WM_ssim, 1000*GM_mu, 1000*GM_sigma, GM_rmse, GM_ssim);
end

%% Error maps
figure('Color', 'w', 'Position', [35.8000 131.4000 1.5028e+03 558]);
subplot(2, 1, 1);
imagesc(toShow(:, :, 1, 1), [0 2]);
h = colorbar; h.FontSize = 16;
axis off;
title('T_1 map, SENSE vs. Direct', 'FontSize', 16);
subplot(2, 1, 2);
imagesc(100*toShow(:, :, 2, 1), [0 20]);
h = colorbar; h.FontSize = 16;
axis off;
title('T_1 map error, SENSE vs. Direct', 'FontSize', 16);
print T1 -dtiff

figure('Color', 'w', 'Position', [35.8000 131.4000 1.5028e+03 558]);
subplot(2, 1, 1);
imagesc(toShow(:, :, 1, 2), [0 2]);
h = colorbar; h.FontSize = 16;
axis off;
title('M_0 map, SENSE vs. Direct', 'FontSize', 16);
subplot(2, 1, 2);
imagesc(100*toShow(:, :, 2, 2), [0 20]);
h = colorbar; h.FontSize = 16;
axis off;
title('M_0 map error, SENSE vs. Direct',  'FontSize', 16);
print M0 -dtiff

%%
figure('Color', 'w');
subplot(1, 2, 1); hold on;
title('ESP plot for SENSE', 'FontSize', 16);
for ct = 1:2:size(matnames, 2)/2
    plot(res_esp(:, 2, ct), res_esp(:, 1, ct), 'linewidth', 2);
end
ylim([0 2]);
legend('Full', num2str(R(15)), num2str(R(13)), ...
    num2str(R(11)), num2str(R(9)), ...
    num2str(R(7)), num2str(R(5)), num2str(R(3)), ...
    num2str(R(1)), 'location', 'northeast');

subplot(1, 2, 2); hold on;
title('ESP plot for direct', 'FontSize', 16);
for ct = size(matnames, 2)/2+1:2:size(matnames, 2)
    plot(res_esp(:, 2, ct), res_esp(:, 1, ct), 'linewidth', 2);
end
ylim([0 2]);
legend('Full', num2str(R(15)), num2str(R(13)), ...
    num2str(R(11)), num2str(R(9)), ...
    num2str(R(7)), num2str(R(5)), num2str(R(3)), ...
    num2str(R(1)), 'location', 'northeast');
print ESP -dtiff

figure('Color', 'w');
plot(R, squeeze(res(17:-1:1, 3)), R, squeeze(res(34:-1:18, 3)), R, squeeze(res(17:-1:1, 7)), R, squeeze(res(34:-1:18, 7)), 'linewidth', 2);
ylim([0 1]);
title('RMSE vs R', 'FontSize', 16);
legend('SENSE WM', 'Direct WM', 'SENSE GM', 'Direct GM', 'location', 'northwest');
xlabel('R'); ylabel('RMSE');
print RMSE -dtiff