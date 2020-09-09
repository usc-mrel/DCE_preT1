clear; clc; close all;

%%
load('discreteBrainModel.mat');
Mo = zeros(size(imData));
T1 = zeros(size(imData));

%%
for ct = 1:numel(imData)
    switch imData(ct)
        case 0 % Background
            Mo(ct) = 0;
            T1(ct) = 0.05;
        case 1 % CSF
            Mo(ct) = 1;
            T1(ct) = 2.750;
        case 2 % GM
            Mo(ct) = 0.8;
            T1(ct) = 1.820;
        case 3 % WM
            Mo(ct) = 0.65;
            T1(ct) = 1.084;
        case 4 % Vessel
            Mo(ct) = 1;
            T1(ct) = 1.920;
        case 5 % Skull
            Mo(ct) = 0.1;
            T1(ct) = 0.200;
        case 6 % Fat
            Mo(ct) = 1;
            T1(ct) = 0.471;
        case 7 % Tumor
            Mo(ct) = 0.7;
            T1(ct) = 2;
    end
end

%%
figure(1);
montage(permute(T1, [1 2 4 3]), 'DisplayRange', [0 3], 'Size', [3 4]);
title('T_1 (sec)');
colormap jet;
colorbar;

figure(2);
montage(permute(Mo, [1 2 4 3]), 'DisplayRange', [0 2], 'Size', [3 4]);
title('M_0 (a.u.)');
colormap jet;
colorbar;

figure(3);
montage(permute(double(imData), [1 2 4 3]), 'DisplayRange', [0 7], 'Size', [3 4]);
title('Tissue tupe');
colormap jet;
colorbar;

save('discreteBrainModel.mat', 'T1', 'Mo', '-append');