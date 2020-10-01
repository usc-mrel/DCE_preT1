% Script for DRO M0 and T1 values assignment.
% The DRO model is as described in [1][2]. T1 values are adapted from [3]
% for healthy tissues, and from [4-6] for brain tumor.
% [1] Bosca RJ, Jackson EF. Creating an anthropomorphic digital MR phantom - An extensible tool for comparing and evaluating quantitative imaging algorithms. Phys Med Biol. 2016;61(2):974-982. doi:10.1088/0031-9155/61/2/974
% [2] Bliesener Y, Lingala SG, Haldar JP, Nayak KS. Impact of (k,t) sampling on DCE MRI tracer kinetic parameter estimation in digital reference objects. Magn Reson Med. 2020;83(5):1625-1639. doi:10.1002/mrm.28024
% [3] Stanisz GJ, Odrobina EE, Pun J, et al. T1, T2 relaxation and magnetization transfer in tissue at 3T. Magn Reson Med. 2005;54(3):507-512. doi:10.1002/mrm.20605
% [4] Müller A, Jurcoane A, Kebir S, et al. Quantitative T1-mapping detects cloudy-enhancing tumor compartments predicting outcome of patients with glioblastoma. Cancer Med. 2017;6(1):89-99. doi:10.1002/cam4.966
% [5] Hattingen E, Müller A, Jurcoane A, et al. Value of quantitative magnetic resonance imaging T1-relaxometry in predicting contrast-enhancement in glioblastoma patients. Oncotarget. 2017;8(32):53542-53551. doi:10.18632/oncotarget.18612
% [6] Badve C, Yu A, Dastmalchian S, et al. MR fingerprinting of adult brain tumors: Initial experience. In: American Journal of Neuroradiology. Vol 38. American Society of Neuroradiology; 2017:492-499. doi:10.3174/ajnr.A5035
clear; clc; close all;

%% Load DRO data and initilize M0 and T1 matrices
load('discreteBrainModel.mat');
M0 = zeros(size(imData), 'single');
T1 = zeros(size(imData), 'single');

%% Loop all voxels and assign values
for ct = 1:numel(imData)
    switch imData(ct)
        case 0 % Background
            M0(ct) = 0;
            T1(ct) = 0.05;
        case 1 % CSF
            M0(ct) = 1;
            T1(ct) = 2.750;
        case 2 % GM
            M0(ct) = 0.8;
            T1(ct) = 1.820;
        case 3 % WM
            M0(ct) = 0.65;
            T1(ct) = 1.084;
        case 4 % Vessel
            M0(ct) = 1;
            T1(ct) = 1.920;
        case 5 % Skull
            M0(ct) = 0.1;
            T1(ct) = 0.200;
        case 6 % Fat
            M0(ct) = 1;
            T1(ct) = 0.471;
        case 7 % Tumor
            M0(ct) = 0.7;
            T1(ct) = 2;
    end
end

%% Display
load('../Utils/T1cm.mat');
figure(1);
montage(1e3*permute(T1, [1 2 4 3]), 'DisplayRange', [0 3e3], 'Size', [3 4]);
title('T_1 (ms)');
colormap(T1colormap);
colorbar;

figure(2);
montage(permute(M0, [1 2 4 3]), 'DisplayRange', [0 2], 'Size', [3 4]);
title('M_0 (a.u.)');
colormap(gray);
colorbar;

figure(3);
montage(permute(double(imData), [1 2 4 3]), 'DisplayRange', [0 7], 'Size', [3 4]);
title('Tissue type');
colormap jet;
colorbar;
save('discreteBrainModel.mat', 'T1', 'Mo', '-append');