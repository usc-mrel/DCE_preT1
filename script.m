clear; clc; close all;

addpath(genpath('../DRO stuff'));
addpath(genpath('../DCE_direct_recon'));
addpath(genpath('../STARDCE'));

% Remove some paths to avoid conflicts
rmpath(genpath('../DCE_direct_recon/minFunc_2012'));
rmpath(genpath('../STARDCE/T1estimation/Undersampling/SP'));

%%
load('noisecov.mat');
load('discreteBrainModel.mat');
load('coilSenseMaps.mat');

FA                      = 1.5 * 10.^((0:6)./6);
B1                      = ones(size(Mo));
TR                      = 0.0049;
opt.FA                  = FA;
opt.tr                  = TR;
opt.plot                = 0;
opt.class               = 'single';
opt.FTdim               = [1 2];
opt.FTshift             = 1;

sf1 = 1;
k = genKspace(sf1*Mo, 1./T1, B1, FA, TR, sMaps, [1 2], 1);
img = spgr(Mo, 1./T1, B1, FA, TR);

[np, nv, ns, nt, nr]    = size(k);
opt.B1                  = B1;
opt.S                   = repmat(sMaps, [1 1 ns nt 1]);
opt.size                = [np nv ns nt nr];
opt.MaxIter             = 500;

% Unused spatial constraints
opt.wname               = 'db4'; % wavelet parameters
opt.worder              = {[1 2],[1 2],[1 2]};
opt.lambda              = [0 0 0 0];
opt.disp_err            = 1;

%% Retrospective data generation
SNR = 50;
for realization = 1
    wm = img(:, :, :, 7).*(imData == tissueTypes.WhiteMatter);
    k_noise = applyNoise_new(k, wm(abs(wm) ~= 0), noisecov, SNR);
    for pattern = [0 1 2 3]
        for R = [1 4 7 10 16 22 28 34 40]
            fprintf(['Reconstructing SNR: ' num2str(SNR) ', pattern: ' num2str(pattern) ', R: ' num2str(R) '\n']);
            % Apply k-space under-sampling. "pattern' controls the type of pattern and
            % "R" controls under-sampling level.
            [kU, U] = applyU(k_noise, pattern, R, 1);
            
            opt.U = U;
            sf2 = 1.5922;
            kU = kU * sf2;
            
            % Direct T1 mapping
            Mo                  = zeros(np, nv, ns);
            R1                  = ones(np, nv, ns);
            P                   = Mo;
            P(:, :, :, 2)       = R1;
            
            [Mo, R1]            = P_SEN(P, kU, opt);
            Mo                  = Mo/sf2/sf1;
            T1                  = 1./R1;
            
            switch pattern
                case 0
                    matname = ['D:\mrel_data\sparseT1\DRO\noisy\SNR' num2str(SNR) '\SP_rect_R' num2str(R) '_re' num2str(realization) '.mat'];
                case 1
                    matname = ['D:\mrel_data\sparseT1\DRO\noisy\SNR' num2str(SNR) '\SP_ellip_R' num2str(R) '_re' num2str(realization) '.mat'];
                case 2
                    matname = ['D:\mrel_data\sparseT1\DRO\noisy\SNR' num2str(SNR) '\RGA_ellip_R' num2str(R) '_re' num2str(realization) '.mat'];
                case 3
                    matname = ['D:\mrel_data\sparseT1\DRO\noisy\SNR' num2str(SNR) '\RGA_rect_R' num2str(R) '_re' num2str(realization) '.mat'];
            end
            save(matname, 'Mo', 'R1', 'T1', 'U');
        end
    end
end