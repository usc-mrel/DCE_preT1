% This is the demo script of Sparse Pre-Contrast T1 Mapping for
% High-Resolution Whole-Brain DCE-MRI
clear; clc; close all;

addpath(genpath('./DRO'));
addpath(genpath('./Utils'));

%% Data preparation and reconstruction settings
load('noisecov.mat');
load('discreteBrainModel.mat');
load('coilSenseMaps.mat');

% Set scan parameters and reconstruction related options
FA                      = 1.5 * 10.^((0:6)./6);
B1                      = ones(size(Mo));
TR                      = 0.0049;
opt.FA                  = FA;
opt.tr                  = TR;
opt.class               = 'single';
opt.FTdim               = [1 2];
opt.FTshift             = 1;

% Synthesize fully sampled kspace data and VFA images
k                       = genKspace(Mo, 1./T1, B1, FA, TR, sMaps, [1 2], 1);
img                     = spgr(Mo, 1./T1, B1, FA, TR);

[np, nv, ns, nt, nr]    = size(k);
opt.B1                  = B1;
opt.S                   = repmat(sMaps, [1 1 ns nt 1]);
opt.size                = [np nv ns nt nr];
opt.MaxIter             = 500;

%%
SNR = inf; % SNR level
for realization = 1 % Noise realization loop
    wm                  = img(:, :, :, 7) .* (imData == tissueTypes.WhiteMatter);
    k_noise             = applyNoise(k, wm(abs(wm) ~= 0), noisecov, SNR);
    for pattern = [0] % Pattern loop
        for R = [4] % Undersampling level loop
            fprintf(['Reconstructing SNR: ' num2str(SNR) ', pattern: ' num2str(pattern) ', R: ' num2str(R) '\n']);
            % Apply k-space under-sampling. "pattern' controls the type of pattern and
            % "R" controls under-sampling level.
            [kU, U]     = applyU(k_noise, pattern, R, 1);
            opt.U       = U;
            
%             tmp1        = spgr(ones(1,1,2), ones(1,1,2)/1.5, 1, opt.FA, opt.tr);
%             tmp2        = sum(conj(sMaps).*iFastFT(kU, opt), 5);
%             wm          = tmp2.*(imData==3); % WM signal
%             wm          = wm(abs(wm~=0));
%             sf2         = median(tmp1(:))/median(abs(wm));
            sf2         = 1.5922;
            kU          = kU * sf2; % Data normalization.
            
            % Direct T1 mapping
            Mo          = zeros(np, nv, ns);
            R1          = ones(np, nv, ns);
            P           = cat(4, Mo, R1);
            
            [Mo, R1]    = P_SEN(P, kU, opt);
            Mo          = Mo/sf2;
            T1          = 1./R1;
            
            % Save results
            switch pattern
                case 0
                    matname = ['.\results\SNR' num2str(SNR) '_SP_rect_R' num2str(R) '_re' num2str(realization) '.mat'];
                case 1
                    matname = ['.\results\SNR' num2str(SNR) '_SP_ellip_R' num2str(R) '_re' num2str(realization) '.mat'];
                case 2
                    matname = ['.\results\SNR' num2str(SNR) '_RGA_ellip_R' num2str(R) '_re' num2str(realization) '.mat'];
                case 3
                    matname = ['.\results\SNR' num2str(SNR) '_RGA_rect_R' num2str(R) '_re' num2str(realization) '.mat'];
            end
            save(matname, 'Mo', 'R1', 'T1', 'U');
        end
    end
end