% This is the demo script of Sparse Pre-Contrast T1 Mapping for
% High-Resolution Whole-Brain DCE-MRI
clear; clc; close all;

addpath(genpath('../DRO'));
addpath(genpath('../Utils'));

field_flag = 0;

%% Data preparation and reconstruction settings
load('discreteBrainModel.mat');
load('coilSenseMaps.mat');
if field_flag == 0
    load('noisecov.mat');
else
    load('noisecov_3T.mat');
end

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
SNR = 50;                              % SNR level
realization = 1;                        % Noise realization indicator
pattern = 0;                            %Pattern type indicator
% 0, Cartesian spiral, rectangular footprint.
% 1, Cartesian spiral, elliptical footprint.
% 2, Cartesian-radial Randomized Golden Angle, elliptical footprint.
% 3, Cartesian-radial Randomized Golden Angle, rectangular footprint.
R = 1;                                  % Undersampling level
wm                  = img(:, :, :, 7) .* (imData == tissueTypes.WhiteMatter);
k_noise             = addNoise(k, wm(abs(wm) ~= 0), noisecov, SNR, field_flag);

fprintf(['Reconstructing SNR: ' num2str(SNR) ', pattern: ' num2str(pattern) ', R: ' num2str(R) '\n']);
% Apply k-space under-sampling. "pattern' controls the type of pattern and
% "R" controls under-sampling level.
[kU, U]     = applyU(k_noise, pattern, R, 1);
opt.U       = U;

% tmp1        = spgr(ones(1,1,2), ones(1,1,2)/1.084, 1, opt.FA, opt.tr);
% tmp2        = sum(conj(sMaps).*iFastFT(kU, opt), 5);
% wm          = tmp2.*(imData==3); % WM signal
% wm          = wm(abs(wm)~=0);
% sf2         = median(tmp1(:))/median(abs(wm));
sf2         = 1.5384;
kU          = kU * sf2;                 % Data normalization.

% Direct T1 mapping
m0          = zeros(np, nv, ns);
r1          = ones(np, nv, ns);
P           = cat(4, m0, r1);

[m0, r1]    = P_SEN(P, kU, opt);
m0          = m0/sf2;
t1          = 1./r1;

% Save results
switch pattern
    case 0
        matname = ['../results/15/SP_rect_R' num2str(R) '_re' num2str(realization) '.mat'];
    case 1
        matname = ['../results/15/SP_ellip_R' num2str(R) '_re' num2str(realization) '.mat'];
    case 2
        matname = ['../results/15/RGA_ellip_R' num2str(R) '_re' num2str(realization) '.mat'];
    case 3
        matname = ['../results/15/RGA_rect_R' num2str(R) '_re' num2str(realization) '.mat'];
end
save(matname, 'm0', 'r1', 't1', 'U');