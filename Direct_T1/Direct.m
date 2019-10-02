clear; clc; close all;
% addpath(genpath('/Users/zhibozhu/Documents/DCE_direct_recon'));
% addpath(genpath('/Users/zhibozhu/Documents/MATLAB/DCE/STARDCE-master-NEW'));
% addpath('/Users/zhibozhu/Desktop/parameter sweep');

%% Load files and preprocessing
load('k_x90.mat');
load('S_x90.mat');
load('sMaps.mat');
load('B1_x90.mat');
load('supp.mat');
load('ellipse.mat');
k                       = 500*k/max(abs(k(:)));
[nx, ny, ns, nt, nr]    = size(k);
k                       = k/sqrt(nx*ny);
imgF                    = sum(repmat(conj(sMaps), [1 1 1 1 1]).*ifft2c(k), 5);

opt.Initer              =100;
opt.Outiter             =10;

%% Fully sampled maps
no_ref                  = 1;
if no_ref
    FA                  = 1.5*10.^((0:6)./6);
    opt.FA              = FA;
    opt.tr              = 0.0049;
    opt.B1              = B1;
    opt.class           = 'single';
    opt.no_reg          = 1;
    [T1, Mo]            = iDESPOT1(real(imgF),opt);
    ind                 = (T1 < 0.05) | (Mo < 0);
    T1(ind)             = 0.05;
    Mo(ind)             = 0;
end

%% Setting opt
opt.size                = [nx ny ns nt nr];
FA                      = 1.5*10.^((0:6)./6);
opt.FA                  = FA;
opt.tr                  = 0.0049;
opt.B1                  = B1;
opt.wname               = 'db4'; % wavelet parameters
opt.worder              = {[1 2],[1 2],[1 2]};
opt.lambdaA             = [0 0 0 0];
opt.plot                = 1;
opt.disp_err            = 1;
opt.Mo                  = squeeze(Mo);
opt.T1                  = squeeze(T1);
opt.supp                = supp;
opt.sMaps               = sMaps;
opt.kfov                = repmat(nshift(kmask, [1 2]), [1 1 1 nt nr]);

%% Sampling with radial spokes reaching rectangular
load('Uyz1.mat');
for ct = 1:16 % For 17 subsampling level in factor decreasing order
    U1                  = repmat(Uyz1(:, :, :, :, ct), [1 1 1 1 nr]);
    U1                  = logical(nshift(U1, [1 2]));

    kU                  = k.*U1;
    imgZF               = sum(repmat(conj(sMaps), [1 1 1 1 1]).*ifft2c(kU), 5);

    figure(1);
    montage(abs(rot90(imgZF)), 'DisplayRange', [0 max(abs(imgZF(:)))]);
    title('Zero-filled VFA images');

    % Direct reconstruct Mo and T1 jointly
    Mo                  = ones(nx, ny);
    R1                  = ones(nx, ny);
    P                   = Mo;
    P(:, :, 2)          = R1;
    tic,
    [Mo, R1, iter, fres, err, corner] = P_SEN(P, sMaps, U1, kU, opt);
    toc,
    T1                  = 1./R1;
    
    matname             = ['./Direct_T1/Direct_rect_R' num2str(R(ct))];
    save(matname, 'Mo', 'T1', 'err', 'fres', 'corner');
    if opt.disp_err
        figure(4);
        semilogy(fres, 'linewidth', 2);
        title('Objective function changes over iteration');
        figure(5);
        semilogy(err,  'linewidth', 2);
        title('MSE changes over iteration');
        figure(6);
        plot(corner,  'linewidth', 2);
        title('k-space corner energy evolution');
    end
end

%% Sampling with radial spokes reaching ellipse, for corner energy track
load('Uyz2.mat');
for ct = 9:17 % For 17 subsampling level in factor decreasing order
    U1                  = repmat(Uyz2(:, :, :, :, ct), [1 1 1 1 nr]);
    U1                  = logical(nshift(U1, [1 2]));

    kU                  = k.*U1;
    imgZF               = sum(repmat(conj(sMaps), [1 1 1 1 1]).*ifft2c(kU), 5);

    figure(1);
    montage(abs(rot90(imgZF)), 'DisplayRange', [0 max(abs(imgZF(:)))]);
    title('Zero-filled VFA images');

    % Direct reconstruct Mo and T1 jointly
    Mo                  = ones(nx, ny);
    R1                  = ones(nx, ny);
    P                   = Mo;
    P(:, :, 2)          = R1;
    tic,
    [Mo, R1, iter, fres, err, corner] = P_SEN(P, sMaps, U1, kU, opt);
    toc,
    T1                  = 1./R1;
    
    k_hat               = fft2c(sMaps.*SPGR(Mo, R1, opt.B1, opt.FA, opt.tr));
    img_hat             = sum(conj(sMaps).*ifft2c(opt.kfov.*k_hat), 5);
    [T1, Mo]            = iDESPOT1(real(img_hat), opt);
    matname             = ['./Direct_T1/Direct_ellip_R' num2str(R(ct))];
    save(matname, 'Mo', 'T1', 'err', 'fres', 'corner');
    if opt.disp_err
        figure(4);
        semilogy(fres);
        title('Objective function changes over iteration');
        figure(5);
        semilogy(err);
        title('MSE changes over iteration');
        figure(6);
        plot(corner,  'linewidth', 2);
        title('k-space corner energy evolution');
    end
end