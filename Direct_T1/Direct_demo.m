clear; clc; close all;

%% Load files and preprocessing
load('k_x90.mat');
load('sMaps_old.mat');
load('B1_x90.mat');
load('supp.mat');
load('ellipse.mat');

k                       = 500*k/max(abs(k(:)));
[nx, ny, ns, nt, nr]    = size(k);
k                       = k/sqrt(nx*ny);
% imgF                    = sum(repmat(conj(sMaps), [1 1 1 1 1]).*ifft2c(k), 5);

opt.Initer              = 500;
opt.Outiter             = 10;

%% Fully sampled maps
refmat                  = 'ref.mat';
no_ref                  = ~exist(refmat);
FA                      = 1.5*10.^((0:6)./6);
opt.FA                  = FA;
opt.tr                  = 0.0049;
opt.B1                  = B1;
opt.class               = 'single';
if no_ref
    % Estimate references from fully sampled data using direct estimation
    Mo                  = zeros(nx, ny);
    R1                  = ones(nx, ny);
    P                   = Mo;
    P(:, :, 2)          = R1;
    tic,
    [Mo, R1, iter, trace, cost] = P_SEN(P, sMaps, 1, k, opt);
    toc,
    T1                  = 1./R1;
    
    matname             = './ref.mat';
    save(matname, 'Mo', 'T1');
else
    load(refmat);
end

%% Setting opt
opt.size                = [nx ny ns nt nr];
opt.plot                = 1;
opt.disp_err            = 1;
opt.Mo                  = squeeze(Mo);
opt.T1                  = squeeze(T1);
opt.supp                = supp;
opt.S                   = sMaps;
opt.kfov                = repmat(kmask, [1 1 1 nt nr]);
opt.pattern_type        = 2; % Sampling pattern footprint flag. 1: Rectangular. 2: Ellipsoid
opt.FTdim           = [1 2];
opt.FTshift           = 1;

% Optional constraints.
opt.wname               = 'db4'; % wavelet parameters
opt.worder              = {[1 2],[1 2],[1 2]};
opt.lambdaA             = [0 0 0 0];

imgF                    = sqrt(sum(abs(sqrt(prod(opt.size(opt.FTdim)))*iFastFT(k, [1 2], 1)).^2, 5));

%% Reconstruction script
switch opt.pattern_type
    case 1 % Rectangular
        load('Uyz1.mat');
        Uyz = Uyz1;
        clear Uyz1;
        prefix = './Direct_T1/results/Direct_rect_R';
    case 2 % Ellipsoid
        load('Uyz2.mat');
        Uyz = Uyz2;
        clear Uyz2;
        prefix = './Direct_T1/results/Direct_ellip_R';
end

%% Reconstruction
range = 1:17;
for ct = range % For 17 subsampling level in factor decreasing order
    U1                  = logical(repmat(Uyz(:, :, :, :, ct), [1 1 1 1 nr]));
%     U1                  = double(logical(nshift(U1, [1 2])));

    kU                  = k.*U1;
    if opt.pattern_type == 2 % For ellipsoid footprint
        U1(~opt.kfov & U1 == 0) = 0.1;
    end
    

    opt.U               = U1;
    
%     imgZF               = sum(repmat(conj(sMaps), [1 1 1 1 1]).*ifft2c(kU), 5);
    imgZF               = sum(conj(opt.S).*sqrt(prod(opt.size(opt.FTdim))).*iFastFT(kU, opt.FTdim, opt.FTshift), 5);

    figure(21);
    montage(abs(rot90(imgZF)), 'DisplayRange', [0 max(abs(imgZF(:)))], 'Size', [2 4]);
    daspect([2 1 1]);
    title('Zero-filled VFA images');

    % Direct reconstruct Mo and T1 jointly
    Mo                  = ones(nx, ny);
    R1                  = ones(nx, ny);
    P                   = Mo;
    P(:, :, 2)          = R1;

    tic,
    [Mo, R1, iter, traces, cost] = P_SEN(P, kU, opt);
    toc,
    T1                  = 1./R1;
    
    matname             = [prefix num2str(R(ct))];
    save(matname, 'Mo', 'T1', 'traces');

    if opt.disp_err
        l = length(traces.fres);
        % Objective function value
        figure(24);
        semilogy(1:l, traces.fres, 'linewidth', 2);
        title('Objective function changes over iteration');
        xlabel('Iter #'); ylabel('Objective function value');
    end
end
