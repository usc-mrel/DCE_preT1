clear; clc; close all;

%% Load data
load('k_x90.mat');
load('sMaps_old.mat');
load('B1_x90.mat');

%% Preprocessing
[nx, ny, ns, nt, nr]    = size(k);
imgFA                   = zeros(nx, ny, nt);
img0                    = zeros(nx, ny, nt);
k                       = 500*k/max(abs(k(:)));
k                       = k/sqrt(nx*ny);

%% Setting opt
opt.size                = [nx ny ns nt nr];
opt.B1                  = B1;
opt.class               = 'double';
opt.tr                  = 0.0049;
opt.FA                  = 1.5*10.^((0:6)./6);
opt.plot                = 1;
opt.pattern_type        = 2;
opt.FTdim               = [1 2];
opt.FTshift             = 1;

imgF                    = sqrt(sum(abs(sqrt(prod(opt.size(opt.FTdim)))*iFastFT(k, [1 2], 1)).^2, 5));

%% SENSE-CG reconstruction
%  Full sampling
[T1, Mo]                = iDESPOT1(real(imgF), opt);
save('./SENSE+DESPOT/results/Conventional_ellip_R1.mat', 'T1', 'Mo');

switch opt.pattern_type
    case 1 % Rectangular
        load('Uyz1.mat');
        Uyz = Uyz1;
        clear Uyz1;
        prefix = './SENSE+DESPOT/results/Conventional_rect_R';
    case 2 % Ellipsoid
        load('Uyz2.mat');
        Uyz = Uyz2;
        clear Uyz2;
        prefix = './SENSE+DESPOT/results/Conventional_ellip_R';
end

%%
range = 1:17;
for ct = range
    imgFA               = zeros(nx, ny, nt);
    for int = 1:nt
        U               = repmat(Uyz(:, :, :, int, ct), [1 1 nr]);
        U               = logical(U);
        kU              = U.*squeeze(k(:, :, :, int, :));
        S               = squeeze(sMaps(:, :, :, int, :));
%         img0(:, :, int) = sum(conj(S).*ifft2c(kU), 3);
        img0(:, :, int) = sum(conj(S).*sqrt(prod(opt.size(opt.FTdim))).*iFastFT(kU, opt.FTdim, opt.FTshift), 3);
        imgR            = SENSE_CG_recon(S, U, kU, opt);
        imgFA(:, :, int)= imgR;
    end

    % iDESPOT1
    [T1, Mo]            = iDESPOT1(real(permute(imgFA, [1 2 4 3])), opt);
    matname             = [prefix num2str(R(ct))];
    save(matname, 'Mo', 'T1');
    
    if opt.plot
        figure(11);
        montage(rot90(abs(permute(img0, [1 2 4 3]))), 'DisplayRange', [0 max(abs(img0(:)))]);
        daspect([2 1 1]);
        title('Zero-filled VFA images');
        figure(12);
        montage(rot90(abs(permute(imgFA, [1 2 4 3]))), 'DisplayRange', [0 max(abs(imgFA(:)))]);
        daspect([2 1 1]);
        title('Reconstructed VFA images');
        figure(13);
        imagesc(rot90(cat(2, Mo, T1)), [0 6]);
        daspect([2 1 1]);
        title('Final estimation');
        colorbar;
        axis off;
    end
end
