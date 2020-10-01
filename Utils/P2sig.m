function [cost, grad] = P2sig(P, kU, opt)
%   Compute a cost function value and gradient.
%
%   Author: Zhibo Zhu
%   Date: 04/2020
%
%   Usage: [cost, grad] = P2sig(P, kU, opt)
%
%   Input:
%       - P             Parameter vector, [np*nv*ns*nu 1]
%       - kU            (Undersampled) k-space data, [np nv ns nt nr]
%       - opt           Reconstruction option variable, [structure]
%                           opt.FA      Flip angles list in degree
%                           opt.tr      Repitition time in sec
%                           opt.class   Data type
%                           opt.FTdim   Fourier transform dimensions
%                           opt.FTshift Fourier transform shift flag, [0 1]
%                           opt.B1      B1+ map, [np nv ns]
%                           opt.S       Coil sensitivity map, [np nv ns nt nr]
%                           opt.size    Data size
%                           opt.MaxIter Maximum iteration number, [scalar]
%                           opt.U       Sampling pattern matrix, [np nv ns nt nr]
%   Output:
%       - cost          Cost function value, [scalar]
%       - grad          Cost function gradient, [np*nv*ns*nu 1]
%   Others:
%       - np            Number of pointx along x
%       - nv            Number of points along y
%       - ns            Number of slices
%       - nt            Number of contrast, e.g. time fram, flip angle
%       - nr            Number of coils
%       - nu            Number of unknowns

FA                      = opt.FA;
tr                      = opt.tr;
B1                      = opt.B1;
[np, nv, ns, nt, ~]     = size(kU);

% Calculate cost function
P                       = reshape(P, [np nv ns 2]);
M0                      = P(:, :, :, 1);
R1                      = P(:, :, :, 2);
S                       = spgr(M0, R1, B1, FA, tr);
S                       = opt.U .* fFastFT(opt.S.*S, opt);
cost                    = 0.5 * sum(abs(S(:)-kU(:)).^2);

% Calculate gradient
FA                      = FA * pi / 180;
E1                      = exp(-tr.*R1);

g1                      = zeros(np, nv, ns, 2);
dyds                    = sum(conj(opt.S) .* prod(opt.size(opt.FTdim)) .* iFastFT((S-kU), opt), 5);
for it = 1:nt
    g1(:, :, :, 1)      = g1(:, :, :, 1) + conj((1-E1).*sin(B1.*FA(it))./(1-E1.*cos(B1.*FA(it)))) .* dyds(:, :, :, it);
    g1(:, :, :, 2)      = g1(:, :, :, 2) + conj(E1.*(-tr).* M0.*sin(B1.*FA(it)).*(cos(B1.*FA(it))-1)./((1-E1.*cos(B1.*FA(it))).^2)) .* dyds(:, :, :, it);
end

grad                    = g1;
grad(:, :, :, 2)        = real(grad(:, :, :, 2));
grad                    = grad(:);
end