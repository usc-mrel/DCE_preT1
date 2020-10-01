function [M0, R1] = P_SEN(P, kU, opt)
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
%       - M0            Reconstructed amount of magnetization, [np nv ns]
%       - R1            Reconstructed longitudinal relaxation rate, [np nv ns]
%   Others:
%       - np            Number of pointx along x
%       - nv            Number of points along y
%       - ns            Number of slices
%       - nt            Number of contrast, e.g. time fram, flip angle
%       - nr            Number of coils
%       - nu            Number of unknowns

options.MaxIter         = opt.MaxIter;
options.display         = 'off';            % Display intermediate results (2D only)
options.Method          = 'cg';             % CG in minFunc
options.useMex          = 0;
options.inter           = 1;                % Intermediate flag
options.numDiff         = 0;
options.PROGTOL         = 1e-9;
options.maxFunEvals     = 3000;
options.cgUpdate        = 1;

[P, fval, outputinfo]   = argmin(@P2sig, P(:), options, kU, opt);
fres                    = outputinfo.trace;
disp(['fval: ' num2str(fres(end)) ',' num2str(100*fres(end)/(0.5*norm(kU(:))^2)) '%']);
fprintf(outputinfo.msg);

P                       = reshape(P, [opt.size(1) opt.size(2) opt.size(3) 2]);
M0                      = abs(reshape(squeeze(P(:,:,:,1)),[opt.size(1) opt.size(2) opt.size(3)]));
R1                      = real(reshape(squeeze(P(:,:,:,2)),[opt.size(1) opt.size(2) opt.size(3)]));
end