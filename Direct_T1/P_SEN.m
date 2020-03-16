function [Mo, R1, iter, traces, cost] = P_SEN(P, kU, opt)

options.MaxIter         = opt.Initer;
options.display         = 'off';
options.Method          = 'cg'; % CG in minFunc
options.useMex          = 0;
options.inter           = 1;                % Intermediate flag
options.numDiff         = 0;
options.PROGTOL         = 1e-20;
options.maxFunEvals     = 1000;
options.cgUpdate        = 1;

iter                    = 1;
fres                    = []; % trace of function value

opt.lambda1             = opt.lambdaA(1:2);
opt.lambda2             = opt.lambdaA(3:4);
[P, f, exitflag, output] = minFunc(@P2sig, P(:), options, kU, opt);
P                       = real(P);

fres                    = [fres;output.trace.fval];
traces                  = struct('fres', fres);
disp(['fval: ' num2str(fres(end))]);
disp(output.message);
iter                    = iter + 1;

P = reshape(P, [opt.size(1) opt.size(2) 2]);
Mo = reshape(squeeze(P(:,:,1)),[opt.size(1) opt.size(2)]);
R1 = reshape(squeeze(P(:,:,2)),[opt.size(1) opt.size(2)]);
    
if ~isfield(opt, 'plot')
    if opt.plot
        figure(3);
        imagesc(rot90(cat(2, Mo, 1./R1)), [0 6]); daspect([2 1 1]);
        title('Final estimation');
        colorbar;
        axis off;
    end
end

%% For parameter sweep, 01/08/2019
S                       = SPGR(Mo, R1, opt.B1, opt.FA, opt.tr);
S                       = opt.U.*(1/sqrt(prod(opt.size(opt.FTdim)))).*fFastFT(opt.S.*repmat(S, [1 1 1 1 opt.size(5)]), opt.FTdim, opt.FTshift);
cost1                   = 0.5*sum(abs(S(:)-kU(:)).^2);

if opt.lambda1(1) ~= 0 || opt.lambda1(2) ~= 0
    TV = T(P);
    TV_Mo = TV(:, :, :, 1);
    TV_R1 = TV(:, :, :, 2);
    cost2 = [sum(abs(TV_Mo(:))) sum(abs(TV_R1(:)))];
else
    cost2 = 0;
end

if opt.lambda2(1) ~= 0  || opt.lambda2(2) ~= 0
    [wc_Mo, fsize_Mo] = compWx(Mo, opt);
    [wc_R1, fsize_R1] = compWx(R1, opt);
    cost2 = [sum(abs(wc_Mo(:))) sum(abs(wc_R1(:)))];
else
    cost2 = 0;
end

cost = struct('cost1', cost1, 'cost2', cost2);
end