function [Mo, R1, iter, fres, err, corner_energy] = P_SEN(P, sMaps, U1, kU, opt)

options.MaxIter         = opt.Initer;
options.display         = 'off';
options.Method          = 'cg'; % CG in minFunc
options.useMex          = 0;
options.inter           = 0;
options.numDiff         = 0;
options.PROGTOL         = 1e-15;

iter                    = 1;
fres                    = []; % trace of function value
err                     = []; % trace of MSE
corner_energy           = []; % Trace of k-space corner energy;

opt.lambda1             = opt.lambdaA(1);
opt.lambda2             = opt.lambdaA(2);
[P , f, exitflag, output] = minFunc(@P2sig, P(:), options, sMaps, U1, kU, opt);
P                       = real(P);

fres                    = [fres;output.trace.fval];
err                     = [err; output.trace.mse];
corner_energy           = [corner_energy; output.trace.corner_energy];
disp(['fval: ' num2str(fres(end))]);
disp(output.message);
iter                    = iter + 1;

if opt.plot
    P = reshape(P, [opt.size(1) opt.size(2) 2]);
    Mo = reshape(squeeze(P(:,:,1)),[opt.size(1) opt.size(2)]);
    R1 = reshape(squeeze(P(:,:,2)),[opt.size(1) opt.size(2)]);
    figure(3);
    imagesc(rot90(cat(2, Mo, 1./R1)), [0 6]);
    title('Final estimation');
    colorbar;
    axis image off;
end

end

