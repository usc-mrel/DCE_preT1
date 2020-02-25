function [Mo, R1, iter, traces, cost] = P_SEN(P, kU, opt)

options.MaxIter         = opt.Initer;
options.display         = 'off';
options.Method          = 'cg'; % CG in minFunc
options.useMex          = 0;
options.inter           = 1;                % Intermediate flag
options.numDiff         = 0;
options.PROGTOL         = 1e-20;
options.maxFunEvals     = 3000;
options.cgUpdate        = 1;

iter                    = 1;
fres                    = []; % trace of function value

opt.lambda1             = opt.lambdaA(1);
opt.lambda2             = opt.lambdaA(2);
[P, f, exitflag, output] = minFunc(@P2sig, P(:), options, kU, opt);
P                       = real(P);

fres                    = [fres;output.trace.fval];
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
    daspect([2 1 1]);
    colorbar;
    axis image off;
end

end

