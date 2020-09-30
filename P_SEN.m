function [Mo, R1] = P_SEN(P, kU, opt)

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
Mo                      = abs(reshape(squeeze(P(:,:,:,1)),[opt.size(1) opt.size(2) opt.size(3)]));
R1                      = real(reshape(squeeze(P(:,:,:,2)),[opt.size(1) opt.size(2) opt.size(3)]));
end