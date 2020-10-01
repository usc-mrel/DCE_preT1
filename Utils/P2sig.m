function [cost, grad] = P2sig(P, kU, opt)

FA                      = opt.FA;
tr                      = opt.tr;
B1                      = opt.B1;
[np, nv, ns, nt, ~]     = size(kU);

% Calculate cost function
P                       = reshape(P, [np nv ns 2]);
Mo                      = P(:, :, :, 1);
R1                      = P(:, :, :, 2);
S                       = spgr(Mo, R1, B1, FA, tr);
S                       = opt.U .* fFastFT(opt.S.*S, opt);
cost                    = 0.5 * sum(abs(S(:)-kU(:)).^2);

% Calculate gradient
FA                      = FA * pi / 180;
E1                      = exp(-tr.*R1);

g1                      = zeros(np, nv, ns, 2);
dyds                    = sum(conj(opt.S) .* prod(opt.size(opt.FTdim)) .* iFastFT((S-kU), opt), 5);
for it = 1:nt
    g1(:, :, :, 1)      = g1(:, :, :, 1) + conj((1-E1).*sin(B1.*FA(it))./(1-E1.*cos(B1.*FA(it)))) .* dyds(:, :, :, it);
    g1(:, :, :, 2)      = g1(:, :, :, 2) + conj(E1.*(-tr).* Mo.*sin(B1.*FA(it)).*(cos(B1.*FA(it))-1)./((1-E1.*cos(B1.*FA(it))).^2)) .* dyds(:, :, :, it);
end

grad                    = g1;
grad(:, :, :, 2)        = real(grad(:, :, :, 2));
grad                    = grad(:);
end