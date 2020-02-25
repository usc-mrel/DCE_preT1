function [cost, grad] = P2sig(P, kU, opt)

FA                      = opt.FA;
tr                      = opt.tr;
B1                      = opt.B1;
[nx, ny, ns, nt, nr]    = size(kU);

%% Calculate cost function
P                       = reshape(P, [nx ny 2]);
Mo                      = squeeze(P(:, :, 1));
R1                      = squeeze(P(:, :, 2));
S                       = SPGR(Mo, R1, B1, FA, tr);
% S                       = U1.*fft2c(repmat(sMaps, [1 1 1 1 1]).*repmat(S, [1 1 1 1 nr]));
S                       = opt.U.*(1/sqrt(prod(opt.size(opt.FTdim)))).*fFastFT(opt.S.*repmat(S, [1 1 1 1 nr]), opt.FTdim, opt.FTshift);
cost1                   = 0.5*sum(abs(S(:)-kU(:)).^2);

if opt.lambda1 ~= 0 % calculate lambda1*||TVx||1
    TV = compTx(Mo,opt);
    cost2 = opt.lambda1*sum(abs(TV(:)));
else
    cost2 = 0;
end

if opt.lambda2 ~= 0  % calculate lambda2*||Wx||1
    [wc,fsize]=compWx(Mo,opt);
    cost3=opt.lambda2*sum(abs(wc(:)));
else
    cost3=0;
end

cost                    = cost1+cost2+cost3;

%% Calculate gradient
FA                      = FA*pi/180;
E1                      = exp(-tr.*R1);

g1                      = zeros(nx, ny, 2);
% dyds                    = ((sum(repmat(conj(sMaps),[1 1 1 1 1]).*ifft2c(S-kU), 5)));
dyds                    = sum(conj(opt.S).*sqrt(prod(opt.size(opt.FTdim))).*iFastFT((S-kU), opt.FTdim, opt.FTshift), 5);
for it = 1:nt
    g1(:, :, 1) = g1(:, :, 1) + (1-E1).*sin(B1.*FA(it))./(1-E1.*cos(B1.*FA(it))) .* dyds(:, :, :, it);
    g1(:, :, 2) = g1(:, :, 2) + E1.*(-tr).* Mo.*sin(B1.*FA(it)).*(cos(B1.*FA(it))-1)./((1-E1.*cos(B1.*FA(it))).^2) .* dyds(:, :, :, it);
end

if opt.lambda1 ~= 0
    u = 1e-8;
    W = sqrt(conj(TV).*TV+u);
    W = 1./W;
    g2 = opt.lambda1.*compThx(W.*TV,opt);
else
    g2=0;
end

if opt.lambda2 ~= 0
    %[wc,fsize]=compWx(Mo,opt);
    u = 1e-8;
    W1 = sqrt(conj(wc).*wc+u);
    W1 = 1./W1;
    g3 = opt.lambda2.*compWhx(W1.*wc,opt,fsize);
else
    g3 = 0;
end

grad                    = g1+g2+g3;
grad                    = real(grad(:));