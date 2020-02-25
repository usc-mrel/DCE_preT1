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

if opt.lambda1(1) ~= 0 || opt.lambda1(2) ~= 0 % calculate lambda1*||TVx||1
    TV = T(P);
    TV_Mo = TV(:, :, :, 1);
    TV_R1 = TV(:, :, :, 2);
    cost2 = opt.lambda1(1)*sum(abs(TV_Mo(:))) + opt.lambda1(2)*sum(abs(TV_R1(:)));
else
    cost2 = 0;
end

if opt.lambda2(1) ~= 0  || opt.lambda2(2) ~= 0% calculate lambda2*||Wx||1
    [wc_Mo, fsize_Mo] = compWx(Mo, opt);
    [wc_R1, fsize_R1] = compWx(R1, opt);
    cost3 = opt.lambda2(1)*sum(abs(wc_Mo(:))) + opt.lambda2(2)*sum(abs(wc_R1(:)));
else
    cost3 = 0;
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

if opt.lambda1(1) ~= 0 || opt.lambda1(2) ~= 0
    u = 1e-6;
    W_Mo = sqrt(conj(TV_Mo).*TV_Mo + u);
    W_R1 = sqrt(conj(TV_R1).*TV_R1 + u);
    W_Mo = 1./W_Mo;
    W_R1 = 1./W_R1;
    g2 = -cat(3, opt.lambda1(1)*Th(W_Mo.*TV_Mo), opt.lambda1(2)*Th(W_R1.*TV_R1));
else
    g2=0;
end

if opt.lambda2(1) ~= 0 || opt.lambda2(2) ~= 0
    u = 1e-6;
    W1_Mo = sqrt(conj(wc_Mo).*wc_Mo + u);
    W1_R1 = sqrt(conj(wc_R1).*wc_R1 + u);
    W1_Mo = 1./W1_Mo;
    W1_R1 = 1./W1_R1;
    g3 = cat(3, opt.lambda2(1)*compWhx(W1_Mo.*wc_Mo,opt,fsize_Mo), ...
        opt.lambda2(2)*compWhx(W1_R1.*wc_R1, opt, fsize_R1));
else
    g3 = 0;
end

grad                    = g1+g2+g3;
grad                    = real(grad(:));