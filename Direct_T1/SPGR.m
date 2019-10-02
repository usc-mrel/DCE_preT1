function S = SPGR(Mo,R1,B1,FA,TR)

E1 = exp(-TR.*R1);
FA = FA*pi/180;
S = zeros([size(Mo) 1 length(FA)]);
for nt = 1:length(FA)
    S(:,:,:,nt) = Mo.*sin(B1.*FA(nt)).*(1-E1)./(1-E1.*cos(B1.*FA(nt)));
end