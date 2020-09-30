function k = genKspace(Mo, R1, B1, FA, TR, sMaps, FTdim, FTshift)
img = spgr(Mo, R1, B1, FA, TR);
[np, nv, ns, nt] = size(img);
nr = size(sMaps, 5);
imgS = repmat(sMaps, [1 1 ns nt 1]) .* repmat(img, [1 1 1 1 nr]);
k = fFastFT(imgS, FTdim, FTshift);
end