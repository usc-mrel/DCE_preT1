function k = genKspace(M0, R1, B1, FA, TR, sMaps, FTdim, FTshift)
%   Generate synthesized k-space data.
%
%   Author: Zhibo Zhu
%   Date: 04/2020
%
%   Usage: k = genKspace(Mo, R1, B1, FA, TR, sMaps, FTdim, FTshift)
%
%   Input:
%       - M0            Amount of magnetization, [np nv ns]
%       - R1            The logitudinal relaxation rate, [np nv ns]
%       - B1            B1+ map, [np nv ns]
%       - FA            Flip angles list in degrees, [1 nt]
%       - TR            Repitition time in second, [scalar]
%       - sMaps         Coil sensitivity map, [np nv ns nt nr]
%       - FTdim         Fourier transform dimensions
%       - FTshift       Fourier transform shift flag, [0 1]
%   Output:
%       - k             k-space data, [np nv ns nt nr]
%   Others:
%       - np            Number of pointx along x
%       - nv            Number of points along y
%       - ns            Number of slices
%       - nt            Number of contrast, e.g. time fram, flip angle
%       - nr            Number of coils

img = spgr(M0, R1, B1, FA, TR);
[np, nv, ns, nt] = size(img);
nr = size(sMaps, 5);
imgS = repmat(sMaps, [1 1 ns nt 1]) .* repmat(img, [1 1 1 1 nr]);
k = fFastFT(imgS, FTdim, FTshift);
end