function S = spgr(M0,R1,B1,FA,TR)
%   Generate synthesized VFA image data.
%
%   Author: Zhibo Zhu
%   Date: 04/2020
%
%   Usage: S = spgr(M0,R1,B1,FA,TR)
%
%   Input:
%       - M0            Amount of magnetization, [np nv ns]
%       - R1            The logitudinal relaxation rate, [np nv ns]
%       - B1            B1+ map, [np nv ns]
%       - FA            Flip angles list in degrees, [1 nt]
%       - TR            Repitition time in second, [scalar]
%   Output:
%       - S             VFA image data, [np nv ns nt]
%   Others:
%       - np            Number of pointx along x
%       - nv            Number of points along y
%       - ns            Number of slices
%       - nt            Number of contrast, e.g. time fram, flip angle

E1                      = exp(-TR.*R1);
FA                      = FA*pi/180;
S                       = zeros([size(M0) length(FA)]);
for nt = 1:length(FA)
    S(:,:,:,nt)         = M0.*sin(B1.*FA(nt)).*(1-E1)./(1-E1.*cos(B1.*FA(nt)));
end
end