%% Demo script for Marc's Cartesian spiral
%% It generates exactly the same pattern as Marc's Sparse DCE
zy_views = 240;
zy_slices = 120;
SPviews = zy_views * zy_slices * 4;
SPpnts = 1/2 * (zy_slices + zy_views)/4;
SProtations = 0.333;
SPdensity = 0.65;
[petab, U1, U2] = genSP(zy_views, zy_slices, SPviews, SPpnts, SProtations, SPdensity, 1, 1);