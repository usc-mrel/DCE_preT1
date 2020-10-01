function [petab, U1, U2] = genSP(yres, zres, N, pts, target_rot, density, calibr, seed, ellip_flag, debug)
%   Cartesian spiral pattern generation. The code is adapted from R. Marc
%   Lebel's Cartesian spiral trajectory generation.
%
%   Authro: Zhibo
%   Date: 04/2020
%
%   Output:         petab       PE table, [N 2]
%                   U1          Pattern inclusing dummy pulses, eddy current correction, calibration region and PE
%                   U2          Pattern consisting of PE only
%   Input:          ky          ky FOV, [scalar]
%                   kz          kz FOV, [scalar]
%                   N           Number of samples, [scalar]
%                   pts         Points per spiral arm, [scalar]
%                   target_rot  Number of target rotation for spiral arm, [scalar]
%                   density     Spiral density factor, [scalar between 0 and 1]
%                   calibr      Central calibration region flag. [0 or 1]
%                                   0, Without calibraion
%                                   1, With caliration
%                   seed        Random seed, [scalr]
%                   ellip_flag  Elliptical footprint flag, [0 1]
%                                   0, Without elliptical footprint
%                                   1, With elliptical footprint
%                   debug       Debug mode flag, [0 1]

if nargin < 10
    debug = 1;
end

if nargin < 9
    ellip_flag = 0;
end

if nargin < 8
    seed = 1;
end

if nargin < 6
    error('Error: Not enough input parameters!')
end

NFULLC                  = 40;
NSTST                   = 512;
NEDDYR                  = 3;
NEDDYE                  = 8;
MGAN                    = 0.618033988749894902525738871191;

fibonacci               = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765];
fib_base                = 0;
% seed                    = 1;
actual_rot              = 0;

diff_rot                = 1000000;
nfully                  = NFULLC;
nfullz                  = NFULLC;

ovs                     = ceil(3*((yres + zres)/2)/pts);
pts                     = ovs * pts;
ptsr                    = pts-1;

for indx = 1:19
    Nact                = ptsr*fibonacci(indx)*floor(N*ovs/(ptsr*fibonacci(indx)) + 1);
    arms                = Nact/ptsr;
    
    slope_rot           = mod(MGAN*arms, 1);
    if (slope_rot > 0.5)
        slope_rot       = 1 - slope_rot;
    end
    
    cur_rot             = slope_rot*ptsr;
    if (abs(cur_rot - target_rot) < diff_rot)
        diff_rot        = abs(cur_rot - target_rot);
        fib_base        = fibonacci(indx);
        actual_rot      = cur_rot;
    end
end

Nact                    = pts*fib_base*floor(N*ovs/(pts*fib_base) + 1);
Nactr                   = ptsr*fib_base*floor(N*ovs/(pts*fib_base) + 1);
arms                    = floor(Nactr/ptsr);
armsneeded              = ceil(2*N/(pts/ovs + 1));

if debug
    fprintf("\nMaking phyllotaxic spiral:\n");
    fprintf("\tFibonacci number: %d\n",fib_base);
    fprintf("\tInput points: %d\n",N);
    fprintf("\tArm oversampling: %d\n",ovs);
    fprintf("\tTotal points computed: %d\n",Nact);
    fprintf("\tTotal points computed (reduced): %d\n",Nactr);
    fprintf("\tPoints per spiral: %d\n",pts);
    fprintf("\tSpiral arms: %d\n",arms);
    fprintf("\tRotations per arm: %f\n",actual_rot);
end

ky                      = zeros(size((armsneeded*(pts/ovs + 1) + nfully*nfullz)));
kz                      = zeros(size((armsneeded*(pts/ovs + 1) + nfully*nfullz)));

sloc                    = 1;

for indx = 1:NSTST
    ky(sloc)            = yres/2 + 1;
    kz(sloc)            = zres/2 + 1;
    sloc                = sloc + 1;
end

for indx = 0:NEDDYR - 1
    for indx2 = 0:NEDDYE - 1
        for indx3 = 0:NEDDYE - 1
            kz(sloc)    = 1;
            ky(sloc)    = yres/2 + 1;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1 - NEDDYE/2 + indx2;
            ky(sloc)    = yres/2 + 1 - NEDDYE/2 + indx3;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1;
            ky(sloc)    = yres/2 + 1;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1;
            ky(sloc)    = yres/2 + 1;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1 - NEDDYE/2 + indx2;
            ky(sloc)    = yres/2 + 1 - NEDDYE/2 + indx3;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1;
            ky(sloc)    = yres/2 + 1;
            sloc        = sloc + 1;
        end
    end
end

for indx = 0:NEDDYR - 1
    for indx2 = 0:NEDDYE - 1
        for indx3 = 0:NEDDYE - 1
            kz(sloc)    = zres/2 + 1;
            ky(sloc)    = 1;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1 - NEDDYE/2 + indx2;
            ky(sloc)    = yres/2 + 1 - NEDDYE/2 + indx3;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1;
            ky(sloc)    = yres/2 + 1;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1;
            ky(sloc)    = yres/2 + 1;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1 - NEDDYE/2 + indx2;
            ky(sloc)    = yres/2 + 1 - NEDDYE/2 + indx3;
            sloc        = sloc + 1;
            
            kz(sloc)    = zres/2 + 1;
            ky(sloc)    = yres/2 + 1;
            sloc        = sloc + 1;
        end
    end
end
if calibr
    sloc_tmp = sloc;
end

for indx = 0:nfully - 1
    kyt = indx  - nfully/2 + 1 + yres/2;
    for indx2 = 0:nfullz - 1
        kzt             = indx2 - nfullz/2 + 1 + zres/2;
        ky(sloc)        = kyt;
        kz(sloc)        = kzt;
        sloc            = sloc + 1;
    end
end
if ~calibr
    sloc_tmp = sloc;
end

for indx = 0:armsneeded -1
    seed                = myrand(seed);
    rnum                = (mod(seed, 1000)) / 1000;
    shft                = floor(rnum*ovs);
    for indx2 = (pts-1):-1:0
        if (indx2 == 0)
            ploc        = 0;
        else
            ploc        = (indx2 - 1)*arms + indx;
        end
        
        cur_phs         = 2*pi*MGAN*ploc;
        cur_rad         = sqrt(2)*power(ploc/Nactr, density);
        
        kyt             = cur_rad * cos(cur_phs);
        kzt             = cur_rad * sin(cur_phs);
        
        kyt             = floor(yres*(kyt/2 + 0.5) + 1);
        kzt             = floor(zres*(kzt/2 + 0.5) + 1);
        
        if (indx2 ==0 || (((mod(indx2 +shft, ovs)) == 0) && kyt >= 1 && kyt <= yres && kzt>= 1 && kzt <= zres))
            if ellip_flag && ((kyt-0.5-yres/2)^2/(0.5+yres/2)^2 + (kzt-0.5-zres/2)^2/(0.5+zres/2)^2 > 1)
                continue;
            end
            ky(sloc)    = floor(kyt);
            kz(sloc)    = floor(kzt);
            
            sloc        = sloc + 1;
        end
    end
end

petab(:, 1)             = ky(1:N);
petab(:, 2)             = kz(1:N);

U1 = zeros(yres, zres);
U2 = zeros(yres, zres);
for ct = 1:N
    U1(ky(ct), kz(ct))   = U1(ky(ct), kz(ct)) + 1;
    if ct >= sloc_tmp
        U2(ky(ct), kz(ct))   = U2(ky(ct), kz(ct)) + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function out = myrand(in)
        out             = mod(in*16807, 2147483647);
    end
end