clear; clc;

%%
nx      = 240;
ny      = 120;
nt      = 7;
A       = [630; 720; 840; 1008; 1120; 1260; 1440; 1680; 2016; 2240; 2520; 2880; 3360; 4032; 5040; 6720];
B       = [2520; 2880; 3360; 4032; 4480; 5040; 5760; 6720; 8064; 8960; 10080; 11520; 13440; 16128; 20160; 26880];
R       = 28800*7./(6*A + B);


%% Radial spokes reach the boundaries of rectangular
Uyz1    = zeros(nx, ny, 1, nt, length(R));

for ct = 1:length(R)
    for int = 1:6
        [~, U]                  = genRGA(nx, ny, nx, ny, A(ct), bin2dec('0000'), 0, 1, 0.3, 0, 0, 0, 0.1*(int-1));
        Uyz1(:, :, :, int, ct)  = U;
    end
    [~, U]                      = genRGA(nx, ny, nx, ny, B(ct), bin2dec('0000'), 0, 1, 0.3, 0, 0, 0, 0.6);
    Uyz1(:, :, :, end, ct)      = U;
end

save('Uyz1.mat', 'Uyz1', 'R');

%% Radial spokes reach the ellipse bounded by rectangular
Uyz2    = zeros(nx, ny, 1, nt, length(R));

for ct = 1:length(R)
    for int = 1:6
        [~, U]                  = genRGA(nx, ny, nx, ny, A(ct), bin2dec('0001'), 0, 1, 0.3, 0, 0, 0, 0.1*(int-1));
        Uyz2(:, :, :, int, ct)  = U;
    end
    [~, U]                      = genRGA(nx, ny, nx, ny, B(ct), bin2dec('0001'), 0, 1, 0.3, 0, 0, 0, 0.6);
    Uyz2(:, :, :, end, ct)      = U;
end

save('Uyz2.mat', 'Uyz2', 'R');