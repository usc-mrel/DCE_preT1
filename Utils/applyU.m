function [kU, U] = applyU(k, pattern, R, realization)
[np, nv, ns, nt, nr] = size(k);
if R ~= 1
    if pattern == 0 || pattern == 1 % Cartesian spiral
        U = makeSP();
    else % RGA
        U = makeRGA();
    end
else
    U = 1;
end

kU = k .* U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function U = makeSP()
        U = zeros(np, nv, ns, nt, nr);
        kfov = zeros(nv, ns);
        pts = 1/2 * (np + nv) / 4;
        rot = 0.333;
        density = 0.65;
        N = np * nv * nt/R;
        A = round(N / 10);
        B = 4 * A;
        if pattern == 0 % Cartesian spiral rect
            ellip_flag = 0;
            kfov = ones(np, nv);
        else % Cartesian spiral elliptic
            ellip_flag = 1;
            for m = 1:np
                for n = 1:nv
                    if ((m-0.5-np/2)^2/(0.5+np/2)^2 + (n-0.5-nv/2)^2/(0.5+nv/2)^2 <= 1)
                        kfov(m, n) = 1;
                    end
                end
            end
        end
        for ct = 1:7
            if ct == 1
                [~, ~, U1] = genSP(np, nv, 4417+A, pts, rot, density, 0, ct+(realization-1)*10, ellip_flag, 0);
                U1 = single(logical(U1));
                U1(~kfov) = 0.1;
                U(:, :, :, ct, :) = repmat(U1, [1 1 ns 1 nr]);
            elseif ct ~= 7
                [~, ~, U1] = genSP(np, nv, 4417+A, pts, rot, density, 0, ct+(realization-1)*10, ellip_flag, 0);
                U1 = single(logical(U1));
                U1(~kfov) = 0.1;
                U(:, :, :, ct, :) = repmat(U1, [1 1 ns 1 nr]);
            else
                [~, ~, U1] = genSP(np, nv, 4417+B, pts, rot, density, 0, ct+(realization-1)*10, ellip_flag, 0);
                U1 = single(logical(U1));
                U1(~kfov) = 0.1;
                U(:, :, :, ct, :) = repmat(U1, [1 1 ns 1 nr]);
            end
            clear U1;
        end
    end

    function U = makeRGA()
        U = zeros(np, nv, ns, nt, nr);
        kfov = zeros(np, nv);
        if pattern == 2 % RGA elliptic
            bincode = '0011';
            for m = 1:np
                for n = 1:nv
                    if ((m-0.5-np/2)^2/(0.5+np/2)^2 + (n-0.5-nv/2)^2/(0.5+nv/2)^2 <= 1)
                        kfov(m, n) = 1;
                    end
                end
            end
        elseif pattern == 3 % RGA rect
            bincode = '0010';
            kfov = ones(np, nv);
        end
        N = np * nv * nt / R;
        A = round(N / 10);
        B = 4 * A;
        for ct = 1:7
            if ct ~= 7
                [~, U1] = genRGA(np, nv, np, nv, A, bin2dec(bincode), 0, 1, 0.3, 0, 0, 0, 0.1*(ct-1), 0, (realization-1)*pi/5);
                U1 = single(logical(U1));
                U1(~kfov) = 0.1;
                U(:, :, :, ct, :) = repmat(U1, [1 1 ns 1 nr]);
            else
                [~, U1] = genRGA(np, nv, np, nv, B, bin2dec(bincode), 0, 1, 0.3, 0, 0, 0, 0.1*(ct-1), 0, (realization-1)*pi/5);
                U1 = single(logical(U1));
                U1(~kfov) = 0.1;
                U(:, :, :, ct, :) = repmat(U1, [1 1 ns 1 nr]);
            end
        end
        clear U1;
    end
end