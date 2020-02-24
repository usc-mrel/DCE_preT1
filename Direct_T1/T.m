function v = T(u)
% Help function computing 2D forward finite difference
% u:    np-nv-nu
% TV:   np-nv-2-nu

[nx, ny, nu] = size(u);
v = zeros(nx, ny, 2, nu);

%%
v(1:end-1, :, 1, :) = u(2:end, :, :) - u(1:end-1, :, :);

%%
v(:, 1:end-1, 2, :) = u(:, 2:end, :) - u(:, 1:end-1, :);