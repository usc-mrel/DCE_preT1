function u = Th(v)
% Helping function computing 2D forward finite difference adjoint
% v:    np-nv-2-nu
% u:    np-nv-nu

[nx, ny, ~, nu] = size(v);
u = zeros(nx, ny, nu);

%%
u(1, :, :) = squeeze(v(1, :, 1, :));
u(end, :, :) = squeeze(-v(end-1, :, 1, :));
u(2:end-1, :, :) = squeeze(v(2:end-1, :, 1, :) - v(1:end-2, :, 1, :));

%%
u(:, 1, :) = squeeze(u(:, 1, :)) + squeeze(v(:, 1, 2, :));
u(:, end, :) = squeeze(u(:, end, :)) - squeeze(v(:, end-1, 2, :));
u(:, 2:end-1, :) = squeeze(u(:, 2:end-1, :)) + squeeze((v(:, 2:end-1, 2, :) - v(:, 1:end-2, 2, :)));