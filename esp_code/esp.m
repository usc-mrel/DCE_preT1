function [esp_plot, grid, lambda] = esp(ref, img, voxel_spacing, lambda, center_corr)
% ESP Computes the Fourier Radial Error Spectrum Plot (ESP) as described in [1,2]
%
% [1] T. H. Kim, J. P. Haldar. The Fourier Radial Error Spectrum Plot: A
%      more nuanced quantitative evaluation of image reconstruction
%      quality. IEEE International Symposium on Biomedical Imaging,
%      Washington, DC, 2018.
%
% [2] T. H. Kim, J. P. Haldar. Assessing MR image reconstruction quality
%      using the Fourier Radial Error Spectrum plot. Joint Annual Meeting
%      ISMRM-ESMRMB, Paris, 2018.
%
% This software is available from <a href="matlab:web('http://mr.usc.edu/download/ESP/')">http://mr.usc.edu/download/ESP/</a>.
% As described on that page, use of this software (or its derivatives) in
% your own work requires that you cite the references listed at
% <a href="matlab:web('http://mr.usc.edu/download/ESP/')">http://mr.usc.edu/download/ESP/</a>.
%
% Usage:
% [esp_plot, grid, lambda] = ESP(ref, img, voxel_spacing, lambda, k-space_center_correction);
%
% Inputs:
%
%    *ref:  A 1D, 2D, or 3D gold-standard reference image of size Nz x Ny x Nz
%            (singleton dimensions are allowed)
%
%    *img: An image with the same dimensions as ref that will be evaluated
%             for its accuracy with respect to the gold standard reference image
%
%    *voxel_spacing:  Optional input that provides the voxel spacing of the
%                             image along each dimension (in units of mm).
%                             Example: use [2 1] for a 2D image with 2mm
%                             voxel spacing along the first dimension and
%                             1mm along the second dimension.  This is used
%                             to place appropriate unit labels on the ESP.
%                             If not specified, then the ESP will be
%                             displayed in generic 'Nyquist Units'.
%
%    *lambda: Optional input that specifies the value of the regularization
%                  parameter to use for the smoothing spline fit.  For
%                  example, lambda = 1e-7.  You can also set lambda = 'cv'
%                  if you would like the code to use cross-validation to
%                  automatically identify an appropriate regularization
%                  parameter value (which takes longer than if a numeric
%                  value is provided).  If not provided by the user, the
%                  code will use lambda = 'cv' by default.
%
%     *center_corr: Optional input parameter that takes value equal to
%                         either true or false.  If true, the code will
%                         identify the center of k-space based on the
%                         sample of the Fourier transform of the reference
%                         image that has the largest magnitude. If false,
%                         the code will identify the center of k-space
%                         based on standard Fourier transform conventions.
%                         If not provided by the user, the code will
%                         default to center_corr = true.  This option may
%                         be useful in cases where the center of k-space
%                         might be positioned incorrectly in the k-space
%                         matrix, e.g., due to timing delays or linear
%                         phase properties of the image.
%
% Outputs:
%
%     *If ESP(...) is called with no output arguments, then it will
%     generate a plot of the ESP in the current figure window
%
%     *[esp_plot, grid] = ESP(...) will return the y-axis values of the ESP
%     plot in the vector esp_plot, with the corresponding x-axis values
%     (describing spatial frequency information) in the vector grid.  This
%     allows you to plot the ESP yourself, e.g., by calling
%     plot(grid, esp_plot)
%
%     *[esp_plot, grid, lambda] = ESP(...) will additionally output the
%     value of the regularization parameter that was used for the smoothing
%     spline fit.  This is primarily useful if you will be creating
%     multiple ESPs and want to use the same regularization parameter for
%     all of them.
%
%     Written by: Tae Hyung Kim (taehyung@usc.edu) 04/03/2018

dim = ndims(ref);   % dimension of the image
num_basis = 250;    % number of basis for smoothing spline

% Error check
if dim > 3
    error('Error: this function only supports maximum 3D data.');
end

if dim ~= ndims(img)
    error('Error: The dimension of the reference and the image should match.');
elseif sum(~(size(ref)==size(img)))
    error('Error: the size of the reference and the image should match.');
end

if exist('voxel_spacing') && numel(voxel_spacing)>1 && numel(voxel_spacing) < dim
    error('Error: dimensions of the image and the number of voxel_spacing entries should match.');
end

if exist('voxel_spacing') && numel(voxel_spacing)==0
    clear voxel_spacing;
end

if not(exist('lambda')) || isempty(lambda)
    lambda = 'cv';
end

if isnumeric(lambda)
    loocv = false;
elseif strcmpi(lambda,'cv')
    loocv = true;
else
    error('Error: lambda input has an unexpected value');
end

if not(exist('center_corr')) || isempty(center_corr)
    center_corr = true;
elseif not(logical('center_corr'))
    error('Error: center_corr input has an unexpected value');
end

% k-space spacing (for Fourier radius)
[nx ny nz] = size(ref);
if exist('voxel_spacing')
    k_spacing = 1./(voxel_spacing .* size(ref));
    nyquist_units = false;
else
    k_spacing = ones(1,dim);
    nyquist_units = true;
end

%% Fourier representation of reference
ideal = ftn(ref);
[~,idx] = max(vect(ideal));
[i, j, k] = ind2sub([nx ny nz], idx);

% for k-space center correction
if center_corr
    delx = (floor(nx/2)+1) - i;
    dely = (floor(ny/2)+1) - j;
    delz = (floor(nz/2)+1) - k;
else
    delx = 0;
    dely = 0;
    delz = 0;
end
y_ideal = vect(abs(ideal).^2);

%% Fourier representation of error
recon = ftn(img);
y_err = vect(abs(recon - ideal).^2);

y = [y_ideal y_err];

%% k-space grid and Fourier radius
[y_grid,x_grid,z_grid] = meshgrid(-floor(ny/2):floor(ny/2) - ~rem(ny,2), -floor(nx/2):floor(nx/2) - ~rem(nx,2), -floor(nz/2):floor(nz/2) - ~rem(nz,2));
x_grid = (x_grid + delx)*(k_spacing(1));
if dim >=2
    y_grid  = (y_grid + dely)*(k_spacing(2));
    if dim >=3
        z_grid = (z_grid + delz)*(k_spacing(3));
    end
end
r = vect(sqrt(x_grid.^2 + y_grid.^2 + z_grid.^2)); % Fourier radius
r = floor(r*1e10)/1e10;                            % rounding radius values to prevent numerical issues

%% define grid for spline fitting
grid = linspace(0, max(r)*0.85, 100).';

%% Combine data with the the same radius (reduce the matrix size for spline fitting)
r_uniq = unique(r(:));      % unique radius values
repeated_r_count = hist(r, r_uniq).';   % combine samples with the same radius
num_r_uniq = numel(r_uniq);

[~, idx] = sort(r);
y_sort = y(idx,:);
y_reduced = zeros(num_r_uniq, size(y,2));

% combine k-space samples with the same radius
curr_idx = 1;
for k = 1 : num_r_uniq
    next_idx = curr_idx + repeated_r_count(k);
    y_reduced(k,:) = sum(y_sort(curr_idx:next_idx-1,:),1);
    curr_idx = next_idx;
end
y_reduced = y_reduced ./ repmat(sqrt(repeated_r_count),[1 size(y_reduced,2)]);

y_ideal = y_reduced(:,1);
y_err = y_reduced(:,2);

% data to fit
d = [y_reduced; zeros(num_basis, size(y_reduced,2))];

%% Defind forward model

% We want to find the solution of min_x ||Dx - y||_2^2 + lambda||Lx||_2^2,
% where x is coefficients for the spline basis functions
% by x = (D'*D + lambda*omega)^-1 * D'y   (omega = L'L)

% location of basis functions
knots = linspace(r_uniq(1),r_uniq(end),num_basis)';

% define basis function / forward model matrix
D = zeros(num_r_uniq, num_basis);
D(:,1) = ones(size(r_uniq));
D(:,2) = r_uniq;
for k = 1:num_basis-2
    D(:,k+2) = ( subplus(r_uniq-knots(k)).^3 - subplus(r_uniq-knots(num_basis)).^3) ./ (knots(num_basis) - knots(k)) ...
        -  ( subplus(r_uniq-knots(num_basis-1)).^3 - subplus(r_uniq-knots(num_basis)).^3) ./ (knots(num_basis) - knots(num_basis-1)) ;
end
D = D .* repmat(sqrt(repeated_r_count), [1 num_basis]); % scaling corresponding to the number of samples with the same radius

%% Regularization (for smoothness)
row_num = repmat((1:(num_basis-2))', [1 num_basis-2]);
col_num = repmat(1:(num_basis-2), [num_basis-2 1]);
i = min(row_num,col_num);
j = max(row_num,col_num);
ei = knots(i);
ej = knots(j);
ek = knots(num_basis);
ek_1 = knots(num_basis-1);

omega = zeros(num_basis,num_basis);
omega(3:end,3:end) = 6*(ej-ek).*(3*ei-ej-2*ek)./(ek-ei) - 6*(ek_1-ek).*(3*ei-ek_1-2*ek)./(ek-ei) - 6*(ek_1-ek).*(3*ej-ek_1-2*ek)./(ek-ej) + 12*(ek-ek_1);

[V, E] = eig(omega);
% L: Tikhonov matrix
L = sqrt(E)*V';

%% Leave one out cross validation (LOOCV) if necessary
if loocv
    lambda_candidate = logspace(-12,-3,10); % lambda candadates
    loocv_ideal = zeros(size(lambda_candidate));    % MSE estimate of reference image fitting
    loocv_err = zeros(size(lambda_candidate));      % MSE esitmate of error image fitting
    
    for k = 1:numel(lambda_candidate)
        A =[D; sqrt(lambda_candidate(k))*L] ;
        pinv_A = pinv(A);
        y_hat = D*(pinv_A*d);
        y_ideal_hat = y_hat(:,1);
        y_err_hat = y_hat(:,2);
        
        diag_S = sum(D.*pinv_A(:,1:end-num_basis).',2);
        Q = numel(y_ideal_hat);
        loocv_ideal(k) = sum(((y_ideal-y_ideal_hat)./(1- diag_S)).^2)/Q;    % MSE estimate of reference image fitting
        loocv_err(k) = sum(((y_err-y_err_hat)./(1- diag_S)).^2)/Q;          % MSE esitmate of error image fitting
    end
    [~, min_idx_ideal] = min(loocv_ideal);      % find a lambda with minimum MSE estimate (for reference fitting)
    [~, min_idx_err] = min(loocv_err);          % find a lambda with minimum MSE estimate (for error image fitting)
    lambda_idx = min(min_idx_ideal,min_idx_err);
    lambda = lambda_candidate(lambda_idx);
end

%% Smoothing spline

% Forward model with regularizer
A = [D; sqrt(lambda)*L];

% Solution: coefficients for spline basis functions
theta = A\d;

% spline fitting
fr = 0;
fr = fr + ones(size(grid))*theta(1,:) + grid*theta(2,:);
for k = 1:num_basis-2
    fr = fr + ((subplus(grid-knots(k)).^3 - subplus(grid-knots(num_basis)).^3) ./ (knots(num_basis) - knots(k)) ...
        - (subplus(grid-knots(num_basis-1)).^3 - subplus(grid-knots(num_basis)).^3) ./ (knots(num_basis) - knots(num_basis-1)))*theta(k+2,:) ;
end

f_ideal = fr(:,1);  % spline fitting for the reference
f_err = fr(:,2);    % spline fitting for the error image

%% ESP plot
esp_plot = real(sqrt(f_err./f_ideal));
if not(nargout)
    if nyquist_units
        plot(grid, esp_plot,'linewidth',3);xlim([0 max(grid)]); xlabel('Spatial Frequency (Nyquist units)'); ylabel('Relative Error'); title('ESP');  % hold on
    else
        plot(grid, esp_plot,'linewidth',3);xlim([0 max(grid)]); xlabel('Spatial Frequency (mm^{-1})'); ylabel('Relative Error'); title('ESP');  % hold on
    end
end
end

%%
function [ out ] = ftn( in )
out = fftshift(fftn(ifftshift(in)));
end

%%
function [ out ] = vect( in )
out = in(:);
end