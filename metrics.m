%% mean and std
[mu_ref, sigma_ref] = mean_std(ref_region);
[mu_test, sigma_test] = mean_std(test_region);

%% MSE
mse = immse(test_region, ref_region);

%% SSIM
ssimval = ssim(test_region, ref_region);