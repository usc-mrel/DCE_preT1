function [mu, sigma] = mean_std(region)

region = region(region~=0);
mu = mean(region);
sigma = std(region);