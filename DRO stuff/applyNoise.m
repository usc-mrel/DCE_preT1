function k = applyNoise_new(k, img_precontrast, noisecov, SNR)
if isinf(SNR)
    return;
end
[np, nv, ns, nt, nr] = size(k);
currentSig = max(img_precontrast(:));
currentSNR = currentSig / mean(sqrt(diag(noisecov)));
SNR = currentSNR;
scaleNoise = (currentSNR/SNR)^2 / 2 * (np*nv); % 2D kspace noise
noisecov = noisecov * scaleNoise;

noise = mvnrnd(zeros(1, nr), noisecov, np*nv*ns*nt ) + 1i * mvnrnd(zeros(1, nr), noisecov, np*nv*ns*nt);
noise = reshape(noise, [np nv ns nt nr]);

img_noise = iFastFT(noise, [1 2], 1);
tmp = reshape(img_noise, [], 8);
cov_tmp = cov(tmp);
snr = currentSig/mean(sqrt(diag(cov_tmp)));

k = k + noise;