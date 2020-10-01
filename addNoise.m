function k = addNoise(k, img_precontrast, noisecov, SNR)

if isinf(SNR) % Noiseless
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

k = k + noise;