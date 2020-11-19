function k = addNoise(k, img, noisecov, SNR, field_flag)
%   Add synthesized noise to k-space data.
%
%   Author: Zhibo Zhu
%   Date: 04/2020
%
%   Usage: k = addNoise(k, img_precontrast, noisecov, SNR)
%
%   Input:
%       - k             Noiseless k-space data, [np nv ns nt nr]
%       - img           Refernce signal, [Column vector]
%       - noisecov      Noise covariance matrix, [nr nr]
%       - SNR           Desired SNR level, [scalar]
%       - field_flag    Binary flag for field strength, [scalar]
%   Output:
%       - k             Noisy k-space data, [np nv ns nt nr]
%   Others:
%       - np            Number of pointx along x
%       - nv            Number of points along y
%       - ns            Number of slices
%       - nt            Number of contrast, e.g. time fram, flip angle
%       - nr            Number of coils

if isinf(SNR) % Noiseless
    return;
end

[np, nv, ns, nt, nr] = size(k);
if field_flag == 0 % noisecov is not 3T
    currentSig = max(img(:));
    currentSNR = currentSig / mean(sqrt(diag(noisecov)));
    scaleNoise = (currentSNR/SNR)^2 / 2 * (np*nv); % 2D kspace noise
else % noisecov is 3T
    scaleNoise = 1;
end
noisecov = noisecov * scaleNoise;

noise = mvnrnd(zeros(1, nr), noisecov, np*nv*ns*nt ) + 1i * mvnrnd(zeros(1, nr), noisecov, np*nv*ns*nt);
noise = reshape(noise, [np nv ns nt nr]);

k = k + noise;