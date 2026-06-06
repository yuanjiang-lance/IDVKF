%
% Add Gaussian white noise into signal
%
% Author: Yuan JIANG
% Time: 2026-06-06

clear, clc, close all;
load data_origin;
seed = 42;
rng(seed);
%%
SNR = 5;
SigFn = awgn(SigF, SNR, 'measured');
SigFn_extend = [SigFn, conj(fliplr(SigFn(2: ceil(Nt/2))))];
SigTn = real(ifft(SigFn_extend));

sigshow(SigTn, Fs);     % noisy signal plot
sigshow(SigT, Fs);      % non-noise signal plot

window = 128;
Spec = mySFFT(SigFn(:), Fs, Nt, window);
stftshow(t, f, Spec);   % noisy signal TF spectrum

save(['data_' num2str(SNR) 'dB'], ...
    'GD', 'IA', 'SigF', 'SigT', 'SigF_mode', 'SigT_mode', ...
    'Fs', 'T', 'Nt', 'Nf', 't', 'f', 'SigFn', 'SigTn');