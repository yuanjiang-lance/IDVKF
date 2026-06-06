%
% Generate a numerical 3-component dispersive signal (no noise)
%
% Author: Yuan JIANG
% Time: 2026-06-06

clear, clc, close all;
%%
Fs = 2000;
T = 10;
Nt = T * Fs;
Nf = floor(Nt/2) + 1;
t = (0: Nt-1)/Fs;
f = (0: Nf-1)/T;
%%
GD1 = -1.6e-5*f.^2 + 0.015*f + 2;
GD2 = -3e-8*f.^3 + 5.8e-5*f.^2 - 0.025*f + 5;
GD3 = 1.25e-8*f.^3 - 1.3e-5*f.^2 - 0.003*f + 8;
GD = [GD1; GD2; GD3];

phase1 = exp(-1j*2*pi * (-1.6e-5*f.^3/3 + 0.015*f.^2/2 + 2*f + 0.2));
phase2 = exp(-1j*2*pi * (-3e-8*f.^4/4 + 5.8e-5*f.^3/3 - 0.025*f.^2/2 + 5*f + 0.7));
phase3 = exp(-1j*2*pi * (1.25e-8*f.^4/4 - 1.3e-5*f.^3/3 - 0.003*f.^2/2 + 8*f + 1.3));

IA1 = 1 + 5 * normpdf(f/20, 15, 2);
IA2 = 1 + 0.3 * sin(2*pi*f.^2/20000);
IA3 = 0.9 * exp(5e-4*f);
IA = [IA1; IA2; IA3];

SigF1 = IA1 .* phase1;
SigF2 = IA2 .* phase2;
SigF3 = IA3 .* phase3;

SigF1_extend = [SigF1, conj(fliplr(SigF1(2: ceil(Nt/2))))];
SigF2_extend = [SigF2, conj(fliplr(SigF2(2: ceil(Nt/2))))];
SigF3_extend = [SigF3, conj(fliplr(SigF3(2: ceil(Nt/2))))];

SigT1 = ifft(SigF1_extend);
SigT2 = ifft(SigF2_extend);
SigT3 = ifft(SigF3_extend);

SigF_mode = [SigF1; SigF2; SigF3];
SigT_mode = real([SigT1; SigT2; SigT3]);
SigT = real(SigT1 + SigT2 + SigT3);
SigF = SigF1 + SigF2 + SigF3;

sigshow(SigT, Fs);  % signal plot

figure(2)           % GD plot
plot(GD1,f,GD2,f,GD3,f);
axis([0 10 0 1000])

figure(3)           % signal component plot
subplot(311),plot(f,real(SigF1));
subplot(312),plot(f,real(SigF2));
subplot(313),plot(f,real(SigF3));

window = 128;
Spec = mySFFT(SigF(:), Fs, Nt, window);
stftshow(t, f, Spec);   % TF spectrum

save data_origin.mat GD IA SigF SigT SigF_mode SigT_mode Fs T Nt Nf t f;