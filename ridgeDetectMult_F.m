function [tindex, ridgeMag] = ridgeDetectMult_F(SigF, Fs, delta, beta, bw, Nt, WinLen, tolout)
% 
% Extract ridges for multi-component dispersive signals without knowing the component number.
% 
% -------- Input --------
%  SigF : Frequency domain signal
%  Fs : Sample Frequency
%  delta : maximum allowable time variation between two consecutive points for ridge detection
%  beta : controls the smooth degree; the curves will be smoother if beta is smaller
%  bw : the time bandwidth of the TF filter (unit/sec);
%  Nt : Signal length in time domain
%  WinLen : Window length
%
% ------- Output --------
%  tindex : obtained ridge index
%  ridgeMag : the corresponding ridge magnitude
%
% Author: Yuan JIANG
% Time: 2024-05-30

LoopFlag = 1;
tindex = []; ridgeMag = [];
SigRes = SigF;
Nf = length(SigF);
f = (0: Nf-1) / (Nf-1) * Fs/2;
T = 2 * Nf / Fs;

while LoopFlag

% ---------- extract ridge -------------
[Spec, ~, t] = mySFFT(SigRes(:), Fs, Nt, WinLen);
index = DFRE_F(Spec, delta);
GD = IFsmooth(t(index), beta);
phase = cumtrapz(f, GD);

% ---------- signal extraction ----------
dSig = SigRes .* exp(1j * 2*pi * phase);
extractSig = low_filter(dSig, T, bw/2);
extractSig = extractSig .* exp(-1j * 2*pi * phase);

tindexSimple = zeros(1, length(SigF));
ridgeMagSimple = zeros(1, length(SigF));
for j = 1: length(SigF)
    [~, tindexSimple(j)] = min(abs(t - GD(1,j)));
    ridgeMagSimple(j) = abs(Spec(j, tindexSimple(j)));
end

% ------- iteration criterion --------
if (norm(extractSig) / norm(SigF))^2 < tolout
    LoopFlag = 0;
else
    tindex = [tindex; tindexSimple];
    ridgeMag = [ridgeMag; ridgeMagSimple];
    SigRes = SigRes - extractSig;
end

end

end

function outSig = low_filter(Sig, T, cuttime)
%
% FIR low-pass filtering in frequency domain
%
% ----------- Input -------------
%  Sig: measured signal for filtering
%  T: time duration
%  cuttime: cutoff time of the filter
%
% ----------- Output -------------
%  outSig: output signal

n0 = length(Sig);
n = floor(n0 * 0.8);    % length of filter
w = 2 * cuttime / T;    % normalized cutoff time bandwidth

if mod(n, 2) == 0
    L = n;
else
    L = n + 1;
end
b = fir1(L, w, 'low');
Sig = conv(b, Sig);
outSig = Sig(L/2 + 1: L/2 + n0);    %correct the phase shift of the filter
end