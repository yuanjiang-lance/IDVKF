function out_Sig = low_filter(Sig,T,cuttime)
%
% FIR low-pass filtering in frequency domain
%
% ------------ Input -------------
%  Sig: measured signal for filtering
%  T: time duration
%  cuttime: cutoff time of the filter
%
% ------------ Output ------------
%  out_Sig: filtered signal
%
% Author: Yuan JIANG
% Time: 2026-06-06

n0 = length(Sig);
n = floor(n0*0.8); % length of the filter
w1 = 2*cuttime/T; % normalized cutoff time bandwidth 

if mod(n,2) == 0
L = n;
else
L = n + 1;
end
b = fir1(L,w1,'low');
Sig1 = conv(b,Sig);
out_Sig = Sig1(L/2+1:L/2+n0); % correct the phase shift of the filter