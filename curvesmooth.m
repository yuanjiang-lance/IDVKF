function outf = curvesmooth(f, beta)
%
% Non-parameterized IF/GD smoothing (low-pass filter)
%
% ----------- Input -----------
%  f: input instantaneous frequencies (IFs) or group delays (GDs)
%  beta: controls the smooth degree; the curves will be smoother if beta is smaller
%
% ----------- Output ----------
%  outf: output smoothed IF or GD curves
% 
% Author: Yuan JIANG
% Time: 2026-06-06


[K,N] = size(f); % K is the number of the components£¬N is thenumber of the samples
e = ones(N,1);
e2 = -2*e;
oper = spdiags([e e2 e], 0:2, N-2, N); % the modified second-order difference matrix
opedoub = oper'*oper;%
outf = zeros (K,N);
for i = 1:K
    outf(i,:) = (2/beta*opedoub + speye(N))\f(i,:).'; % smooth the instantaneous frequency curves by low pass filtering
end
