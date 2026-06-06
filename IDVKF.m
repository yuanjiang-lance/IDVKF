function [Sigest, GDest, IAest, r, Sigiter, GDiter, IAiter] = IDVKF(Sig, T, iniGD, r0, beta, winLen, tolout, maxit)
% 
% Iterative Dispersive Vold-Kalman Filter (IDVKF)
%
% ------------- Input --------------
%  Sig : measured signal, one row/colum vector
%  T : signal time duration (s)
%  iniGD : initial group delays (GDs), each GD lies in onw row
%          ATTENTION: The length of iniGD and Sig must be equal
%  r0 : initial bandwidth controller, larger r0 results in narrower bandwidth
%       one number if all components share the same bendwidth, or a vector
%       whose length is same to the row number of iniGD
%  beta : GD smooth degree controlling parameter, smaller beta results in smoother IF result
%         one number if all components share the same parameter, or a vector
%         whose length is same to the row number of iniGD
%  tolout: iteration stopping criterion
%  maxit: maximum iteration number to avoid dead loop, default 300
%
% ------------- Output --------------
%  Sigest: decomposed dispersive signal components (n_comp * sigLen)
%  GDest: estimated GD (n_comp * sigLen)
%  IAest: estimated envelope (n_comp * sigLen)
%  r: bandwidth controller record (n_comp * sigLen-2 * iter)
%  Sigiter: decomposed dispersive signal components record at each
%           iteration (n_comp * sigLen * iter)
%  GDiter: estimated GD record at each iteration (n_comp * sigLen * iter)
%  IAiter: estimated envelope record at each iteration (n_comp * sigLen * iter)
%
% Author: Yuan JIANG
% Time: 2026-06-06


%% Initialization
if nargin < 7, maxit = 300; end
if size(Sig, 1) > size(Sig, 2)
    Sig = Sig.';
end
if size(iniGD, 1) > size(iniGD, 2), iniGD = iniGD.'; end
if size(Sig, 2) ~= size(iniGD, 2)
    error('Length of signal does not equal the length of initial GD');
end
[num, N] = size(iniGD); % num: number of components; N: length of signal
if isscalar(r0)
    r0 = r0 * ones(1, N-2);
elseif isvector(r0)
    r0 = r0(:) .* ones(num, N-2);
end

if length(beta) == 1
    beta = beta * ones(1, num);
elseif length(beta) ~= num
    error('Length of GD smoothness parameter beta must equal to 1 or the number of components');
end
%%
f = (0: N-1) / T;
e = ones(N, 1);
e2 = -2*e;
Oper = spdiags([e e2 e], 0:2, N-2, N);
Operdoub = Oper' * Oper;
GDiter = zeros(num, N, maxit+1);
GDiter(:, :, 1) = iniGD;
Sigiter = zeros(num, N, maxit);
IAiter = zeros(num, N, maxit);
r = zeros(num, N-2, maxit+1);
r(:, :, 1) = r0;

it = 0;
sDif = tolout + 1;

while (sDif > tolout && it <= maxit)
    
it = it + 1;
PHIdoub = [];
Theta = [];

for i = 1: num
    % initialize phase matrix and operator matrix
    phase = exp(-1j*2*pi*cumtrapz(f, GDiter(i, :, it)));
    Theta = [Theta, spdiags(phase(:), 0, N, N)];
    ROper = spdiags(r(i, :, it)', 0, N-2, N-2) * Oper;
    ROperdoub = ROper' * ROper;
    PHIdoub = blkdiag(PHIdoub, ROperdoub);
end

% calculate envelope vector
a = (PHIdoub + Theta'*Theta) \ (Theta'*Sig(:));

for i = 1: num
    % signal reconstruction
    ai = a((i-1)*N+1: i*N);
    y = Theta(:, (i-1)*N+1: i*N) * ai;
    Sigiter(i, :, it) = y;
    IAiter(i, :, it) = ai;
    
    % update demodulated signals
    deltaphase = unwrap(angle(ai));
    deltaGD = Differ(deltaphase, 1/T) / (2*pi);
    deltaGD = (1/beta(i) * Operdoub + speye(N)) \ deltaGD(:);
    GDiter(i, :, it+1) = GDiter(i, :, it) - deltaGD.';
    for mid = 1: N-2
        left = max(mid - floor(winLen/2) + 1, 1);
        right = min(mid + floor(winLen/2) + 1, N);
        y_win = y(left:right);
        Sig_win = Sig(left:right);
        r(i, mid, it+1) = r(i, mid, it) * real(y_win'*Sig_win(:)) / (y_win'*y_win);
    end
end

% ------------ iteration cretirion --------------
if it > 1
    sDif = 0;
    for i = 1: num
        sDif = sDif + (norm(Sigiter(i, :, it) - Sigiter(i, :, it-1)) / norm(Sigiter(i, :, it-1))) ^2;
    end
end

end

Sigiter = Sigiter(:, :, 1:it);
GDiter = GDiter(:, :, 1:it);
IAiter = IAiter(:, :, 1:it);
r = r(:, :, 1:it);
Sigest = Sigiter(:, :, end);
GDest = GDiter(:, :, end);
IAest = IAiter(:, :, end);
