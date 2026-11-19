function y = power_law(K, N, alpha, h, mu)
% Generate K independent Gaussian noise sequences of length N.
%
% Two-sided PSD:
%       S(f) = h * |f|^alpha
%
% Frequency is normalized (cycles/sample).
%
% INPUT:
%   K     : number of realizations
%   N     : samples per realization
%   alpha : spectral exponent (default 0)
%           -2 : 1/f^2 noise
%           -1 : 1/f noise
%            0 : white noise
%            1 : f noise
%            2 : f^2 noise
%   h     : PSD coefficient (default 1)
%   mu    : mean (default 0)
%
% OUTPUT:
%   y : K x N matrix, one realization per row

    if nargin < 3
        alpha = 0;
    endif

    if nargin < 4
        h = 1;
    endif

    if nargin < 5
        mu = 0;
    endif

    % Validate inputs
    if ~ismember(alpha, [-2, -1, 0, 1, 2])
        error('alpha must be -2, -1, 0, 1, or 2.');
    endif

    if h < 0
        error('h must be nonnegative.');
    endif

    % Generate Gaussian white noise
    w = randn(K, N);

    % Fourier transform along each row
    W = fft(w, [], 2);

    % Two-sided normalized frequency grid
    f = 0:N-1;
    f = min(f, N-f) / N;

    % Spectral amplitude multiplier
    H = zeros(1, N);

    nonzero = (f > 0);
    H(nonzero) = sqrt(h) * f(nonzero).^(alpha/2);

    % DC component
    if alpha == 0
        H(1) = sqrt(h);
    else
        H(1) = 0;
    endif

    % Apply spectral shaping
    Y = W .* H;

    % Inverse FFT
    y = mu + real(ifft(Y, [], 2));

endfunction

