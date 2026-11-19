
function b = b_values_x(x, n, Q)
% Estimate spectral coefficients b_n of x from its AVAR.
%
% INPUT:
%   x : K x N matrix, one signal per row
%   n : vector of nonnegative integer indices
%   Q : optional maximum AVAR averaging factor
%
% OUTPUT:
%   b : K x length(n) matrix
%
% Requires AVAR_values_x.m and AVAR_values_y.m.

    [K, N] = size(x);
    n = n(:);
    M = length(n);

    Qmax = floor((N - 1) / 2);

    if nargin < 3
        Q = Qmax;
    endif

    if Q < 1 || Q > Qmax || Q ~= floor(Q)
        error('Q must be an integer between 1 and %d.', Qmax);
    endif

    if any(n < 0 | n ~= floor(n))
        error('n must contain nonnegative integers.');
    endif

    b = zeros(K, M);

    if isempty(n)
        return;
    endif

    if any(n > Q)
        error('Requested n values must not exceed Q.');
    endif

    % Compute AVAR for q = 1, ..., Q.
    % Process in blocks to avoid allocating a huge AVAR matrix.
    block_size = 256;

    weighted_sum = zeros(K, 1);

    % Store only the AVAR values needed for b_n, n >= 1
    positive_n = unique(n(n >= 1));

    q_needed = [];

    for j = 1:length(positive_n)
        q = positive_n(j);

        while q <= Q
            q_needed(end+1, 1) = q;
            q = 2*q;
        endwhile
    endfor

    q_needed = unique(q_needed);
    AVAR_needed = zeros(K, length(q_needed));

    % Calculate b0 and collect the required AVAR values
    for q_start = 1:block_size:Q

        q_end = min(q_start + block_size - 1, Q);
        q_block = q_start:q_end;

        A = AVAR_values_x(x, q_block);

        % Accumulate q^2 * AVAR(q)
        weighted_sum = weighted_sum ...
                     + A * (q_block(:).^2);

        % Save selected AVAR values
        [tf, loc] = ismember(q_needed, q_block);

        if any(tf)
            AVAR_needed(:,tf) = A(:,loc(tf));
        endif

    endfor

    % Cesaro estimate of b0
    b0 = weighted_sum / (3*Q);

    % Assemble requested coefficients
    for i = 1:M

        ni = n(i);

        if ni == 0
            b(:,i) = b0;
        else
            q = ni;
            dyadic_sum = zeros(K, 1);

            while q <= Q
                idx = find(q_needed == q, 1);
                dyadic_sum = dyadic_sum + AVAR_needed(:,idx);
                q = 2*q;
            endwhile

            b(:,i) = 2*b0 - (ni^2/2)*dyadic_sum;
        endif

    endfor

endfunction

