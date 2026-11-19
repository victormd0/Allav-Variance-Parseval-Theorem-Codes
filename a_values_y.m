
function a = a_values_y(y, n)
% Estimate spectral coefficients a_n from the AVAR of y.
%
% INPUT:
%   y : K x N matrix (one signal per row)
%   n : vector of nonnegative integer indices
%
% OUTPUT:
%   a : K x length(n) matrix
%       a(j,i) estimates a_{n(i)} for signal j.
%
% Requires AVAR_values_y.m.

    [K, N] = size(y);
    n = n(:);
    M = length(n);

    qmax = floor(N / 2);

    if any(n < 0 | n ~= floor(n))
        error('n must contain nonnegative integers.');
    endif

    if isempty(n)
        a = zeros(K, 0);
        return;
    endif

    % Determine which S_m are needed
    bases = unique([1; n(n >= 2)-1; ...
                   n(n >= 2); n(n >= 1)+1]);

    if max(bases) > qmax
        error('N is too small for the requested n values.');
    endif

    B = length(bases);

    % Generate all dyadic scales q = 2^ell * m
    q_all = [];
    base_id = [];

    for j = 1:B
        q = bases(j);

        while q <= qmax
            q_all(end+1, 1) = q;
            base_id(end+1, 1) = j;

            q = 2*q;
        endwhile
    endfor

    % Remove duplicate q values
    [q_unique, ~, q_index] = unique(q_all);

    % Compute each AVAR only once
    AVAR = AVAR_values_y(y, q_unique);

    % W maps AVAR(q) into S_m = sum_l AVAR(2^l*m)
    W = sparse(q_index, base_id, 1, ...
               length(q_unique), B);

    S = AVAR * W;

    % Assemble the requested coefficients
    a = zeros(K, M);

    for i = 1:M
        ni = n(i);

        if ni == 0
            j1 = find(bases == 1);
            a(:,i) = 0.5 * S(:,j1);

        elseif ni == 1
            j1 = find(bases == 1);
            j2 = find(bases == 2);

            a(:,i) = 2*S(:,j2) - S(:,j1);

        else
            jm = find(bases == ni-1);
            j0 = find(bases == ni);
            jp = find(bases == ni+1);

            a(:,i) = 0.5*(ni+1)^2*S(:,jp) ...
                    - ni^2*S(:,j0) ...
                    + 0.5*(ni-1)^2*S(:,jm);
        endif
    endfor

endfunction
