function AVAR_values_y = AVAR_values_y(y, q)
# Take the matrix y and calculate the AVAR values of each row
# for q(1), q(2), ... q(M).
# Each row of y is treated as an independent signal.
# Return a K x M matrix.

    [K, N] = size(y);
    M = length(q);

    AVAR_values_y = zeros(K, M);

    % Cumulative sum along each row
    S = [zeros(K,1), cumsum(y, 2)];

    for i = 1:M

        qi = q(i);

        % Check if q is valid
        if qi <= 0
            error("Os valores de q devem ser positivos.");
        endif

        if 2 * qi > N
            error("q = %d não eh compatível com N = %d.", qi, N);
        endif

        % Moving averages of q consecutive samples
        media = (S(:, qi+1:end) - S(:, 1:end-qi)) / qi;

        % Differences between averages separated by q
        diferencas = media(:, qi+1:end) - media(:, 1:end-qi);

        % Allan Variance
        AVAR_values_y(:, i) = 0.5 * mean(diferencas.^2, 2);

    endfor

endfunction
