
function Test_b_values_x_gen(K, N, alpha, h, n, plot_type)

    % ============================================================
    % Default parameters and validation
    % ============================================================

    if nargin < 5
        n = 0:100;
    endif

    if nargin < 6
        plot_type = 'plot';
    endif

    if ~ismember(alpha, [-2, -1, 0, 1, 2])
        error('alpha must be -2, -1, 0, 1, or 2.');
    endif

    if ~strcmp(plot_type, 'plot') && ...
       ~strcmp(plot_type, 'loglog')
        error('plot_type must be "plot" or "loglog".');
    endif

    if any(n < 0) || any(n ~= floor(n)) || any(n >= N)
        error('n must contain integers between 0 and N-1.');
    endif

    n = n(:).';

    % ============================================================
    % Generate x with PSD S_x(f) = h*|f|^alpha
    % ============================================================

    x = power_law(K, N, alpha, h);

    % ============================================================
    % Method 1: Your AVAR inversion formula
    % ============================================================

    % b_values_x must calculate b_n from AVAR of diff(x)
    b = b_values_x(x, n);
    b_AVAR = mean(b, 1);

    % ============================================================
    % Method 2: Autocorrelation of x
    % ============================================================

    x_centered = x - mean(x, 2);

    b_autocorrelation = zeros(size(n));

    for i = 1:length(n)

        ni = n(i);

        R = mean(mean( ...
            x_centered(:, 1:N-ni) .* ...
            x_centered(:, 1+ni:N), 2));

        if ni == 0
            b_autocorrelation(i) = R;
        else
            b_autocorrelation(i) = 2 * R;
        endif

    endfor

    % ============================================================
    % Method 3: Theoretical b_n
    % ============================================================

    b_theoretical = zeros(size(n));

    if alpha >= 0

        for i = 1:length(n)

            ni = n(i);

            if ni == 0

                b_theoretical(i) = ...
                    h * 2^(-alpha) / (alpha + 1);

            elseif alpha == 0

                b_theoretical(i) = 0;

            elseif alpha == 1

                b_theoretical(i) = ...
                    h * ((-1)^ni - 1) / (pi^2 * ni^2);

            elseif alpha == 2

                b_theoretical(i) = ...
                    h * (-1)^ni / (pi^2 * ni^2);

            endif

        endfor

    else

        % Finite FFT theoretical coefficients
        k = 1:N-1;
        f = min(k, N-k)/N;
        S = h * f.^alpha;

        for i = 1:length(n)

            ni = n(i);

            R = sum(S .* cos(2*pi*k*ni/N))/N;

            if ni == 0
                b_theoretical(i) = R;
            else
                b_theoretical(i) = 2 * R;
            endif

        endfor

    endif

    % ============================================================
    % Absolute differences
    % ============================================================

    % AVAR vs theoretical
    absolute_theory = abs(b_AVAR - b_theoretical);

    % AVAR vs autocorrelation
    absolute_auto = abs(b_AVAR - b_autocorrelation);

    % ============================================================
    % Relative differences (%)
    % ============================================================

    relative_theory = NaN(size(n));
    relative_auto = NaN(size(n));

    % Denominator: theoretical b_n
    valid_theory = (b_theoretical ~= 0);

    relative_theory(valid_theory) = ...
        100 * absolute_theory(valid_theory) ...
        ./ abs(b_theoretical(valid_theory));

    % Denominator: autocorrelation b_n
    valid_auto = (b_autocorrelation ~= 0);

    relative_auto(valid_auto) = ...
        100 * absolute_auto(valid_auto) ...
        ./ abs(b_autocorrelation(valid_auto));

    % ============================================================
    % Figure 1: All three b_n values
    % ============================================================

    figure;

    if strcmp(plot_type, 'loglog')

        v1 = (n > 0) & (abs(b_AVAR) > 0);
        v2 = (n > 0) & (abs(b_autocorrelation) > 0);
        v3 = (n > 0) & (abs(b_theoretical) > 0);

        loglog(n(v1), abs(b_AVAR(v1)), ...
               'bo-', 'LineWidth', 1.5);
        hold on;

        loglog(n(v2), abs(b_autocorrelation(v2)), ...
               'g^-', 'LineWidth', 1.5);

        loglog(n(v3), abs(b_theoretical(v3)), ...
               'rs--', 'LineWidth', 2);

        ylabel('|b_n|');

    else

        plot(n, b_AVAR, 'bo-', 'LineWidth', 1.5);
        hold on;

        plot(n, b_autocorrelation, 'g^-', ...
             'LineWidth', 1.5);

        plot(n, b_theoretical, 'rs--', ...
             'LineWidth', 2);

        ylabel('b_n');

    endif

    xlabel('n');

    title(sprintf('b_n Comparison: alpha=%g, h=%g', ...
                  alpha, h));

    legend('From AVAR', 'From Autocorrelation', ...
           'Theoretical');

    grid on;

    % ============================================================
    % Figure 2: Absolute difference (AVAR vs theory)
    % ============================================================

    figure;

    plot(n, absolute_theory, 'bo-', ...
         'LineWidth', 1.5);

    xlabel('n');
    ylabel('Absolute Difference');

    title(sprintf('AVAR vs Theory: Absolute Difference, alpha=%g', ...
                  alpha));

    grid on;

    % ============================================================
    % Figure 3: Relative difference (AVAR vs theory)
    % ============================================================

    figure;

    plot(n, relative_theory, 'bo-', ...
         'LineWidth', 1.5);

    xlabel('n');
    ylabel('Relative Difference (%)');

    title(sprintf('AVAR vs Theory: Relative Difference, alpha=%g', ...
                  alpha));

    grid on;

    % ============================================================
    % Figure 4: Absolute difference (AVAR vs autocorrelation)
    % ============================================================

    figure;

    plot(n, absolute_auto, 'go-', ...
         'LineWidth', 1.5);

    xlabel('n');
    ylabel('Absolute Difference');

    title(sprintf('AVAR vs Autocorrelation: Absolute Difference, alpha=%g', ...
                  alpha));

    grid on;

    % ============================================================
    % Figure 5: Relative difference (AVAR vs autocorrelation)
    % ============================================================

    figure;

    plot(n, relative_auto, 'go-', ...
         'LineWidth', 1.5);

    xlabel('n');
    ylabel('Relative Difference (%)');

    title(sprintf('AVAR vs Autocorrelation: Relative Difference, alpha=%g', ...
                  alpha));

    grid on;

    % ============================================================
    % Print comparison table
    % ============================================================

    fprintf('\nb_n comparison: K=%d, N=%d, alpha=%g, h=%g\n', ...
            K, N, alpha, h);

    fprintf('%8s %18s %18s %18s %16s %16s\n', ...
            'n', 'AVAR b_n', 'Autocorr b_n', ...
            'Theoretical b_n', ...
            'Rel. Theory (%)', 'Rel. Auto (%)');

    for i = 1:length(n)

        fprintf('%8d %18.8e %18.8e %18.8e %16.5f %16.5f\n', ...
                n(i), b_AVAR(i), ...
                b_autocorrelation(i), b_theoretical(i), ...
                relative_theory(i), relative_auto(i));

    endfor

endfunction

