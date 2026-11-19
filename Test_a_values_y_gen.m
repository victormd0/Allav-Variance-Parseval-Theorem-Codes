
function Test_a_values_y_gen(K, N, alpha, h, n, plot_type)

    % ============================================================
    % Default parameters
    % ============================================================

    if nargin < 5
        n = 0:100;
    endif

    if nargin < 6
        plot_type = 'plot';
    endif

    % Validate inputs
    if ~ismember(alpha, [-2, -1, 0, 1, 2])
        error('alpha must be -2, -1, 0, 1, or 2.');
    endif

    if ~strcmp(plot_type, 'plot') && ~strcmp(plot_type, 'loglog')
        error('plot_type must be "plot" or "loglog".');
    endif

    if any(n < 0) || any(n ~= floor(n)) || any(n >= N)
        error('n must contain integers between 0 and N-1.');
    endif

    n = n(:).';

    % ============================================================
    % Generate K realizations of power-law noise
    % ============================================================

    y = power_law(K, N, alpha, h);

    % ============================================================
    % Method 1: Your AVAR inversion formula
    % ============================================================

    a = a_values_y(y, n);
    a_numerical = mean(a, 1);

    % ============================================================
    % Method 2: Autocorrelation
    % ============================================================

    y_centered = y - mean(y, 2);

    a_autocorrelation = zeros(size(n));

    for i = 1:length(n)

        ni = n(i);

        R = mean(mean( ...
            y_centered(:, 1:N-ni) .* ...
            y_centered(:, 1+ni:N), 2));

        if ni == 0
            a_autocorrelation(i) = R;
        else
            a_autocorrelation(i) = 2 * R;
        endif

    endfor

    % ============================================================
    % Method 3: Theoretical a_n
    % ============================================================

    a_theoretical = zeros(size(n));

    if alpha >= 0

        % Infinite-discrete analytical coefficients
        for i = 1:length(n)

            ni = n(i);

            if ni == 0

                a_theoretical(i) = ...
                    h * 2^(-alpha) / (alpha + 1);

            elseif alpha == 0

                a_theoretical(i) = 0;

            elseif alpha == 1

                a_theoretical(i) = ...
                    h * ((-1)^ni - 1) / (pi^2 * ni^2);

            elseif alpha == 2

                a_theoretical(i) = ...
                    h * (-1)^ni / (pi^2 * ni^2);

            endif

        endfor

    else

        % Finite-FFT theoretical coefficients
        k = 1:N-1;
        f = min(k, N-k)/N;

        S = h * f.^alpha;

        for i = 1:length(n)

            ni = n(i);

            R = sum(S .* cos(2*pi*k*ni/N))/N;

            if ni == 0
                a_theoretical(i) = R;
            else
                a_theoretical(i) = 2 * R;
            endif

        endfor

    endif

    % ============================================================
    % Absolute differences
    % ============================================================

    % Your formula vs theory
    absolute_theory = abs(a_numerical - a_theoretical);

    % Your formula vs autocorrelation
    absolute_auto = abs(a_numerical - a_autocorrelation);

    % ============================================================
    % Relative differences (%)
    % ============================================================

    relative_theory = NaN(size(n));
    relative_auto = NaN(size(n));

    % Denominator: theoretical a_n
    valid_theory = (a_theoretical ~= 0);

    relative_theory(valid_theory) = ...
        100 * absolute_theory(valid_theory) ...
        ./ abs(a_theoretical(valid_theory));

    % Denominator: autocorrelation a_n
    valid_auto = (a_autocorrelation ~= 0);

    relative_auto(valid_auto) = ...
        100 * absolute_auto(valid_auto) ...
        ./ abs(a_autocorrelation(valid_auto));

    % ============================================================
    % Figure 1: All three a_n values
    % ============================================================

    figure;

    if strcmp(plot_type, 'loglog')

        valid_num = (n > 0) & (abs(a_numerical) > 0);
        valid_auto = (n > 0) & (abs(a_autocorrelation) > 0);
        valid_theo = (n > 0) & (abs(a_theoretical) > 0);

        loglog(n(valid_num), abs(a_numerical(valid_num)), ...
               'bo-', 'LineWidth', 1.5);
        hold on;

        loglog(n(valid_auto), abs(a_autocorrelation(valid_auto)), ...
               'g^-', 'LineWidth', 1.5);

        loglog(n(valid_theo), abs(a_theoretical(valid_theo)), ...
               'rs--', 'LineWidth', 2);

        ylabel('|a_n|');

    else

        plot(n, a_numerical, 'bo-', 'LineWidth', 1.5);
        hold on;

        plot(n, a_autocorrelation, 'g^-', ...
             'LineWidth', 1.5);

        plot(n, a_theoretical, 'rs--', ...
             'LineWidth', 2);

        ylabel('a_n');

    endif

    xlabel('n');

    title(sprintf('a_n Comparison: K=%d, N=%d, alpha=%g, h=%g', ...
                  K, N, alpha, h));

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

    fprintf('\na_n comparison: K=%d, N=%d, alpha=%g, h=%g\n', ...
            K, N, alpha, h);

    fprintf('%8s %18s %18s %18s %16s %16s\n', ...
            'n', 'AVAR a_n', 'Autocorr a_n', ...
            'Theoretical a_n', ...
            'Rel. Theory (%)', 'Rel. Auto (%)');

    for i = 1:length(n)

        fprintf('%8d %18.8e %18.8e %18.8e %16.5f %16.5f\n', ...
                n(i), a_numerical(i), ...
                a_autocorrelation(i), a_theoretical(i), ...
                relative_theory(i), relative_auto(i));

    endfor

endfunction

