
function Test_b_values_x_vs_N(K, N, alpha, h, n, plot_type)

    % ============================================================
    % Default parameters
    % ============================================================

    if nargin < 6 || isempty(plot_type)
        plot_type = 'plot';
    endif

    % ============================================================
    % Validate inputs
    % ============================================================

    if ~isscalar(K) || K < 1 || K ~= floor(K)
        error('K must be a positive integer.');
    endif

    if isempty(N) || ~isvector(N) || ...
       any(N < 4) || any(N ~= floor(N))
        error('N must be a vector of integers >= 4.');
    endif

    if ~isscalar(alpha) || ...
       ~ismember(alpha, [-2, -1, 0, 1, 2])
        error('alpha must be -2, -1, 0, 1, or 2.');
    endif

    if ~isscalar(h) || h <= 0
        error('h must be positive.');
    endif

    if ~isscalar(n) || n < 0 || n ~= floor(n)
        error('n must be a nonnegative integer.');
    endif

    if ~any(strcmp(plot_type, ...
           {'plot', 'semilog', 'loglog'}))
        error("plot_type must be 'plot', 'semilog', or 'loglog'.");
    endif

    % Sort signal lengths
    N = sort(N(:).');
    M = length(N);

    % b_values_x requires n <= Q
    % where Q = floor((N-1)/2)

    if any(n > floor((N-1)/2))
        error('Each N must satisfy n <= floor((N-1)/2).');
    endif

    % ============================================================
    % Allocate results
    % ============================================================

    mean_AVAR = zeros(1, M);
    mean_auto = zeros(1, M);
    b_theory = zeros(1, M);

    std_AVAR = zeros(1, M);
    std_auto = zeros(1, M);

    % ============================================================
    % Loop over signal lengths
    % ============================================================

    for i = 1:M

        Ni = N(i);

        fprintf('Processing N = %d (%d/%d)\n', Ni, i, M);

        % --------------------------------------------------------
        % Generate K independent realizations of x
        % --------------------------------------------------------

        x = power_law(K, Ni, alpha, h);

        % --------------------------------------------------------
        % Method 1: AVAR inversion
        % --------------------------------------------------------

        b = b_values_x(x, n);

        b_AVAR = b(:, 1);

        mean_AVAR(i) = mean(b_AVAR);
        std_AVAR(i) = std(b_AVAR);

        % --------------------------------------------------------
        % Method 2: Direct autocorrelation
        % --------------------------------------------------------

        x_centered = x - mean(x, 2);

        R = mean( ...
            x_centered(:, 1:Ni-n) .* ...
            x_centered(:, 1+n:Ni), 2);

        if n == 0
            b_auto = R;
        else
            b_auto = 2 * R;
        endif

        mean_auto(i) = mean(b_auto);
        std_auto(i) = std(b_auto);

        % --------------------------------------------------------
        % Theoretical coefficient
        % --------------------------------------------------------

        if alpha >= 0

            % Infinite-discrete theoretical coefficient

            if n == 0

                b_theory(i) = ...
                    h * 2^(-alpha) / (alpha + 1);

            elseif alpha == 0

                b_theory(i) = 0;

            elseif alpha == 1

                b_theory(i) = ...
                    h * ((-1)^n - 1) / (pi^2*n^2);

            elseif alpha == 2

                b_theory(i) = ...
                    h * (-1)^n / (pi^2*n^2);

            endif

        else

            % Finite-FFT theoretical coefficient

            j = 1:Ni-1;
            f = min(j, Ni-j)/Ni;

            S = h * f.^alpha;

            R_theory = ...
                sum(S .* cos(2*pi*j*n/Ni))/Ni;

            if n == 0
                b_theory(i) = R_theory;
            else
                b_theory(i) = 2*R_theory;
            endif

        endif

    endfor

    % ============================================================
    % Absolute differences between the mean estimates
    % ============================================================

    diff_AVAR_theory = abs(mean_AVAR - b_theory);

    diff_auto_theory = abs(mean_auto - b_theory);

    diff_AVAR_auto = abs(mean_AVAR - mean_auto);

    % ============================================================
    % Print results table
    % ============================================================

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('b_%d: K=%d, alpha=%g, h=%g\n', ...
            n, K, alpha, h);
    fprintf('============================================================\n\n');

    fprintf('%10s %16s %16s %16s\n', ...
            'N', 'AVAR mean', 'Auto mean', 'Theory');

    fprintf('%s\n', repmat('-', 1, 64));

    for i = 1:M

        fprintf('%10d %16.8e %16.8e %16.8e\n', ...
                N(i), mean_AVAR(i), ...
                mean_auto(i), b_theory(i));

    endfor

    fprintf('\n');

    % ============================================================
    % Plot all four graphs
    % ============================================================

    figure;

    % ------------------------------------------------------------
    % Plot 1: Mean coefficients
    % ------------------------------------------------------------

    subplot(2,2,1);

    plot(N, mean_AVAR, 'bo-', 'LineWidth', 1.5);
    hold on;

    plot(N, mean_auto, 'gs-', 'LineWidth', 1.5);

    plot(N, b_theory, 'r--', 'LineWidth', 2);

    xlabel('N');
    ylabel('b_n');
    title('Mean estimates vs N');

    legend('AVAR inversion', ...
           'Autocorrelation', ...
           'Theoretical value');

    grid on;

    % ------------------------------------------------------------
    % Plot 2: AVAR vs theory
    % ------------------------------------------------------------

    subplot(2,2,2);

    plot(N, diff_AVAR_theory, 'bo-', ...
         'LineWidth', 1.5);

    xlabel('N');
    ylabel('Absolute difference');
    title('|AVAR - Theory|');

    grid on;

    % ------------------------------------------------------------
    % Plot 3: Autocorrelation vs theory
    % ------------------------------------------------------------

    subplot(2,2,3);

    plot(N, diff_auto_theory, 'gs-', ...
         'LineWidth', 1.5);

    xlabel('N');
    ylabel('Absolute difference');
    title('|Autocorrelation - Theory|');

    grid on;

    % ------------------------------------------------------------
    % Plot 4: AVAR vs autocorrelation
    % ------------------------------------------------------------

    subplot(2,2,4);

    plot(N, diff_AVAR_auto, 'mo-', ...
         'LineWidth', 1.5);

    xlabel('N');
    ylabel('Absolute difference');
    title('|AVAR - Autocorrelation|');

    grid on;

    % ============================================================
    % Set axis scales
    % ============================================================

    for i = 1:4

        subplot(2,2,i);

        if strcmp(plot_type, 'plot')

            set(gca, 'XScale', 'linear');
            set(gca, 'YScale', 'linear');

        elseif strcmp(plot_type, 'semilog')

            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'linear');

        elseif strcmp(plot_type, 'loglog')

            set(gca, 'XScale', 'log');

            if i == 1

                % The main graph may contain negative b_n
                set(gca, 'YScale', 'linear');

            else

                % Absolute-difference graphs
                set(gca, 'YScale', 'log');

            endif

        endif

    endfor

    % ============================================================
    % Overall title
    % ============================================================

    if exist('sgtitle', 'file')

        sgtitle(sprintf( ...
            'b_%d: K=%d, alpha=%g, h=%g', ...
            n, K, alpha, h));

    endif

endfunction

