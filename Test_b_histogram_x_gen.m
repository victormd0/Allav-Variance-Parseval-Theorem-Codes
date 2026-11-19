
function Test_b_histogram_x_gen(K, N, alpha, h, n, nbins, Q)

    % ============================================================
    % Default parameters
    % ============================================================

    if nargin < 6 || isempty(nbins)
        nbins = 30;
    endif

    if nargin < 7 || isempty(Q)
        Q = floor((N - 1)/2);
    endif

    % ============================================================
    % Validate inputs
    % ============================================================

    if ~ismember(alpha, [-2, -1, 0, 1, 2])
        error('alpha must be -2, -1, 0, 1, or 2.');
    endif

    if ~isscalar(K) || K < 1 || K ~= floor(K)
        error('K must be a positive integer.');
    endif

    if ~isscalar(N) || N < 4 || N ~= floor(N)
        error('N must be an integer greater than 3.');
    endif

    if ~isscalar(h) || h <= 0
        error('h must be positive.');
    endif

    if ~isscalar(n) || n < 0 || n ~= floor(n)
        error('n must be a nonnegative integer.');
    endif

    if ~isscalar(nbins) || nbins < 1 || nbins ~= floor(nbins)
        error('nbins must be a positive integer.');
    endif

    Qmax = floor((N - 1)/2);

    if ~isscalar(Q) || Q < 1 || Q > Qmax || Q ~= floor(Q)
        error('Q must be an integer between 1 and %d.', Qmax);
    endif

    if n > Q
        error('n must not exceed Q.');
    endif

    % ============================================================
    % Generate K independent realizations of x
    % ============================================================

    x = power_law(K, N, alpha, h);

    % ============================================================
    % Method 1: Your AVAR inversion formula
    % ============================================================

    b = b_values_x(x, n, Q);

    b_AVAR = b(:, 1);

    % ============================================================
    % Method 2: Direct autocorrelation of x
    % ============================================================

    x_centered = x - mean(x, 2);

    R = mean( ...
        x_centered(:, 1:N-n) .* ...
        x_centered(:, 1+n:N), 2);

    if n == 0
        b_auto = R;
    else
        b_auto = 2 * R;
    endif

    % ============================================================
    % Method 3: Theoretical b_n
    % ============================================================

    if alpha >= 0

        % Infinite-discrete theoretical coefficients

        if n == 0

            b_theory = h * 2^(-alpha) / (alpha + 1);

        elseif alpha == 0

            b_theory = 0;

        elseif alpha == 1

            b_theory = ...
                h * ((-1)^n - 1) / (pi^2*n^2);

        elseif alpha == 2

            b_theory = ...
                h * (-1)^n / (pi^2*n^2);

        endif

    else

        % Finite-FFT theoretical coefficients
        j = 1:N-1;
        f = min(j, N-j)/N;

        S = h * f.^alpha;

        R_theory = sum(S .* cos(2*pi*j*n/N))/N;

        if n == 0
            b_theory = R_theory;
        else
            b_theory = 2 * R_theory;
        endif

    endif

    % ============================================================
    % Statistics
    % ============================================================

    mean_AVAR = mean(b_AVAR);
    mean_auto = mean(b_auto);

    std_AVAR = std(b_AVAR);
    std_auto = std(b_auto);

    % Pearson correlation coefficient
    if K >= 2 && std_AVAR > 0 && std_auto > 0

        C = corrcoef(b_AVAR, b_auto);
        rho = C(1, 2);

    else

        rho = NaN;

    endif

    % ============================================================
    % Print statistical comparison table
    % ============================================================

    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('b_%d: K=%d, N=%d, alpha=%g, h=%g, Q=%d\n', ...
            n, K, N, alpha, h, Q);
    fprintf('========================================================\n\n');

    fprintf('%-20s %18s %18s\n', ...
            'Method', 'Mean', 'Std. Deviation');

    fprintf('%s\n', repmat('-', 1, 60));

    fprintf('%-20s %18.10e %18.10e\n', ...
            'AVAR inversion', mean_AVAR, std_AVAR);

    fprintf('%-20s %18.10e %18.10e\n', ...
            'Autocorrelation', mean_auto, std_auto);

    fprintf('%-20s %18.10e %18s\n', ...
            'Theoretical', b_theory, 'N/A');

    fprintf('%s\n', repmat('-', 1, 60));

    fprintf('\nPearson correlation (AVAR vs Auto): %.10f\n', rho);

    fprintf('Bias of AVAR estimator:            %.10e\n', ...
            mean_AVAR - b_theory);

    fprintf('Bias of autocorrelation estimator: %.10e\n\n', ...
            mean_auto - b_theory);

    % ============================================================
    % Common histogram bins
    % ============================================================

    all_values = [b_AVAR; b_auto; b_theory];

    xmin = min(all_values);
    xmax = max(all_values);

    if xmin == xmax
        xmin = xmin - 1;
        xmax = xmax + 1;
    endif

    margin = 0.05 * (xmax - xmin);

    xmin = xmin - margin;
    xmax = xmax + margin;

    edges = linspace(xmin, xmax, nbins+1);
    centers = (edges(1:end-1) + edges(2:end))/2;

    counts_AVAR = histc(b_AVAR, edges);
    counts_auto = histc(b_auto, edges);

    % Include observations on the final edge
    counts_AVAR(end-1) = ...
        counts_AVAR(end-1) + counts_AVAR(end);

    counts_auto(end-1) = ...
        counts_auto(end-1) + counts_auto(end);

    counts_AVAR = counts_AVAR(1:end-1);
    counts_auto = counts_auto(1:end-1);

    ymax = max([counts_AVAR(:); counts_auto(:); 1]) * 1.15;

    % ============================================================
    % Figure: Three histograms + realization scatter
    % ============================================================

    figure;

    % ------------------------------------------------------------
    % Panel 1: AVAR inversion histogram
    % ------------------------------------------------------------

    subplot(2, 2, 1);

    bar(centers, counts_AVAR, 0.85, ...
        'FaceColor', [0.2 0.4 0.9]);

    hold on;

    plot([b_theory b_theory], [0 ymax], 'r--', ...
         'LineWidth', 2);

    xlabel('b_n');
    ylabel('Frequency');
    title('AVAR inversion estimates');
    xlim([xmin xmax]);
    ylim([0 ymax]);
    grid on;

    % ------------------------------------------------------------
    % Panel 2: Autocorrelation histogram
    % ------------------------------------------------------------

    subplot(2, 2, 2);

    bar(centers, counts_auto, 0.85, ...
        'FaceColor', [0.2 0.7 0.3]);

    hold on;

    plot([b_theory b_theory], [0 ymax], 'r--', ...
         'LineWidth', 2);

    xlabel('b_n');
    ylabel('Frequency');
    title('Autocorrelation estimates');
    xlim([xmin xmax]);
    ylim([0 ymax]);
    grid on;

    % ------------------------------------------------------------
    % Panel 3: Both histograms together
    % ------------------------------------------------------------

    subplot(2, 2, 3);

    bar(centers, counts_AVAR, 0.85, ...
        'FaceColor', [0.2 0.4 0.9]);

    hold on;

    bar(centers, counts_auto, 0.45, ...
        'FaceColor', [0.2 0.7 0.3]);

    plot([b_theory b_theory], [0 ymax], 'r--', ...
         'LineWidth', 2);

    xlabel('b_n');
    ylabel('Frequency');
    title('Comparison of distributions');

    xlim([xmin xmax]);
    ylim([0 ymax]);

    legend('AVAR inversion', 'Autocorrelation', ...
           'Theoretical value');

    grid on;

    % ------------------------------------------------------------
    % Panel 4: Estimates by realization
    % ------------------------------------------------------------

    subplot(2, 2, 4);

    k = (1:K).';

    plot(k, b_AVAR, 'bo', ...
         'MarkerSize', 5);

    hold on;

    plot(k, b_auto, 'g^', ...
         'MarkerSize', 5);

    plot([1 max(2,K)], [b_theory b_theory], ...
         'r--', 'LineWidth', 2);

    xlabel('Realization k');
    ylabel('b_n');
    title('Estimates by realization');

    legend('AVAR inversion', 'Autocorrelation', ...
           'Theoretical value');

    xlim([1 max(2,K)]);
    grid on;

    % ============================================================
    % Overall title
    % ============================================================

    if exist('sgtitle', 'file')
        sgtitle(sprintf( ...
            'b_%d: K=%d, N=%d, alpha=%g, h=%g, Q=%d', ...
            n, K, N, alpha, h, Q));
    endif

endfunction

