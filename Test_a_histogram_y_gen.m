
function Test_a_histogram_y_gen(K, N, alpha, h, n, nbins)

    % ============================================================
    % Default parameters
    % ============================================================

    if nargin < 6
        nbins = 10;
    endif

    % ============================================================
    % Validate inputs
    % ============================================================

    if ~ismember(alpha, [-2, -1, 0, 1, 2])
        error('alpha must be -2, -1, 0, 1, or 2.');
    endif

    if K < 1 || K ~= floor(K)
        error('K must be a positive integer.');
    endif

    if N < 4 || N ~= floor(N)
        error('N must be an integer greater than 3.');
    endif

    if h <= 0
        error('h must be positive.');
    endif

    if ~isscalar(n) || n < 0 || n ~= floor(n) || n >= N
        error('n must be an integer between 0 and N-1.');
    endif

    if ~isscalar(nbins) || nbins < 1 || nbins ~= floor(nbins)
        error('nbins must be a positive integer.');
    endif

    % ============================================================
    % Generate K independent realizations
    % ============================================================

    y = power_law(K, N, alpha, h);

    % ============================================================
    % Method 1: Your AVAR inversion formula
    % ============================================================

    a = a_values_y(y, n);
    a_AVAR = a(:, 1);

    % ============================================================
    % Method 2: Direct autocorrelation
    % ============================================================

    y_centered = y - mean(y, 2);

    R = mean( ...
        y_centered(:, 1:N-n) .* ...
        y_centered(:, 1+n:N), 2);

    if n == 0
        a_auto = R;
    else
        a_auto = 2 * R;
    endif

    % ============================================================
    % Method 3: Theoretical coefficient
    % ============================================================

    if alpha >= 0

        if n == 0

            a_theory = h * 2^(-alpha) / (alpha + 1);

        elseif alpha == 0

            a_theory = 0;

        elseif alpha == 1

            a_theory = h * ((-1)^n - 1) / (pi^2*n^2);

        elseif alpha == 2

            a_theory = h * (-1)^n / (pi^2*n^2);

        endif

    else

        % Finite-FFT theoretical coefficient
        j = 1:N-1;
        f = min(j, N-j)/N;

        S = h * f.^alpha;

        R_theory = sum(S .* cos(2*pi*j*n/N))/N;

        if n == 0
            a_theory = R_theory;
        else
            a_theory = 2 * R_theory;
        endif

    endif

    % ============================================================
    % Statistics
    % ============================================================

    mean_AVAR = mean(a_AVAR);
    mean_auto = mean(a_auto);

    std_AVAR = std(a_AVAR);
    std_auto = std(a_auto);

    % Pearson correlation coefficient
    if K >= 2 && std_AVAR > 0 && std_auto > 0

        C = corrcoef(a_AVAR, a_auto);
        rho = C(1, 2);

    else

        rho = NaN;

    endif

    % ============================================================
    % Print statistical comparison table
    % ============================================================

    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('a_%d: K=%d, N=%d, alpha=%g, h=%g\n', ...
            n, K, N, alpha, h);
    fprintf('========================================================\n\n');

    fprintf('%-20s %18s %18s\n', ...
            'Method', 'Mean', 'Std. Deviation');

    fprintf('%s\n', repmat('-', 1, 60));

    fprintf('%-20s %18.10e %18.10e\n', ...
            'AVAR inversion', mean_AVAR, std_AVAR);

    fprintf('%-20s %18.10e %18.10e\n', ...
            'Autocorrelation', mean_auto, std_auto);

    fprintf('%-20s %18.10e %18s\n', ...
            'Theoretical', a_theory, 'N/A');

    fprintf('%s\n', repmat('-', 1, 60));

    fprintf('\nPearson correlation (AVAR vs Auto): %.10f\n', rho);

    fprintf('Bias of AVAR estimator:           %.10e\n', ...
            mean_AVAR - a_theory);

    fprintf('Bias of autocorrelation estimator: %.10e\n\n', ...
            mean_auto - a_theory);

    % ============================================================
    % Common histogram bins
    % ============================================================

    all_values = [a_AVAR; a_auto; a_theory];

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

    counts_AVAR = histc(a_AVAR, edges);
    counts_auto = histc(a_auto, edges);

    % Include values falling on the final edge
    counts_AVAR(end-1) = counts_AVAR(end-1) + counts_AVAR(end);
    counts_auto(end-1) = counts_auto(end-1) + counts_auto(end);

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

    plot([a_theory a_theory], [0 ymax], 'r--', ...
         'LineWidth', 2);

    xlabel('a_n');
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

    plot([a_theory a_theory], [0 ymax], 'r--', ...
         'LineWidth', 2);

    xlabel('a_n');
    ylabel('Frequency');
    title('Autocorrelation estimates');
    xlim([xmin xmax]);
    ylim([0 ymax]);
    grid on;

    % ------------------------------------------------------------
    % Panel 3: Both distributions
    % ------------------------------------------------------------

    subplot(2, 2, 3);

    bar(centers, counts_AVAR, 0.85, ...
        'FaceColor', [0.2 0.4 0.9]);

    hold on;

    bar(centers, counts_auto, 0.45, ...
        'FaceColor', [0.2 0.7 0.3]);

    plot([a_theory a_theory], [0 ymax], 'r--', ...
         'LineWidth', 2);

    xlabel('a_n');
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

    plot(k, a_AVAR, 'bo', ...
         'MarkerSize', 5);

    hold on;

    plot(k, a_auto, 'g^', ...
         'MarkerSize', 5);

    plot([1 max(2,K)], [a_theory a_theory], ...
         'r--', 'LineWidth', 2);

    xlabel('Realization k');
    ylabel('a_n');
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
            'a_%d: K=%d, N=%d, alpha=%g, h=%g', ...
            n, K, N, alpha, h));
    endif

endfunction

