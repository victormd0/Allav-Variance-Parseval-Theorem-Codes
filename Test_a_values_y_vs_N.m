
function Test_a_values_y_vs_N(K, N, alpha, h, n, plot_type)

    if nargin < 6
        plot_type = 'plot';
    endif

    if ~ismember(plot_type, {'plot', 'semilog', 'loglog'})
        error("plot_type must be 'plot', 'semilog', or 'loglog'.");
    endif

    % ============================================================
    % Validate inputs
    % ============================================================

    if ~isscalar(K) || K < 1 || K ~= floor(K)
        error('K must be a positive integer.');
    endif

    if isempty(N) || any(N < 4) || any(N ~= floor(N))
        error('N must be a vector of integers >= 4.');
    endif

    if ~isscalar(alpha) || ~ismember(alpha, [-2,-1,0,1,2])
        error('alpha must be -2, -1, 0, 1, or 2.');
    endif

    if ~isscalar(h) || h <= 0
        error('h must be positive.');
    endif

    if ~isscalar(n) || n < 0 || n ~= floor(n)
        error('n must be a nonnegative integer.');
    endif

    % Sort N
    N = sort(N(:).');
    M = length(N);

    if any(2*N*0 + 2*n > N)
        error('Each N must satisfy 2*n <= N.');
    endif

    % ============================================================
    % Allocate results
    % ============================================================

    mean_AVAR = zeros(1,M);
    mean_auto = zeros(1,M);
    a_theory = zeros(1,M);

    % ============================================================
    % Loop over signal lengths
    % ============================================================

    for i = 1:M

        Ni = N(i);

        fprintf('Processing N = %d (%d/%d)\n', Ni, i, M);

        % Generate K realizations of length Ni
        y = power_law(K, Ni, alpha, h);

        % --------------------------------------------------------
        % Method 1: AVAR inversion
        % --------------------------------------------------------

        a = a_values_y(y, n);

        mean_AVAR(i) = mean(a(:,1));

        % --------------------------------------------------------
        % Method 2: Direct autocorrelation
        % --------------------------------------------------------

        y_centered = y - mean(y, 2);

        R = mean( ...
            y_centered(:,1:Ni-n) .* ...
            y_centered(:,1+n:Ni), 2);

        if n == 0
            a_auto = R;
        else
            a_auto = 2*R;
        endif

        mean_auto(i) = mean(a_auto);

        % --------------------------------------------------------
        % Theoretical coefficient
        % --------------------------------------------------------

        if alpha >= 0

            if n == 0

                a_theory(i) = ...
                    h * 2^(-alpha)/(alpha+1);

            elseif alpha == 0

                a_theory(i) = 0;

            elseif alpha == 1

                a_theory(i) = ...
                    h*((-1)^n-1)/(pi^2*n^2);

            elseif alpha == 2

                a_theory(i) = ...
                    h*(-1)^n/(pi^2*n^2);

            endif

        else

            % Finite-FFT theoretical coefficient
            j = 1:Ni-1;
            f = min(j,Ni-j)/Ni;

            S = h*f.^alpha;

            R_theory = ...
                sum(S.*cos(2*pi*j*n/Ni))/Ni;

            if n == 0
                a_theory(i) = R_theory;
            else
                a_theory(i) = 2*R_theory;
            endif

        endif

    endfor

    % ============================================================
    % Absolute differences
    % ============================================================

    diff_AVAR_theory = abs(mean_AVAR - a_theory);

    diff_auto_theory = abs(mean_auto - a_theory);

    diff_AVAR_auto = abs(mean_AVAR - mean_auto);

    % ============================================================
    % Print results table
    % ============================================================

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('a_%d: K=%d, alpha=%g, h=%g\n', n, K, alpha, h);
    fprintf('============================================================\n');

    fprintf('%8s %16s %16s %16s\n', ...
        'N', 'AVAR mean', 'Auto mean', 'Theory');

    fprintf('%s\n', repmat('-',1,60));

    for i = 1:M
        fprintf('%8d %16.8e %16.8e %16.8e\n', ...
            N(i), mean_AVAR(i), mean_auto(i), a_theory(i));
    endfor

    % ============================================================
    % Plot all four graphs
    % ============================================================

    figure;

    % ------------------------------------------------------------
    % Plot 1: Mean coefficients
    % ------------------------------------------------------------

    subplot(2,2,1);

    plot(N, mean_AVAR, 'bo-', 'LineWidth',1.5);
    hold on;

    plot(N, mean_auto, 'gs-', 'LineWidth',1.5);

    plot(N, a_theory, 'r--', 'LineWidth',2);

    xlabel('N');
    ylabel('a_n');
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
         'LineWidth',1.5);

    xlabel('N');
    ylabel('Absolute difference');
    title('|AVAR - Theory|');
    grid on;

    % ------------------------------------------------------------
    % Plot 3: Autocorrelation vs theory
    % ------------------------------------------------------------

    subplot(2,2,3);

    plot(N, diff_auto_theory, 'gs-', ...
         'LineWidth',1.5);

    xlabel('N');
    ylabel('Absolute difference');
    title('|Autocorrelation - Theory|');
    grid on;

    % ------------------------------------------------------------
    % Plot 4: AVAR vs autocorrelation
    % ------------------------------------------------------------

    subplot(2,2,4);

    plot(N, diff_AVAR_auto, 'mo-', ...
         'LineWidth',1.5);

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
                % Main comparison graph
                set(gca, 'YScale', 'linear');
            else
                % Three absolute-difference graphs
                set(gca, 'YScale', 'log');
            endif

        endif

    endfor

    % ============================================================
    % Overall title
    % ============================================================

    if exist('sgtitle','file')
        sgtitle(sprintf( ...
            'a_%d: K=%d, alpha=%g, h=%g', ...
            n, K, alpha, h));
    endif

endfunction

