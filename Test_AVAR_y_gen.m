
function Test_AVAR_y_gen(K, N, alpha, h, q)
% Test Allan variance of power-law Gaussian noise.
%
% INPUT:
%   K     : number of independent realizations
%   N     : samples per realization
%   alpha : spectral exponent (-2, -1, 0, 1, 2)
%   h     : two-sided PSD coefficient
%   q     : vector of averaging factors (optional)
%
% Requires:
%   power_law.m
%   AVAR_values_y.m

    if nargin < 5
        q = 2.^(0:floor(log2(N/2)));
    endif

    q = q(:)';

    if any(q < 1 | q ~= floor(q) | 2*q > N)
        error('q must contain positive integers with 2*q <= N.');
    endif

    % Generate power-law noise
    y = power_law(K, N, alpha, h);

    % Numerical Allan variance
    AVAR = AVAR_values_y(y, q);
    AVAR_numerical = mean(AVAR, 1);

    % Theoretical Allan variance from discrete PSD
    k = 1:N-1;
    f = min(k, N-k) / N;

    S = h * f.^alpha;

    AVAR_theoretical = zeros(size(q));

    denominator = sin(pi*f).^2;

    for i = 1:length(q)

        qi = q(i);

        H = 2 * sin(pi*qi*f).^4 ...
              ./ (qi^2 * denominator);

        AVAR_theoretical(i) = sum(S .* H) / N;

    endfor

    % Relative error
    relative_error = abs(AVAR_numerical - AVAR_theoretical) ...
                     ./ AVAR_theoretical;

    relative_error_percent = 100 * relative_error;

    % Plot numerical and theoretical AVAR
    figure;

    loglog(q, AVAR_numerical, 'o-', 'LineWidth', 2);
    hold on;

    loglog(q, AVAR_theoretical, '--', 'LineWidth', 2);

    xlabel('Averaging factor q');
    ylabel('Allan Variance');
    legend('Numerical', 'Theoretical');
    title(sprintf('AVAR of Power-Law Noise: alpha=%g, h=%g', ...
                  alpha, h));
    grid on;

    % Plot relative error
    figure;

    semilogx(q, relative_error_percent, 'o-', ...
             'LineWidth', 2);

    xlabel('Averaging factor q');
    ylabel('Relative error (%)');
    title(sprintf('Relative AVAR Error: alpha=%g', alpha));
    grid on;

    % Display results
    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('              Allan Variance y Test Results\n');
    fprintf('========================================================\n');
    fprintf('K = %d, N = %d, alpha = %g, h = %g\n', ...
            K, N, alpha, h);
    fprintf('--------------------------------------------------------\n');

    fprintf('%8s %18s %18s %14s\n', ...
            'q', 'Numerical AVAR', 'Theoretical AVAR', ...
            'Error (%)');

    fprintf('--------------------------------------------------------\n');

    for i = 1:length(q)
        fprintf('%8d %18.8e %18.8e %13.4f%%\n', ...
                q(i), AVAR_numerical(i), ...
                AVAR_theoretical(i), ...
                relative_error_percent(i));
    endfor

    fprintf('========================================================\n');

endfunction

