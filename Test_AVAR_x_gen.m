
function Test_AVAR_x_gen(K, N, alpha, h, q)

    if nargin < 5
        q = 2.^(0:floor(log2((N-1)/2)));
    endif

    if ~ismember(alpha, [-2, -1, 0, 1, 2])
        error('alpha must be -2, -1, 0, 1, or 2.');
    endif

    if any(q < 1) || any(q ~= floor(q)) || any(2*q > N-1)
        error('q must contain positive integers with 2*q <= N-1.');
    endif

    % Generate K realizations of x
    x = power_law(K, N, alpha, h);

    % Numerical AVAR of the differences of x
    AVAR = AVAR_values_x(x, q);
    AVAR_numerical = mean(AVAR, 1);

    % Theoretical AVAR for the finite FFT model
    k = 1:N-1;
    f = min(k, N-k)/N;

    Sx = h * f.^alpha;

    AVAR_theoretical = zeros(1, length(q));

    for i = 1:length(q)
        qi = q(i);

        H = 8 * sin(pi*qi*f).^4 / qi^2;

        AVAR_theoretical(i) = sum(Sx .* H)/N;
    endfor

    % Relative error
    relative_error = 100 * abs(AVAR_numerical - AVAR_theoretical) ...
                         ./ AVAR_theoretical;

    % Plot AVAR
    figure;
    loglog(q, AVAR_numerical, 'bo-', 'LineWidth', 1.5);
    hold on;
    loglog(q, AVAR_theoretical, 'r--', 'LineWidth', 2);
    grid on;

    xlabel('q');
    ylabel('AVAR');
    title(sprintf('AVAR x: K=%d, N=%d, alpha=%g, h=%g', ...
                  K, N, alpha, h));
    legend('Numerical AVAR', 'Theoretical AVAR');

    % Plot relative error
    figure;
    semilogx(q, relative_error, 'ko-', 'LineWidth', 1.5);
    grid on;

    xlabel('q');
    ylabel('Relative error (%)');
    title(sprintf('AVAR x relative error: alpha=%g, h=%g', ...
                  alpha, h));

    % Print results
    fprintf('\nAVAR x: K=%d, N=%d, alpha=%g, h=%g\n', ...
            K, N, alpha, h);

    fprintf('%8s %20s %20s %15s\n', ...
            'q', 'Numerical AVAR', 'Theoretical AVAR', 'Error (%)');

    for i = 1:length(q)
        fprintf('%8d %20.10e %20.10e %15.6f\n', ...
                q(i), AVAR_numerical(i), ...
                AVAR_theoretical(i), relative_error(i));
    endfor

endfunction

