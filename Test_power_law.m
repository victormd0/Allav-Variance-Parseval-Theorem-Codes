
function Test_power_law(K, N, alpha, h)

    % Generate power-law noise
    y = power_law(K, N, alpha, h);

    % Compute FFT and average the periodograms
    Y = fft(y, [], 2);
    PSD = mean(abs(Y).^2, 1) / N;

    % Positive frequencies (excluding DC and Nyquist)
    m = 1:floor((N-1)/2);
    f = m / N;

    % Numerical and theoretical PSD
    PSD_numerical = PSD(m+1);
    PSD_theoretical = h * f.^alpha;

    % Plot
    figure;

    loglog(f, PSD_numerical, 'b', 'LineWidth', 1);
    hold on;

    loglog(f, PSD_theoretical, 'r--', 'LineWidth', 2);

    xlabel('Frequency (cycles/sample)');
    ylabel('Two-sided PSD');
    title(sprintf('Power-law noise: K=%d, N=%d, alpha=%g, h=%g', ...
                  K, N, alpha, h));

    legend('Numerical PSD', 'Theoretical PSD');
    grid on;

endfunction

