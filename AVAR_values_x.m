function AVAR_values_x = AVAR_values_x(x, q)
# Take the matrix x and calculate the AVAR values of each row
# for q(1), q(2), ... q(M).#
# AVAR_values_x(j,i) = AVAR of x(j,:) for q(i).

    % Convert x[n] to y[n]
    y = diff(x, 1, 2);

    % Calculate AVAR using the previously defined function
    AVAR_values_x = AVAR_values_y(y, q);

endfunction
