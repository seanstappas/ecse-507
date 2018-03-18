function display_solution(f, x, k, x_ks, x_range, y_range)
disp("Solution x:");
disp(x);
disp("Solution f(x):");
disp(f(x));
disp("Number of iterations:");
disp(k);
fcontour(f, x_range, y_range);