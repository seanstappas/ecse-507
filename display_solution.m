function display_solution(f, x, k, x_ks, x_range)
disp("Solution x:");
disp(x);
disp("Solution f(x):");
disp(f(x));
disp("Number of iterations:");
disp(k);
f_vector = @(x1,x2) f([x1, x2]);
fcontour(f_vector, x_range);
colorbar