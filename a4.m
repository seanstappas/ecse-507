delta_values = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6];

for delta = delta_values
    disp("=======================")
    disp("delta value:")
    disp(delta)
    Q = diag([1 1 1 delta]);
    c = [1 1 1 1];
    gamma = 0;

    syms x1 x2 x3 x4
    f = 0.5 * [x1 x2 x3 x4] * Q * [x1 x2 x3 x4].' + c * [x1 x2 x3 x4].' + gamma;
    grad = gradient(f, [x1 x2 x3 x4]);
    n = norm(grad);

    f_func = matlabFunction(f, 'Vars', {[x1 x2 x3 x4]});
    grad_func = matlabFunction(grad, 'Vars', {[x1 x2 x3 x4]});
    norm_func = matlabFunction(n, 'Vars', {[x1 x2 x3 x4]});

    x = [1 2 3 4];
    iter = 0;
    while norm_func(x) > 1e-5
        g = grad_func(x);
        d = -g;
        t = norm(g)^2 / (g' * Q * g);
        x = x + t * d.';
        iter = iter + 1;
    end

    disp("Exact minimization rule");
    disp("Solution x:");
    disp(x);
    disp("Solution f(x):");
    disp(f_func(x));
    disp("Number of iterations:");
    disp(iter);

    x = [1 2 3 4];
    y = x;
    alpha = 1;
    iter = 0;
    L = norm(Q);
    while norm_func(x) > 1e-5
        x_old = x;
        alpha_old = alpha;

        x = y - (1 / L) * grad_func(y).';
        alpha = (1 + sqrt(1 + 4 * alpha^2)) / 2;
        y = x + (alpha_old - 1) * (x - x_old) / alpha;
        iter = iter + 1;
    end

    disp("Accelerated gradient method");
    disp("Solution x:");
    disp(x);
    disp("Solution f(x):");
    disp(f_func(x));
    disp("Number of iterations:");
    disp(iter);
end