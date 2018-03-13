function [x, k] = newton_inexact(x0, norm_func, epsilon, grad_func, c1, c2, ...
    grad2_func, p, f_func, sigma, beta, rho)

k = 0;
x = x0;
while norm_func(x) > epsilon
    grad_val = grad_func(x);
    eta = min([c1 / (k + 1), c2 * norm_func(x)]);

    % CG method (to solve inexact Newton equation)
    d = cg(grad2_func(x), -grad_val, -grad_val, ...
        eta * norm_func(x));

    if grad_val' * d > - rho * (norm(d)^p)
        d = -grad_val;
    end

    % Armijo update
    t = armijo(f_func, x, sigma, grad_val, d, beta);

    x = x + t*d;
    k = k + 1;
    disp(k);
end
end