function [x, k] = newton_inexact(f, x0, epsilon, c1, c2, p, sigma, beta, rho)

k = 0;
x = x0;
[grad_val, ~] = gradest(f, x);
grad_val = grad_val.';
n = norm(grad_val);
while n > epsilon
    eta = min([c1 / (k + 1), c2 * n]);

    % CG method (to solve inexact Newton equation)
    [hess, ~] = hessian(f, x);
    d = cg(hess, -grad_val, -grad_val, eta * n);

    if grad_val' * d > - rho * (norm(d)^p)
        d = -grad_val;
    end

    % Armijo update
    t = armijo(f, x, sigma, grad_val, d, beta);

    x = x + t*d;
    k = k + 1;
    [grad_val, ~] = gradest(f, x);
    grad_val = grad_val.';
    n = norm(grad_val);
end
end