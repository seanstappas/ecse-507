function [x, k] = newton_inexact(f, x0, epsilon)
beta = 0.5;
sigma = 1e-4;
rho = 1e-8;
p = 2.1;
c1 = 1e-2;
c2 = 1;

k = 0;
x = x0;
grad_val = gradest(f, x).';
n = norm(grad_val);
while n > epsilon
    eta = min([c1 / (k + 1), c2 * n]);

    % CG method (to solve inexact Newton equation)
    d = cg(hessian(f, x), -grad_val, -grad_val, eta * n);

    if grad_val' * d > - rho * (norm(d)^p)
        d = -grad_val;
    end

    % Armijo update
    t = armijo(f, x, sigma, grad_val, d, beta);

    x = x + t*d;
    k = k + 1;
    grad_val = gradest(f, x).';
    n = norm(grad_val);
end
end