function [x, k, x_ks] = newton_global(f, x0, epsilon)
rho = 1e-8;
p = 2.1;
beta = 0.5;
sigma = 1e-4;
k_max = 200;

k = 0;
x = x0;
x_ks = {x0};
grad_val = gradest(f, x).';
n = norm(grad_val);
while n > epsilon && k <= k_max
    hess = hessian(f, x);
    
    d = -grad_val;
    if cond(inv(hess)) > 1e-12
        d_soln = hess \ (-grad_val);
        if grad_val' * d_soln <= - rho * (norm(d_soln) ^ p)
            d = d_soln;
        end
    end
    
    % Armijo update
    t = armijo(f, x, sigma, grad_val, d, beta);
    
    x = x + t*d;
    x_ks{end + 1} = x;
    k = k + 1;
    grad_val = gradest(f, x).';
    n = norm(grad_val);
end
end