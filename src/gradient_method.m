function [x, k] = gradient_method(f, x0, epsilon)
beta = 0.5;
sigma = 1e-4;

k = 0;
x = x0;
grad_val = gradest(f, x).';
n = norm(grad_val);
while n > epsilon
    d = -grad_val;
    
    % Armijo update
    t = armijo(f, x, sigma, grad_val, d, beta);
    
    x = x + t*d;
    k = k + 1;
    grad_val = gradest(f, x).';
    n = norm(grad_val);
end
end