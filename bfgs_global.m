function [x, k, x_ks] = bfgs_global(f, x0, epsilon)
sigma = 1e-4;
rho = 0.9;
H_0 = eye(size(x0, 1));
t_0 = 1;
gamma = 2;

k = 0;
x = x0;
x_ks = {x0};
H = H_0;
grad_val = gradest(f, x).';
n = norm(grad_val);
while n > epsilon
    d = H \ (-grad_val);
    
    phi = @(t) f(x + t * d);
    phi_prime_0 = gradest(phi, 0);
    psi = @(t) phi(t) - phi(0) - sigma * t * phi_prime_0;

    t_i = t_0;

    stop = false;
    while ~stop && psi(t_i) < 0
        if gradest(phi, t_i) >= rho * phi_prime_0
            stop = true;
        else
            t_i = gamma * t_i;
        end
    end

    a = 0;
    b = t_i;
    while ~stop
        if psi(t_i) >= 0
            b = t_i;
        elseif gradest(phi, t_i) < rho * phi_prime_0
            a = t_i;
        else
            stop = true;
        end
        t_i = a + (b - a) / 2;
    end

    
    x_old = x;
    x = x + t_i*d;
    x_ks{end + 1} = x;
    grad_val_new = gradest(f, x).';
    s = x - x_old;
    y = grad_val_new - grad_val;
    H = H + (y * y') / (y' * s) - (H * (s * s') * H) / (s' * H * s);
    k = k + 1;
    
    k = k + 1;
    grad_val = grad_val_new;
    n = norm(grad_val);
end
end