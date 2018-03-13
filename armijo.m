function t = armijo(f_func, x_k, sigma, grad_val, d, beta)
l = 0;
t = 1;
f_val = f_func(x_k);
rhs = sigma * grad_val.' * d;
while f_func(x_k + t*d) > f_val + t*rhs
    l = l + 1;
    t = beta ^ l;
end
end