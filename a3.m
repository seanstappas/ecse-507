tic;
syms x1 x2

f = 100*(x2 - x1^2)^2 + (1 - x1)^2;
grad = gradient(f, [x1, x2]);
n = norm(grad);

f_func = matlabFunction(f, 'Vars', {[x1; x2]});
grad_func = matlabFunction(grad, 'Vars', {[x1; x2]});
norm_func = matlabFunction(n, 'Vars', {[x1; x2]});

x = [-1.2; 1];

beta = 0.5;
sigma = 1e-4;
epsilon = 1e-4;
iter = 0;

while norm_func(x) > epsilon
    grad_val = grad_func(x);
    d = -grad_val;
    l = 0;
    t = 1;
    f_val = f_func(x);
    rhs = sigma * grad_val.' * d;
    while f_func(x + t*d) > f_val + t*rhs
        l = l + 1;
        t = beta ^ l;
    end
    x = x + t*d;
    iter = iter + 1;
end

toc;
disp("Final x: ")
disp(x)
disp("Final f(x): ")
disp(f_func(x))
disp("Number of iterations: ")
disp(iter)