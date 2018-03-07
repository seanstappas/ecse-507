x = sym('x', [1, 100]).';

f_rosenbrock = 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
f_bateman = 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;

epsilon = 1e-6;
beta = 0.5;
sigma = 1e-4;
rho = 1e-8;
p = 2.1;
c1 = 1e-2;
c2 = 1;

functions = [f_rosenbrock f_bateman];
function_names = ["Rosenbrock", "Bateman"];
x0s = [[[-1.2; 1; 0]; zeros(1, 100 - 3).'] ...
    [[0.05; 0.1; 0.4]; zeros(100 - 3, 1)]];

for i = [1 2]
    f = functions(i);
    
    grad = gradient(f, x);
    n = norm(grad);
    grad2 = hessian(f, x);

    f_func = matlabFunction(f, 'Vars', {x});
    grad_func = matlabFunction(grad, 'Vars', {x});
    grad2_func = matlabFunction(grad2, 'Vars', {x});
    norm_func = matlabFunction(n, 'Vars', {x});
    
    x_k = x0s(:,i);
    k = 0;
    while norm_func(x_k) > epsilon
        grad_val = grad_func(x_k);
        eta1 = c1 / (k + 1);
        eta2 = c2 * norm_func(x_k);
        eta = min([eta1 eta2]);
        
        % CG method (to solve inexact Newton equation)
        A = grad2_func(x_k);
        x_cg = zeros(100, 1); % Starting point for d?
        b = -grad_val;
        g = A * x_cg - b;
        d_cg = -g;
        epsilon_cg = eta * norm_func(x_k);
        while norm(g) > epsilon_cg
            g_norm_2 = norm(g)^2;
            t_cg = g_norm_2 / (d_cg' * A * d_cg);

            x_cg = x_cg + t_cg * d_cg;
            g = g + t_cg * A * d_cg;
            beta_cg = norm(g)^2 / g_norm_2;
            d_cg = -g + beta_cg * d_cg;
        end
        
        d = x_cg;
        if grad_val' * d > - rho * (norm(d)^p)
            d = -grad_val;
        end
        
        % Armijo update
        l = 0;
        t = 1;
        f_val = f_func(x_k);
        rhs = sigma * grad_val.' * d;
        while f_func(x_k + t*d) > f_val + t*rhs
            l = l + 1;
            t = beta ^ l;
        end
        
        x_k = x_k + t*d;
        k = k + 1;
    end

    disp(function_names(i));
    disp("Solution x:");
    disp("x1 = " + x_k(1));
    disp("x2 = " + x_k(2));
    if i == 2
        disp("x3 = " + x_k(3));
    end
    disp("Solution f(x):");
    disp("f(x) = " + f_func(x_k));
    disp("Number of iterations:");
    disp(k);
    disp("-----------------");
end
