syms x1 x2 x3 t1 t

f_rosenbrock = 100 * (x2 - x1^2)^2 + (1 - x1)^2;

f_bateman = x3 * (exp(-x1 * t1) - exp(-x2 * t1));
bateman_ti_yi = [[15; 0.038] [25; 0.085] [35; 0.1] [45; 0.103] ...
    [55; 0.093] [65; 0.095] [75; 0.088] [85; 0.08] [105; 0.073] ...
    [245; 0.038] [305; 0.028] [365; 0.02]];
f_bateman_min = 0;
for ti_yi = bateman_ti_yi
    ti = ti_yi(1);
    yi = ti_yi(2);
    f_bateman_min = f_bateman_min + (subs(f_bateman, t1, ti) - yi)^2;
end
f_bateman_min = f_bateman_min / 2;

epsilon = 1e-6;
beta = 0.5;
sigma = 1e-4;
rho = 1e-8;
p = 2.1;
c1 = 1e-2;
c2 = 1;

functions = [f_rosenbrock f_bateman_min];
function_names = ["Rosenbrock", "Bateman"];
x0s = [[-1.2; 1; 0] [0.05; 0.1; 0.4]];

for i = [1 2]
    f = functions(i);
    
    grad = gradient(f, [x1; x2; x3]);
    n = norm(grad);
    grad2 = hessian(f, [x1; x2; x3]);

    f_func = matlabFunction(f, 'Vars', {[x1; x2; x3]});
    grad_func = matlabFunction(grad, 'Vars', {[x1; x2; x3]});
    grad2_func = matlabFunction(grad2, 'Vars', {[x1; x2; x3]});
    norm_func = matlabFunction(n, 'Vars', {[x1; x2; x3]});
    
    x = x0s(:,i);
    k = 0;
    while norm_func(x) > epsilon
        grad_val = grad_func(x);
        eta1 = c1 / (k + 1);
        eta2 = c2 * norm_func(x);
        eta = min([eta1 eta2]);
        
        % CG method (to solve inexact Newton equation)
        A = grad2_func(x);
        x_cg = [0; 0; 0]; % Starting point for d?
        b = -grad_val;
        g = A * x_cg - b;
        d_cg = -g;
        epsilon_cg = eta * norm_func(x);
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
        f_val = f_func(x);
        rhs = sigma * grad_val.' * d;
        while f_func(x + t*d) > f_val + t*rhs
            l = l + 1;
            t = beta ^ l;
        end
        
        x = x + t*d;
        k = k + 1;
    end

    disp(function_names(i));
    disp("Solution x:");
    disp("x1 = " + x(1));
    disp("x2 = " + x(2));
    if i == 2
        disp("x3 = " + x(3));
    end
    disp("Solution f(x):");
    disp("f(x) = " + f_func(x));
    disp("Number of iterations:");
    disp(k);
    disp("-----------------");
end

f_bateman_optimal = subs(f_bateman, [x1; x2; x3], x);
fplot(f_bateman_optimal);
xlim([0 365]);
ylim([0 0.1]);