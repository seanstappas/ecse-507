rho = 1e-8;
p = 2.1;
beta = 0.5;
sigma = 1e-4;
epsilon = 1e-6;
k_max = 200;

syms x1 x2

f_rosenbrock = 100 * (x2 - x1^2)^2 + (1 - x1)^2;
f_brown = (x1 - 1e6)^2 + (x2 - 2e-6)^2 + (x1 * x2 - 2)^2;
functions = [f_rosenbrock f_brown];
x0s = [[1.2; 1] [1; 1]];
function_names = ["Rosenbrock", "Brown"];

for i = [1 2]
    f = functions(i);
    
    grad = gradient(f, [x1; x2]);
    n = norm(grad);
    grad2 = hessian(f, [x1; x2]);

    f_func = matlabFunction(f, 'Vars', {[x1; x2]});
    grad_func = matlabFunction(grad, 'Vars', {[x1; x2]});
    grad2_func = matlabFunction(grad2, 'Vars', {[x1; x2]});
    norm_func = matlabFunction(n, 'Vars', {[x1; x2]});

    x = x0s(:,i);
    k = 0;
    while norm_func(x) > epsilon && k <= k_max
        grad_val = grad_func(x);
        grad2_val = grad2_func(x);
        d = -grad_val;
        if cond(inv(grad2_val)) > 1e-12
            d_soln = grad2_val \ (-grad_val);
            if grad_val' * d_soln <= - rho * (norm(d_soln) ^ p)
                d = d_soln;
            end
        end
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
    disp("Solution f(x):");
    disp("f(x) = " + f_func(x));
    disp("Number of iterations:");
    disp(k);
    disp("-----------------");
end
