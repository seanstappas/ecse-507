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
sigma = 1e-4;
rho = 0.9;
H_0 = eye(3);

gamma = 2;
t_0 = 1;
tau_1 = 1/4;
tau_2 = 1/4;

functions = [f_rosenbrock f_bateman_min];
function_names = ["Rosenbrock", "Bateman"];
x0s = [[1.2; 1; 0] [0.05; 0.1; 0.4]];

for i = [1 2]
    f = functions(i);
    
    grad = gradient(f, [x1; x2; x3]);
    n = norm(grad);

    f_func = matlabFunction(f, 'Vars', {[x1; x2; x3]});
    grad_func = matlabFunction(grad, 'Vars', {[x1; x2; x3]});
    norm_func = matlabFunction(n, 'Vars', {[x1; x2; x3]});
    
    H = H_0;
    x = x0s(:,i);
    k = 0;
    while norm_func(x) > epsilon
        grad_val_old = grad_func(x);
        d = H \ (-grad_val_old);
        
        phi = subs(f, [x1; x2; x3], x + t * d);
        phi_prime = diff(phi);
        phi_func = matlabFunction(phi, 'Vars', t);
        phi_prime_func = matlabFunction(phi_prime, 'Vars', t);
        phi_prime_0 = phi_prime_func(0);
        
        psi = phi - phi_func(0) - sigma * phi_prime_0;
        psi_func = matlabFunction(psi, 'Vars', t);
        
        t_i = t_0;
        
        stop = false;
        while ~stop && psi_func(t_i) < 0
            if phi_prime_func(t_i) >= rho * phi_prime_0
                stop = true;
            else
                t_i = gamma * t_i;
            end
        end
        
        a = 0;
        b = t_i;
        while ~stop
            if psi_func(t_i) >= 0
                b = t_i;
            elseif phi_prime_func(t_i) < rho * phi_prime_0
                a = t_i;
            else
                stop = true;
            end
            t_i = a + (b - a) / 2;
        end
        
        x_old = x;
        x = x + t_i*d;
        s = x - x_old;
        y = grad_func(x) - grad_val_old;
        H = H + (y * y') / (y' * s) - (H * (s * s') * H) / (s' * H * s);
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