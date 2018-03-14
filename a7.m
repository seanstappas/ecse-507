function_names = ["Rosenbrock", "f(x) (n = 10)", "f(x) (n = 100)"];
functions = {f_rosenbrock() a7_function(10) a7_function(100)};
x0s = {[-1.2; 1] a7_function_x0(10) a7_function_x0(100)};

epsilon = 1e-6;
beta = 0.5;
sigma = 1e-4;
rho = 1e-8;
p = 2.1;
c1 = 1e-2;
c2 = 1;

for i = 1:3
    f = functions{i};
    
    x0 = x0s{i};
    [x_k, k] = newton_inexact(f, x0, epsilon, c1, c2, p, sigma, beta, rho);

    disp(function_names(i));
    disp("Solution x:");
    disp(x_k);
    disp("Solution f(x):");
    disp("f(x) = " + f(x_k));
    disp("Number of iterations:");
    disp(k);
    disp("-----------------");
end
