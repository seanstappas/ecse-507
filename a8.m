epsilon = 1e-6;

disp("1. Rosenbrock function");
disp("2. Himmelblau function");
function_index = input('Choose the function to be minimized. ');
switch function_index
    case 1
        f = f_rosenbrock();
        x0 = [-1.2; 1];
        x_range = [-1.5, 1.5, -0.5, 2];
    case 2
        f = f_himmelblau();
        x0 = [4; 4];
        x_range = [-6, 6, -6, 6];
    otherwise
        disp("Invalid function.");
        return;
end
disp(function_index);

disp("1. Gradient method");
disp("2. Globalized Newton's method");
disp("3. Globalized BFGS method");
disp("4. Globalized inexact Newton's method");
disp("5. Fletcher-Reeves method");
disp("6. All methods (for comparison)");
method_index = input('Choose the minimization method. ');
switch method_index
    case 1
        disp("Gradient method");
        [x, k, x_ks] = gradient_method(f, x0, epsilon);
        display_solution(f, x, k, x_ks, x_range);
    case 2
        disp("Globalized Newton's method");
        [x, k, x_ks] = newton_global(f, x0, epsilon);
        display_solution(f, x, k, x_ks, x_range);
    case 3
        disp("Globalized BFGS method");
        [x, k, x_ks] = bfgs_global(f, x0, epsilon);
        display_solution(f, x, k, x_ks, x_range);
    case 4
        disp("Globalized inexact Newton's method");
        [x, k, x_ks] = newton_inexact(f, x0, epsilon);
        display_solution(f, x, k, x_ks, x_range);
    case 5
        disp("Fletcher-Reeves method");
    case 6
        disp("----------------------------------")
        disp("Gradient method");
        [x, k, x_ks] = gradient_method(f, x0, epsilon);
        display_solution(f, x, k, x_ks, x_range);
        disp("----------------------------------")
        disp("Globalized Newton's method");
        [x, k, x_ks] = newton_global(f, x0, epsilon);
        display_solution(f, x, k, x_ks, x_range);
        disp("----------------------------------")
        disp("Globalized BFGS method");
        [x, k, x_ks] = bfgs_global(f, x0, epsilon);
        display_solution(f, x, k, x_ks, x_range);
        disp("----------------------------------")
        disp("Globalized inexact Newton's method");
        [x, k, x_ks] = newton_inexact(f, x0, epsilon);
        display_solution(f, x, k, x_ks, x_range);
        disp("----------------------------------")
        disp("Fletcher-Reeves method");
        disp("----------------------------------")
    otherwise
        disp('Invalid minimization method.');
end