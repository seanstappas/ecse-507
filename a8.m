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
disp("5. All methods (for comparison)");
method_index = input('Choose the minimization method. ');
switch method_index
    case 1
        disp("Gradient method");
        [x, k] = gradient_method(f, x0, epsilon);
        display_solution(f, x, k);
    case 2
        disp("Globalized Newton's method");
        [x, k] = newton_global(f, x0, epsilon);
        display_solution(f, x, k);
    case 3
        disp("Globalized BFGS method");
        [x, k] = bfgs_global(f, x0, epsilon);
        display_solution(f, x, k);
    case 4
        disp("Globalized inexact Newton's method");
        [x, k] = newton_inexact(f, x0, epsilon);
        display_solution(f, x, k);
    case 5
        disp("----------------------------------")
        disp("Gradient method");
        [x, k] = gradient_method(f, x0, epsilon);
        display_solution(f, x, k);
        disp("----------------------------------")
        disp("Globalized Newton's method");
        [x, k] = newton_global(f, x0, epsilon);
        display_solution(f, x, k);
        disp("----------------------------------")
        disp("Globalized BFGS method");
        [x, k] = bfgs_global(f, x0, epsilon);
        display_solution(f, x, k);
        disp("----------------------------------")
        disp("Globalized inexact Newton's method");
        [x, k] = newton_inexact(f, x0, epsilon);
        display_solution(f, x, k);
    otherwise
        disp('Invalid minimization method.');
end