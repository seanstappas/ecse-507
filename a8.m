disp("1. Rosenbrock function");
disp("2. Himmelblau function");
function_index = input('Choose the function to be minimized. ');
switch function_index
    case 1
        f = f_rosenbrock();
    case 2
        f = f_himmelblau();
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
method_index = input('Choose the minimization method. ');
switch method_index
    case 1
        % a3
    case 2
        % a5
    case 3
        % a6
    case 4
        % a7
    case 5
    otherwise
        disp('Invalid minimization method.');
        return;
end
disp(method_index);