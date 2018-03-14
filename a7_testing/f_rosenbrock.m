function f_rosenbrock = f_rosenbrock()
f_rosenbrock = @(x) 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
end