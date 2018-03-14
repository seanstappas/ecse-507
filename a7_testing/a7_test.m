rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
% The gradient should be zero (within floating point noise)
[grad,err] = gradest(rosen, [1, 1])

x = [-1.2; 1; 0];
f = f_rosenbrock();
[grad_val, ~] = gradest(f, x)