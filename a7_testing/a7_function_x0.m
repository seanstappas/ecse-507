function x0 = a7_function_x0(n)
x0 = zeros(n, 1);
for i = 1:n
    x0(i) = 1 - i / n;
end
end