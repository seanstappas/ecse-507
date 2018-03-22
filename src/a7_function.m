function f = a7_function(n)
f = @(x) 0;
for i = 1:n
    F_i = @(x) x(i) - 1;
    f = @(x) f(x) + F_i(x) ^ 2;
end
F_n1 = @(x) 0;
for j = 1:n
    F_n1 = @(x) F_n1(x) + j * (x(j) - 1);
end

f = @(x) f(x) + F_n1(x) ^ 2;
F_n2 = @(x) F_n1(x)^2;
f = @(x) f(x) + F_n2(x)^ 2;

end