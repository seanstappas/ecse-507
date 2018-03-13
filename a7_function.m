function f = a7_function(x, n)
f = 0;
for i = 1:n
    f = f + (x(i) - 1)^2;
end
F_n1 = 0;
for j = 1:n
    F_n1 = F_n1 + j * (x(j) - 1);
end

f = f + F_n1 ^ 2;
f = f + F_n1 ^ 4;
end