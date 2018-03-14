function f = a7_function(x, n)
f = 0;
for i = 1:n
    F_i = x(i) - 1;
    f = f + F_i ^ 2;
end
F_n1 = 0;
for j = 1:n
    F_n1 = F_n1 + j * (x(j) - 1);
end

f = f + F_n1 ^ 2;
F_n2 = F_n1^2;
f = f + F_n2^ 2;

end