function x = cg(A, x0, b, epsilon)
x = x0;
g = A * x - b;
d_cg = -g;
while norm(g) > epsilon
    g_norm_2 = norm(g)^2;
    t_cg = g_norm_2 / (d_cg' * A * d_cg);

    x = x + t_cg * d_cg;
    g = g + t_cg * A * d_cg;
    beta_cg = norm(g)^2 / g_norm_2;
    d_cg = -g + beta_cg * d_cg;
end
end