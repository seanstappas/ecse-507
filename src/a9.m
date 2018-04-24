gamma = -sqrt(2);

fun = @(x) -(x(1) + 1)^2 - (x(2) + 1)^2;

lb = [-Inf,-Inf];
ub = [gamma,Inf];  

A = [];
b = [];
Aeq = [];
beq = [];  

x0 = [gamma,0];  

nonlcon = @circlecon;
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)   
