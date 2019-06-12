function [obj,y_star]=optdual(model)
m = length(model.mets); 
n = length(model.rxns); 
p = 4;
b = model.b; 
A = model.S;
c = model.c; 
low = model.lb;
high = model.ub;
cvx_begin
    variable y(m)
    variable l(n)
    variable g(n)
    maximize (-b'*(y)+ l'*low - g'*high)
    subject to 
        A'*y-l== -g 
        l >= zeros(n,1)
        g >= zeros(n,1)
cvx_end
obj = cvx_optval;
y_star = y;
