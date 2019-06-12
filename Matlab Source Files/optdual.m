function [obj,y_star]=optdual(model)
m = length(model.mets); 
n = length(model.rxns); 
p = 4;
b = model.b; 
A = model.S;
c = model.c; 
lb = model.lb; % lower bound
ub = model.ub; % upper bound
cvx_begin
    variable v(m)
    variable l(n)
    variable g(n)
    minimize (-b'*(v)+ l'*lb - g'*ub)
    subject to 
        A'*v-l+c == -g 
        l >= zeros(n,1)
        g >= zeros(n,1)
cvx_end
obj = cvx_optval;
y_star = y;
