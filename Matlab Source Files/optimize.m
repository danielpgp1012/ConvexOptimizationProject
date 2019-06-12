function [obj,x_star]=optimize(model)
m = length(model.mets); 
n = length(model.rxns); 
p = 4;
b = model.b; A = model.S;
c = model.c; 
cvx_begin
    variable x(n)
    maximize( c'*(x))
    subject to
        A*x == b
        x >= model.lb
        x <= model.ub
cvx_end 
obj = cvx_optval;
x_star = x;




        
        
    