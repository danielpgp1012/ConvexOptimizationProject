%Interior point method using "linprog"
function [obj,x_star]=interior_point()
load('ecoli_core_model.mat');
options_ipm = optimoptions('linprog','Algorithm','interior-point');
A=model.S; A=full(A); % Creating a dense matrix from a sparse matrix
b = model.b; c = model.c; 
lb = model.lb; ub = model.ub;
x0=lb;
func = @(x) -c'*x;
fprintf("IPM solution: \n");
tic;
x_star=linprog(-c,A,b,A,b,lb,ub,x0,options_ipm); % Solving the problem using IPM
toc;
obj=c'*x_star;
end