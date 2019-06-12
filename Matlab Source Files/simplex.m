%Dual-implex method using "linprog"
function [obj,x_star]=simplex()
load('ecoli_core_model.mat');
options_sim = optimoptions('linprog','Algorithm','dual-simplex');
A=model.S; A=full(A); % Creating a full matrix from a sparse matrix
b = model.b; c = model.c; 
lb = model.lb; ub = model.ub;
x0=lb;
tic;
x=linprog(-c,A,b,A,b,lb,ub,x0,options_sim);
toc;
fprintf("Dual simplex solution: \n");
obj=c'*x;
x_star=x; 
end