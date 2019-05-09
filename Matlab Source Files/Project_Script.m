%% init
initCobraToolbox
model = load('ecoli_core_model.mat');
%%
m = length(model.mets); n = length(model.rxns); p = 4;
b = model.b; A = model.S;
c = model.c; 
biomass_idx = find(contains(model.rxns,"Biomass_Ecoli_core_N(w/GAM)-Nmet2"));
ATP_idx = find(contains(model.rxns,'ATPM'));
cvx_begin
    variable x(n)
    maximize( c'*x )
    subject to
        A*x == b
        x >= model.lb
        x <= model.ub
cvx_end 
