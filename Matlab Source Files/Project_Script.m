%% init
initCobraToolbox
load('ecoli_core_model.mat');
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

%% Clostridium Model

load ('iCbu641.mat');
PDO_idx = find(ismember(model.rxnNames,'Propane-1,3-diol exchange '));
biomass_idx = find(ismember(model.rxnNames,'Biomass Synthesis '));
alphas = linspace(0.5,0.8,100);
c = model.c;
n = length(model.rxns);
biomass_list = zeros(length(alphas),1);
PDO_list = zeros(length(alphas),1);

%% Filter input substrates
[selExc,selUpt]=findExcRxns(model,0); %outputs exchange reactions

uptakes=model.rxns(selUpt); %selects substrates from exchange reactions

substratesModel=extractSubNetwork(model,uptakes); %create submodel only
% with given substrates of substrates

cReactions=findCarbonRxns(substratesModel,1); %find substrates that contain
cindices=find(ismember(model.rxns,cReactions));
glycerol_idx=find(ismember(model.rxnNames,'Glycerol exchange'));
model.lb(cindices)=zeros(length(cindices),1);
model.lb(glycerol_idx) = -10;
%% 

for i=1:length(alphas)
    alphas(i) = 1;
    c(biomass_idx)=alphas(i);
    c(PDO_idx) = 1-alphas(i);
    
    cvx_begin
        variable x(n)
        maximize( c'*x )
        subject to
            model.S*x == model.b
            x >= model.lb
            x <= model.ub
    cvx_end
    
    PDO_list(i) = x(PDO_idx);
    biomass_list(i) = x(biomass_idx);
end