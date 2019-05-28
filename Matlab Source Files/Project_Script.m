%% init
initCobraToolbox
%%
%% Load model
load('ecoli_core_model.mat');
%% Set consumption of all carbon substrates to 0
[selExc,selUpt]=findExcRxns(model,0); %outputs exchange reactions

uptakes=model.rxns(selExc); %selects substrates from exchange reactions

substratesModel=extractSubNetwork(model,uptakes); %create submodel only
% with given substrates of substrates

cReactions=findCarbonRxns(substratesModel,1); %find substrates that contain
% AT LEAST 1 carbon
cReactions = findRxnIDs(model, cReactions);

%% Sample Optimization with succinate uptake of -10 g/h

biomass_idx = find(contains(model.rxns,"Biomass_Ecoli_core_N(w/GAM)-Nmet2"));
succinate_idx = find(contains(model.rxns,'EX_succ(e)'));
oxygen_idx = find(contains(model.rxns,'EX_o2(e)'));
model.lb(cReactions) = 0; % set consumption of everything to 0
model.lb(succinate_idx) = -10;
fprintf("LP solution: \n");
tic;
[obj,x] = optimize(model);
toc

%% Interior point method
func = @(x) -model.c'*x;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',300);
S=model.S; % Extracting the 
S=full(S); % Creating a full matrix from a sparse matrix
fprintf("IP solution: \n");
tic;
optimal=fmincon(func,x,[],[],S,model.b,model.lb,model.ub,[],options); % Solving the problem using IPM
toc;


%% change oxygen and succinate
succ_uptake = linspace(-30,0,30);
ox_uptake = linspace(-30,0,30);
biomass = zeros(length(ox_uptake),length(succ_uptake));


model=changeRxnBounds(model,cReactions,0,'b'); %set all carbon substrate uptakes to 0
model=changeRxnBounds(model,cReactions,1000,'u'); %set upper boundary to 1000
UpOxygen=linspace(0,30,40); %oxygen vector 
Objective_s=zeros(length(UpSuccinate),length(UpOxygen)); %initialize biomass matrix 

for i=1:length(UpOxygen) %nested for loop looking at biomass output at each glucose and oxygen combination
    model = changeRxnBounds(model,'EX_o2(e)',-UpOxygen(i),'b');
    for j=1:length(UpSuccinate)
        model = changeRxnBounds(model,'EX_succ(e)',-UpSuccinate(j),'b');
        sol_s=optimizeCbModel(model,'max');
        Objective_s(i,j)=sol_s.obj;
    end
end

fig6=figure;
surfl(UpSuccinate,UpOxygen,Objective_s); %3D plot
xlabel('Succinate Uptake [mmol/gDW/h]')
zlabel('Biomass [mmol/gDW/h]')
ylabel('Oxygen Uptake [mmol/gDW/h]');
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