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

cReaction_names=findCarbonRxns(substratesModel,1); %find substrates that contain
% AT LEAST 1 carbon
cReactions = findRxnIDs(model, cReaction_names);

%% Sample Optimization with succinate uptake of -10 g/h

biomass_idx = find(contains(model.rxns,"Biomass_Ecoli_core_N(w/GAM)-Nmet2"));
succinate_idx = find(contains(model.rxns,'EX_succ(e)'));
oxygen_idx = find(contains(model.rxns,'EX_o2(e)'));
model.lb(cReactions) = 0; % set consumption of everything to 0
model.lb(succinate_idx) = -10;
fprintf("LP solution: \n");
tic;
[obj,x] = optimize(model);
toc;

%% Interior point method
%x = model.lb;
%func = @(x) -model.c'*x;
%options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',2000);
%S=model.S; % Extracting the 
%S=full(S); % Creating a full matrix from a sparse matrix
fprintf("IP solution: \n");
[ipm_obj,ipm_x]=interior_point(model); % Solving the problem using IPM

%% Simplex method
fprintf("Simplex solution: \n");
tic;
[sim_obj,sim_x]=simplex(model); % Solving the problem using simplex method
toc;

%% Linprog


%% change oxygen and succinate
%model.ub(cReactions) = 1000;
UpOxygen=linspace(0,30,40); %oxygen vector 
UpSuccinate = linspace(0,30,40); % uptake Succinate
Objective_s=zeros(length(UpSuccinate),length(UpOxygen)); %initialize biomass matrix 
for i=1:length(UpOxygen) %nested for loop looking at biomass output at each glucose and oxygen combination
    model.lb(oxygen_idx) = -UpOxygen(i);
    %model.ub(oxygen_idx) = -UpOxygen(i);
    for j=1:length(UpSuccinate)
        model.lb(succinate_idx)=-UpSuccinate(j);
        %model.ub(succinate_idx)=-UpSuccinate(j);
        [sol_s,]=interior_point(model);
        Objective_s(i,j)=sol_s;
    end
end

fig6=figure;
surfl(UpSuccinate,UpOxygen,Objective_s); %3D plot
xlabel('Succinate Uptake [mmol/gDW/h]')
zlabel('Biomass [mmol/gDW/h]')
ylabel('Oxygen Uptake [mmol/gDW/h]');
