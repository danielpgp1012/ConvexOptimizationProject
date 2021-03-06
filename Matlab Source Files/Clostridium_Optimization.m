%% Load the COBRA Toolbox
addpath("C:\Users\danie\Documents\GitHub\cobratoolbox")
addpath('C:\tomlab\tomlab20181122_1')

initCobraToolbox
%% convert model
%model=sbmlTestModelToMat(pwd, pwd);

%% load model
load ('iCbu641.mat');

%% Change objective
model=changeObjective(model,'Biomass',1);


%% Test biomass at given glycerol uptake
%[selExc,selUpt]=findExcRxns(model,0); %outputs exchange reactions

%uptakes=model.rxns(selUpt); %selects substrates from exchange reactions

%substratesModel=extractSubNetwork(model,uptakes); %create submodel only
% with given substrates of substrates

%cReactions=findCarbonRxns(substratesModel,1); %find substrates that contain

load ('cNames') ; %Filter carbon-containing substrates
indices=find(ismember(model.rxnNames,cNames));
cReactions = model.rxns(indices);
model=changeRxnBounds(model,cReactions,0,'l'); %change boundaries
glycerol_rxn=cReactions(find(ismember(cNames,'Glycerol exchange')));
model=changeRxnBounds(model,glycerol_rxn,-10,'l');

%% find solution
solFBA=optimizeCbModel(model);

%% for loop to plot curv
uptake_flux=linspace(0,100);
biomass_LP=zeros(length(uptake_flux),1);
PDO_LP = biomass_LP;
x0=zeros(length(model.rxns),length(uptake_flux));
PDO_idx = find(ismember(model.rxnNames,'Propane-1,3-diol exchange '));
for i=1:length(uptake_flux)
    model=changeRxnBounds(model,glycerol_rxn,-uptake_flux(i),'l');
    solFBA=optimizeCbModel(model);
    biomass_LP(i)=solFBA.f;
    PDO_LP(i) = solFBA.full(PDO_idx);
    x0(:,i)=solFBA.full;
end
%% Plot LP
figure 
subplot(1,2,1)
plot(uptake_flux,biomass_LP)
xlim([0 100])
xticks(0:20:100)
ylabel('Growth rate (h^{-1})')
xlabel('Glycerol Uptake flux mmol/gDWh')
ylim([0 1.2])
subplot(1,2,2)
plot (uptake_flux,PDO_LP)
ylabel('PDO flux mmol/gDWh')
xlabel('Glycerol Uptake flux mmol/gDWh')
xticks(0:20:100)
%% block formic and hydrogen production

BlockRxnNames={'Hydrogen exchange ','Formate exchange '};
BlockRxnIndices=[find(ismember(model.rxnNames,'Hydrogen exchange ')),find(ismember(model.rxnNames,'Formate exchange '))];
BlockRxnCodes=model.rxns(BlockRxnIndices);
 model=changeRxnBounds(model,BlockRxnCodes,0,'u');
 

biomass_LP2=zeros(length(uptake_flux),1);
PDO_LP2 = biomass_LP2;
for i=1:length(uptake_flux)
    model=changeRxnBounds(model,glycerol_rxn,-uptake_flux(i),'l');
    solFBA=optimizeCbModel(model);
    biomass_LP2(i)=solFBA.f;
    PDO_LP2(i) = solFBA.full(PDO_idx);
end

%% Plot H2 and formic acid block
figure
plot (uptake_flux,biomass_LP2)
ylabel('Growth rate (h^{-1})')
xlabel('Glycerol Uptake flux mmol/gDWh')
 
%% Load NL model
load ('model_NL.mat')
%% Find Enzyme Reactions
T=strfind(model_NL.rxnNames,'ase');
indeces_enzymes=zeros(length(T),1);
k=1;
for i=1:length(T)
    if ~isempty(T{i})
        indeces_enzymes(k)=i;
        k=k+1;
    end
end
indeces_enzymes(k:end)=[];
model_NL.enzyme_rxn=indeces_enzymes;
model_NL.b=zeros(701,1);



%% Conopt Solver Setup

acetate_idx=find(ismember(model_NL.rxnNames,'Acetate exchange '));

C_L = -100;
C_U = 0;
model_NL.lb(acetate_idx)= 0;
model_NL=changeRxnBounds(model_NL,glycerol_rxn,0,'b');
x_0 = x0(:,30);
NLP = conAssign('objfun', 'objGradient', [], [], model_NL.lb, model_NL.ub, 'NLP', x_0,0, [], ...
                           model_NL.A, model_NL.b, model_NL.b,'acetate_constrain',[],[],[],C_L,C_U);

glycerol_index = find(ismember(model_NL.rxns,glycerol_rxn));
uptake_flux=linspace(0,90);
biomass_NL1=zeros(length(uptake_flux),1);
PDO_NL1 = biomass_NL1;
x0=x0*0;

%% For loop
for i=1:length(uptake_flux)
   
    NLP.x_L(glycerol_index) = -uptake_flux(i);
    NLP.x_U(glycerol_index) = -uptake_flux(i);
    solConopt = tomRun('conopt',NLP,1);
    biomass_NL1(i)=solConopt.x_k(17);
    PDO_NL1(i) = solConopt.x_k(PDO_idx);
    NLP.x_0 = solConopt.x_k;
    x0(:,i) = solConopt.x_k;
end
clear NLP;
%% Second NL Function





%% minimize atp requirement

ATP_idx = find(model_NL.A(1,:)==1);
% NLP = conAssign('objfun_2', 'objGradient_2', [], [], model_NL.lb, model_NL.ub, 'NLP', x_0,0, [], ...
%                            model_NL.A, model_NL.b, model_NL.b,'acetate_constrain',[],[],[],C_L,C_U);
 NLP = conAssign('objfun_2', 'objGradient_2', [], [], model_NL.lb, model_NL.ub, 'NLP', x_0,0, [], ...
                            model_NL.A, model_NL.b, model_NL.b);


uptake_flux_NL2=linspace(0,90);
biomass_NL2=zeros(length(uptake_flux_NL2),1);
PDO_NL2 = biomass_NL2; 

for i=1:length(uptake_flux_NL2)
   
    NLP.x_L(glycerol_index) = -uptake_flux_NL2(i);
    NLP.x_U(glycerol_index) = -uptake_flux_NL2(i);
    solConopt = tomRun('conopt',NLP,1);
    biomass_NL2(i)=solConopt.x_k(17);
    PDO_NL2(i) = solConopt.x_k(PDO_idx);
    NLP.x_0 = solConopt.x_k;
    %NLP.x_0 = x0(i);
end

%%
%% Comparison plot
uptake_flux=linspace(0,100);


figure 
subplot(1,2,1)

plot(uptake_flux,biomass_LP,'b',uptake_flux,biomass_LP2,'c',...
    uptake_flux_NL2,biomass_NL1,'r',...
    uptake_flux_NL2,biomass_NL2,'g');
leg=legend('$\mu$','$\mu\,no\,H_2$','$\frac{\mu}{\sum(v_i^2)}$','$\frac{\mu}{(w*\sum(vi^2)+(1-w)v_{ATP}^2)}$');
set(leg,'Interpreter','latex');
ylabel('Growth rate (h^{-1})')
xlabel('Glycerol Uptake flux mmol/gDWh')
ylim([0 1.2])
xlim([0 100])
xticks(0:20:100)
subplot(1,2,2)
plot(uptake_flux,PDO_LP,'b',uptake_flux,PDO_LP2,'c',...
uptake_flux_NL2,PDO_NL1,'r',uptake_flux_NL2,PDO_NL2,'g')
ylabel('PDO flux mmol/gDWh')
xlabel('Glycerol Uptake flux mmol/gDWh')
xticks(0:20:100)
ylim([-2 70])
%%
figure 
plot (uptake_flux(2:end),transpose(PDO_NL2(2:end))./uptake_flux(2:end))
xlabel('Uptake flux Glycerol')
ylabel({'Y_{PDO/S}'},'Interpreter','latex')

%% Second Graph
if (~isempty(NLP))
    clear NLP;
end
%NLP = conAssign('objfun_2', 'objGradient_2', [], [], model_NL.lb, model_NL.ub, 'NLP', x_0,0, [], ...
%                            model_NL.A, model_NL.b, model_NL.b);
NLP = conAssign('objfun', 'objGradient', [], [], model_NL.lb, model_NL.ub, 'NLP', x_0,0, [], ...
                           model_NL.A, model_NL.b, model_NL.b,'acetate_constrain',[],[],[],C_L,C_U);
%NLP = conAssign('objfun', 'objGradient', [], [], model_NL.lb, model_NL.ub, 'NLP', x_0,0, [], ...
%                           model_NL.A, model_NL.b, model_NL.b);

                       
glucose_uptake = linspace(0,14,25);
glycerol_uptake = linspace(0,40,25);

glucose_index = find(ismember(model_NL.rxnNames,'D-Glucose exchange'));

yield_PDO = zeros(length(glucose_uptake),length(glycerol_uptake));
oxygen_index = find(ismember(model_NL.rxnNames,'Oxygen exchange'));
NLP.x_L(oxygen_index) = 0;
NLP.x_U(oxygen_index) = 0;


for i=1:length(glucose_uptake)
   NLP.x_L(glucose_index) = -glucose_uptake(i);
   NLP.x_U(glucose_index) = -glucose_uptake(i);

   
    for j = 1:length(glycerol_uptake)
        NLP.x_L(glycerol_index) = -glycerol_uptake(j);
        NLP.x_U(glycerol_index) = -glycerol_uptake(j);
       
        
        %if (i>1 ||  j>1)
        if (j>1)
            solConopt = tomRun('conopt',NLP,1);
            %yield_PDO(i,j) = solConopt.x_k(PDO_idx)/(-solConopt.x_k(glucose_index)-solConopt.x_k(glycerol_index));
            yield_PDO(i,j) = solConopt.x_k(PDO_idx)/(-solConopt.x_k(glycerol_index));
            yield_PDO(1,1)=0;
            if (~contains(solConopt.ExitText,'Infeasible'))
            NLP.x_0 = solConopt.x_k;
            end
           
        end
        
    end
    
end

%% Plot
[GLC,GLY] = meshgrid(glucose_uptake,glycerol_uptake);
yield_PDO(yield_PDO>1)=1;
surf(GLY,GLC,transpose(yield_PDO))
title('Co-Fermentation yield of PDO')
ylabel('Glucose Uptake Flux (mmol/gDW h)')
xlabel('Glycerol Uptake Flux (mmol/gDW h)')
zlabel({'$Y_(PDO/S)$'},'Interpreter','latex')

%% Left plot

Glu_Gly_ratio = GLC./GLY;
Glu_Gly_ratio = Glu_Gly_ratio(2:end,:);
Glu_Gly_ratio = reshape(Glu_Gly_ratio.',1,[]);
yield_PDO_flattened = reshape((yield_PDO(2:end,:)).',1,[]);

x_bins = [0,0.053,0.104,0.238,1.149];
tol = [1e-3,5e-3,1e-2,5e-2,5e-2,1e-1];
y_bins = zeros(1,length(x_bins));
std_dev = zeros(1,length(x_bins));
for i=1:length(x_bins)
   y_vals = yield_PDO_flattened(abs(Glu_Gly_ratio-x_bins(i))<tol(i));
   y_bins(i) = mean(y_vals);
   std_dev(i) = std(y_vals);
end
figure 
plot (x_bins,y_bins,'*')
ylim([0 1.1])
xlim([-0.1 1.2])
errorbar(x_bins,y_bins,std_dev);



