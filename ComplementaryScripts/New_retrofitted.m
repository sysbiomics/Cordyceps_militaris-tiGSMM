%% Exploring light-exposed metabolic responses of Cordyceps 
% militaris through transcriptome-integrated genome-scale 
% modeling: 
%% The process was carried through the following steps:
% STEP 1: TEMPLATE MODEL PREPARATION
% STEP 2: GAP-FILLING AND ADD NEWLY-IDENTIFIED REACTIONS
% STEP 3: ADD BIOMASS REACTIONS
% STEP 4: MODEL VALIDATION
%% WORKSPACE
cd 'C:\github\panyawarin'
% Initialize COBRA toolbox
initCobraToolbox;
%setRavenSolver('cobra');
%% STEP 1: TEMPLATE MODEL PREPARATION
% RETROFITTED THE GSMM-MODEL
 % The Genome-scale metabolic model(GSMM)of C. militaris,iNR1329 and iPC1469 was used as a template to retrofitted genome-scale
 % metabolic network reconstructions.
 
% Import the model
load('ComplementaryData/revised_iPC1469.mat');% Import template the model
iPC1469Model = model;
iPC1469Cobra = ravenCobraWrapper(iPC1469Model);%Converting RAVEN structure to COBRA
iPC1469Raven = ravenCobraWrapper(iPC1469Cobra);%Converting COBRA structure to RAVEN

%Check the iPC1469 model
ValidateModel = iPC1469Raven;
ValidateModel = setParam(ValidateModel,'lb',{'bmOUT','cordycepinOUT'},[0 0]);
ValidateModel = setParam(ValidateModel,'eq',{'matp'},1);
ValidateModel = setParam(ValidateModel,'ub',{'bmOUT','cordycepinOUT'},[1000 1000]);
ValidateModel = setParam(ValidateModel,'obj',{'bmOUT'},1);
sol1 = solveLP(ValidateModel,1);

sugars = {'glcIN' 'fruIN' 'arabIN' 'xylIN' 'sucIN'};
uptake = [0.1448,0.1472,0.1074,0.0681,0.0815]; %observed from experiment(Raethong  et al., 2020) 
sumax(11) = sol1;
printFluxes(ValidateModel, sumax(11).x,true);


for i = 1:numel(uptake)
    model=setParam(ValidateModel,'ub',sugars,[0 0 0 0 0]);
    model=setParam(model,'ub',sugars(i),uptake(i));
    sugar = sugars(i);
    sol = solveLP(model,1);
    sumax(i) = sol;
    umax = num2str(sol.f*-1);
    fprintf([char(sugar) '\t' num2str(sol.f*-24) ' per hour' '\n']);
end

%% STEP 2: GAP-FILLING AND ADD NEWLY-IDENTIFIED REACTIONS

% 2.1 update grRules
% Manually curated with gene and protein information in C. militaris genome using the changeGeneAssociation function.

update_grmodel = iPC1469Raven
[~, newGrRules]=xlsread('ComplementaryData/supplementary_new.xlsx','updated_grRules');
rxnID = newGrRules(2:end,1);
geneAssoc = newGrRules(2:end,2);

for i = 1:numel(rxnID);
    rxn = rxnID(i);
    new = geneAssoc{i};
    update_grmodel = changeGeneAssoc(update_grmodel,rxn,new,true);
end

update_grmodel=deleteUnusedGenes(update_grmodel);
solgr = solveLP(update_grmodel,1);
solgr.f

% 2.2 remove old biomass reaction
rxnsToRemove = {'biomass',};
update_grmodel = removeReactions(update_grmodel,rxnsToRemove,true,...
            true,true);
           
        
% 2.3 Addition of new metabolites to the templates network 
  % New metabolites listed in new_Metabolites were introduced to the templates network by addMets function.

[~, newMets]=xlsread('ComplementaryData/supplementary_new.xlsx','new_Metabolites');
metsToAdd = struct();
metsToAdd.mets = newMets(2:end,1);
metsToAdd.metNames = newMets(2:end,2);
metsToAdd.metFormulas = newMets(2:end,3);
metsToAdd.compartments = 'c';
metsToAdd.inchis = newMets(2:end,9);
addMetsModel=addMets(update_grmodel,metsToAdd);
addMetsModel = contractModel(addMetsModel);

addMetsModel.equations = constructEquations(addMetsModel);
solmet = solveLP(addMetsModel,1);
solmet.f

% 2.3 Addition of new reactions to the templates network 
  % New reactions listed in new_Reactions were introduced to the templates network by addRxns function.

[~, SheetS]=xlsread('ComplementaryData/supplementary_new.xlsx','new_Reactions');
newRxns = struct();
newRxns.rxns = SheetS(2:end,1);
newRxns.rxnNames = SheetS(2:end,2);
newRxns.equations = SheetS(2:end,3);
newRxns.eccodes = SheetS(2:end,4);
newRxns.grRules = SheetS(2:end,5);
newRxns.rxnNotes = SheetS(2:end,13);

addRxnModel = addRxns(addMetsModel,newRxns,3,'',true,true);
addRxnModel.equations = constructEquations(addRxnModel);
addRxnModel = sortModel(addRxnModel);

model_new1 = contractModel(addRxnModel);
model_new1.lb(isinf(-addRxnModel.lb))=-1000; 
model_new1.ub(isinf(addRxnModel.ub))=1000;

model_new1.equations = constructEquations(model_new1);
solrxn = solveLP(model_new1,1);
printFluxes(model_new1, solrxn.x,true);
iPC1469plus = model_new1;

%% STEP 3 : ADD BIOMASS REACTIONS
%  The biomass reaction was added to the model, which was adopted from the template-based model of iPC1469.

[~, SheetS]=xlsread('ComplementaryData/supplementary_new.xlsx','new_Biomass');
add_Biomass = struct();
add_Biomass.rxns = SheetS(2:end,1);
add_Biomass.rxnNames = SheetS(2:end,2);
add_Biomass.equations = SheetS(2:end,3);
add_Biomass.rxnNotes = SheetS(2:end,13);

model_bm = addRxns(iPC1469plus,add_Biomass,3,'',true,false);
model_bm.equations = constructEquations(model_bm);
model_bm = sortModel(model_bm);
model_bm1 = contractModel(model_bm);

[a, b]=ismember(model_bm1.rxns,add_Biomass.rxns);
I=find(a);
model_bm1.rev(I)=0;
reducedModel1=model_bm1;

for i = 1:numel(add_Biomass.rxns)
    reducedModel1=setParam(reducedModel1,'lb',add_Biomass.rxns(i),[0]);
    reducedModel1=setParam(reducedModel1,'ub',add_Biomass.rxns(i),[1000]);
end

reducedModel1.equations = constructEquations(reducedModel1);
solrxn = solveLP(reducedModel1,1);
printFluxes(reducedModel1, solrxn.x,true);
iPC1469plus = reducedModel1;

%exportForGit(iPC1469plus,'iPC1469pl','',{'xlsx'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPS1474 = iPC1469plus;
iPS1474.id = 'iPS1474';
iPS1474.name = 'iPS1474';
iPS1474.description = 'light-exposed metabolic responses of Cordyceps militaris';
iPS1474.annotation = [];
printModelStats(iPS1474);

%%%%%%%%%%%%%%%%
%exportForGit(iPS1474,'iPS1474_retrofitted','',{'mat', 'txt', 'xlsx', 'xml', 'yml'});

