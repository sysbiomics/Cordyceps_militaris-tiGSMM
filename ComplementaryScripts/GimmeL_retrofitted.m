%% Exploring  light-exposed metabolic responses of Cordyceps 
% militaris through transcriptome-integrated genome-scale 
% modeling: 
%% iPS474:
initCobraToolbox;
%WORKSPACE
cd 'C:\github\panyawarin'
model=importExcelModel('model/xlsx/iPS1474_retrofitted.xlsx')
%setRavenSolver('cobra');
%% iPS1474_light
ValidateModel_light = model;
ValidateModel_light = setParam(ValidateModel_light,'lb',{'bmOUT','cordycepinOUT','carotenoid'},[0 0 0]);
ValidateModel_light = setParam(ValidateModel_light,'eq',{'matp'},1);
ValidateModel_light = setParam(ValidateModel_light,'eq',{'R09782_c'},0);% constaint torulene to main caroteniod pathway
ValidateModel_light = setParam(ValidateModel_light,'ub',{'bmOUT','cordycepinOUT','carotenoid'},[1000 1000 1000]);
ValidateModel_light = setParam(ValidateModel_light,'obj',{'bmOUT'},1);

%Check the iPS1474 model
solL1 = solveLP(ValidateModel_light,1);
sugars = {'glcIN' 'sucIN'};
uptake = [0.1593,0.0845]; %observed from experiment(Thananusak et al., 2020)
sumax(11) = solL1;
printFluxes(ValidateModel_light, sumax(11).x,true);

for i = 1:numel(uptake)
    model_L=setParam(ValidateModel_light,'ub',sugars,0);
    model_L=setParam(model_L,'ub',sugars(i),uptake(i));
    sugar = sugars(i);
    solL = solveLP(model_L,1);
    sumax(i) = solL;
    umax = num2str(solL.f*-1);
    fprintf([num2str(solL.f*-1) ' per hour' '\n']);
end


%% GIMME algorithm 

model_iL = model;
%Converting RAVEN structure to COBRA structure
ModelCb_L = ravenCobraWrapper(model_iL);
changeCobraSolver('glpk');

ModelCb_L = changeRxnBounds(ModelCb_L,'matp',1,'u');
ModelCb_L = changeRxnBounds(ModelCb_L,'matp',1,'l');
ModelCb_L = changeRxnBounds(ModelCb_L,'bmOUT',0,'l');
ModelCb_L = changeRxnBounds(ModelCb_L,'R09782_c',0,'l');% constaint torulene to main caroteniod pathway
ModelCb_L = changeRxnBounds(ModelCb_L,'R09782_c',0,'u');
ModelCb_L = changeRxnBounds(ModelCb_L,'bmOUT',1000,'u');
ModelCb_L = changeRxnBounds(ModelCb_L,'cordycepinOUT',0,'l');
ModelCb_L = changeRxnBounds(ModelCb_L,'cordycepinOUT',1000,'u');
ModelCb_L = changeRxnBounds(ModelCb_L,'carotenoid',1000,'u');
ModelCb_L = changeRxnBounds(ModelCb_L,'carotenoid',0,'l');

%growth on Light_Sucrose
ModelCb_LS = ModelCb_L;
ModelCb_LS = changeRxnBounds(ModelCb_LS,'glcIN',0,'u'); 
ModelCb_LS = changeRxnBounds(ModelCb_LS,'glcIN',0,'l');
ModelCb_LS = changeRxnBounds(ModelCb_LS,'sucIN',0.0845,'u');%growth on Light_Sucrose 0.0122(0.012287134546773)
ModelCb_LS = changeRxnBounds(ModelCb_LS,'sucIN',0,'l'); 
checkObjective(ModelCb_LS);

solutionCb_LS = optimizeCbModel(ModelCb_LS,'max');
solutionCb_LS.f 
printFluxes(ModelCb_L, solutionCb_LS.x,true);

%growth on Light_Glucose
ModelCb_LG = ModelCb_L; 
ModelCb_LG = changeRxnBounds(ModelCb_LG,'sucIN',0,'u');
ModelCb_LG = changeRxnBounds(ModelCb_LG,'sucIN',0,'l');
ModelCb_LG = changeRxnBounds(ModelCb_LG,'glcIN',0.1593,'u');
ModelCb_LG = changeRxnBounds(ModelCb_LG,'glcIN',0,'l');
checkObjective(ModelCb_LG);

solutionCb_LG = optimizeCbModel(ModelCb_LG,'max'); 
solutionCb_LG.f 
printFluxes(ModelCb_LG, solutionCb_LG.x,true);

%% Gimme Light_Sucrose
%GIMME algorithm
   %%% Maps gene expression into reaction
parsedGPR_LS = GPRparser(ModelCb_LS);

%import Gene expression data 1.Gene 2.Exp
expressionCol_LS = selectGeneFromGPR(ModelCb_LS, Gene_LS, Exp_LS, parsedGPR_LS);

 %Simulates fluxes using GIMME
    %%%Gimme of Light_Sucrose condition
LSgmodel_50 = GIMME(ModelCb_LS, expressionCol_LS, 44.82); %Thresholds 50th percentile 

%OptimizeCbModel_Light sucrose
Sol_LSgmodel50 = optimizeCbModel(LSgmodel_50); 
Sol_LSgmodel50.f
printFluxes(LSgmodel_50, Sol_LSgmodel50.x,true);

%exportForGit(LSgmodel_50,'LSgmodel_50','',{'mat'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gimme Light_Glucose
%GIMME algorithm
   %%% Maps gene expression into reaction
parsedGPR_LG = GPRparser(ModelCb_LG);

%import Gene expression data 1.Gene 2.Exp
expressionCol_LG = selectGeneFromGPR(ModelCb_LG, Gene_LG, Exp_LG, parsedGPR_LG);

%Simulates fluxes using GIMME
    %%%Gimme of Light_Glucose condition
LGgmodel_50 = GIMME(ModelCb_LG, expressionCol_LG, 46.7681); % Thresholds 50th percentile

%optimizeCbModel_Light Glucose
Sol_LGgmodel50 = optimizeCbModel(LGgmodel_50);
Sol_LGgmodel50.f
printFluxes(LGgmodel_50, Sol_LGgmodel50.x,true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



