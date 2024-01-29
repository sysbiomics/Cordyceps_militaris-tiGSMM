%% Exploring light-exposed metabolic responses of Cordyceps 
% militaris through transcriptome-integrated genome-scale 
% modeling: 
%% The process was carried through the following steps:
% STEP 1: TEMPLATE MODEL PREPARATION
% STEP 2: GAP-FILLING AND ADD NEWLY-IDENTIFIED REACTIONS
% STEP 3: ADD BIOMASS REACTIONS
% STEP 4: MODEL VALIDATION
%% WORKSPACE:
cd 'C:\github\panyawarin'
% Initialize COBRA toolbox
initCobraToolbox;
%setRavenSolver('cobra');

% Load the model
model=importExcelModel('model/xlsx/iPS1474_retrofitted.xlsx')
%% STEP4 : MODEL VALIDATION
 %%% 4.1: light-programming condition
     % The growth simulation on glucose cultures grown under Light condition
    
ValidateModel_light = model;
ValidateModel_light = setParam(ValidateModel_light,'lb',{'bmOUT','cordycepinOUT','caroteniod'},[0 0 0]);
ValidateModel_light = setParam(ValidateModel_light,'eq',{'matp'},1);
ValidateModel_light = setParam(ValidateModel_light,'eq',{'R09782_c'},0);% constaint torulene to main caroteniod pathway
ValidateModel_light = setParam(ValidateModel_light,'ub',{'bmOUT','cordycepinOUT','caroteniod'},[1000 1000 1000]);
ValidateModel_light = setParam(ValidateModel_light,'obj',{'bmOUT'},1);
solL1 = solveLP(ValidateModel_light,1);

gluIN_L = setParam(ValidateModel_light,'ub',{'glcIN'},0.1593);%observed from experiment light conditions (Thananusak et al., 2020)
gluIN_L = setParam(gluIN_L,'ub',{'sucIN'},0); 
sol_LG = solveLP(gluIN_L,1);
sol_LG.f 

sucIN_L = setParam(ValidateModel_light,'ub',{'sucIN'},0.0845);%observed from experiment light conditions (Thananusak et al., 2020)
sucIN_L = setParam(sucIN_L,'ub',{'glcIN'},0); 
sol_LS = solveLP(sucIN_L,1);
sol_LS.f 

%sugars = {'glcIN' 'sucIN'};
%uptake = [0.1593,0.0845]; 
%sumax(11) = solL1;
%printFluxes(ValidateModel_light, sumax(11).x,true);

%for i = 1:numel(uptake)
   % model=setParam(ValidateModel_light,'ub',sugars,[0 0]);
   % model=setParam(model,'ub',sugars(i),uptake(i));
   % sugar = sugars(i);
   % solL = solveLP(model,1);
   % sumax(i) = solL;
   % umax = num2str(solL.f*-1);
   % fprintf([char(sugar) '\t' num2str(solL.f*-24) ' per day' '\n']);
%end

% 4.2: Synthesis capabilities of Caroteniod in Light condition

finalValidateModel10 = setParam(gluIN_L,'obj',{'cordycepinOUT'},1);
finalValidateModel11 = setParam(gluIN_L,'obj',{'exc sphinganine'},1);
finalValidateModel12 = setParam(gluIN_L,'obj',{'exc phytosphingosine'},1);
finalValidateModel13 = setParam(gluIN_L,'obj',{'exc sphingosine'},1);
finalValidateModel14 = setParam(gluIN_L,'obj',{'caroteniod'},1);

sol10 = solveLP(finalValidateModel10,1);
fprintf(['production rate of cordycepin' '\t' num2str(sol10.f*-1) ' mmol/g DW/h' '\n']);
sol11 = solveLP(finalValidateModel11,1);
fprintf(['production rate of sphinganine' '\t' num2str(sol11.f*-1) ' mmol/g DW/h' '\n']);
sol12 = solveLP(finalValidateModel12,1);
fprintf(['production rate of phytosphingosine' '\t' num2str(sol12.f*-1) ' mmol/g DW/h' '\n']);
sol13 = solveLP(finalValidateModel13,1);
fprintf(['production rate of sphingosine' '\t' num2str(sol13.f*-1) ' mmol/g DW/h' '\n']);
sol14 = solveLP(finalValidateModel14,1);
fprintf(['production rate of caroteniod' '\t' num2str(sol14.f*-1) ' mmol/g DW/h' '\n']);

%%%%%%%%%%%%%%%%%%%%%%%%
