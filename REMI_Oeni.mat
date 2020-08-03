% This script was used to integrate transcriptome and metabolome data for acidic stress conditions to make a context-specific model. 
% This script was modified from tutorial of REMI methods
% Ref.: Pandey V, Hadadi N, Hatzimanikatis V. 2019. Enhanced flux prediction by integrating relative expression 
% and relative metabolite abundance into thermodynamically consistent metabolic models. PLoS Comput Biol 15:e1007036.

addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128')); % add your CPLEX path into matlab environment
initCobraToolbox

% load iQY486 model
load('C:\Users\Administrator\Documents\Tencent Files\296658517\FileRecv\SD-2a_all_open2'); 

% Add relative expression constraint 
% (1) load relative expression data
delimiterIn = '\t';
[A,delimiterOut]=importdata('C:\Users\Administrator\Desktop\result_3h.txt');
genedata = importdata('C:\Users\Administrator\Desktop\result_3h.txt', delimiterIn);
genename = genedata.textdata(:,1);
ratios = genedata.data(:,1);
up_cutoff=2; % fold-change>2 was regraded as up regulated genes
down_cutoff=0.5;% fold-change < 1/2 was regraded as down regulated genes 
% we need to identify genes which are up and down regulated
% input argument : 1) model
%                  2) gene names
%                  3) gene ratios
%                  4) and opertaor (evalute based on GPR)
%                  5) or operator
% length of argument 2 and 3 should be same
[rxns,ratio]=evaluateGPR(iQY486,genename,ratios,@min,@max);
% find up- and down- regulated reactions
indUP=find(ratio>up_cutoff);
ratioUP=ratio(indUP);
indDOWN=find(ratio<down_cutoff);
ratioDOWN=ratio(indDOWN);
regRxnRatio=[ratioUP;ratioDOWN];
regRxns=[rxns(indUP)';rxns(indDOWN)'];

% avoid numerical error (more than 100 fold is taken only 100)
regRxnRatio(regRxnRatio>100)=100;
regRxnRatio(regRxnRatio<1/100)=1/100;

% if we want to add relative constraint into TFA model then we need to add
% net flux variable to the model using 'addNetFluxVariablesNEW'
% and if one want to add into FBA model  then evalute scripts given below
mTmp=iQY486;
mTmp.rev=ones(numel(mTmp.rev),1);
use_m=addUseVariablesDH(mTmp);
netM=addNetFluxVariablesNEW(use_m);
netM1=addNetFluxVariablesNEW(use_m);
netM2=addNetFluxVariablesNEW(use_m);
% To force the growth consistent with prediction adding transcriptome data
% torlerance: 0.1 for no convergence was found in sometimes
netM2.var_lb(find(netM2.f))= 0.269308320451827;
netM2.var_ub(find(netM2.f))= 0.279308320451827;
netM1.var_lb(find(netM1.f))= 1.05839899362023;
netM1.var_ub(find(netM1.f))= 1.06839899362023;

%% We are going to add constraints for only relative expression
% now we add relative expression constarint 
% input argument : 1) model1 represents condition 1
%                  2) model2 represents condition2 
%                  3) regulated rxns 
%                  4) regulated reaction rations

% now we are build gex Model:  which integeartes relative gene expression
[gex]=addRelConsExpression(netM1,netM2,regRxns,regRxnRatio)
sol_gex=solveTFBAmodel(gex)

% maximum consistency score (MCS) will be objective value of the solution.
MCS=sol_gex.val; % this is the maximum consistency score of GeX model 

%% ALTERNATIVE ANALYSIS on a given MCS
% enumerate alternatives on expression variabels or for each MCS score
% there can be many alternative set of constraints and if one want to
% enumerate those alternatives one should use follwing script


path_save='C:\Users\Administrator\Documents\Tencent Files\296658517\FileRecv\new_3h.mat' % This path will save all alternative results
% all binary variables are stored as model.relExp.forB
%                  3) will be model.relExp.forB

% we want to find alternative at maximum consistency then we need to add
% maximum consistency to the lowerbound
coM2=gex;
coM2.var_lb(end)=MCS; % this will force maximum consistency can not be 
time=20;
% input argument : 1) numsol is number of alternatives
%                  2) model
%                  3) index for binary variables 
%                  4) path for saving the resultes
%                  5) time for solver
findAltCombi(time,coM2,coM2.relExp.forB,path_save);

%% ADD RELATIVE METABOLITE on M and Gex Model
%
% first load relative data
% load('/Users/vikash/Desktop/REMI/data/relMetdata.mat')
delimiterIn = '\t';
[A,delimiterOut]=importdata('C:\Users\Administrator\Desktop\metabolite_for_REMI.txt');
metdata = importdata('C:\Users\Administrator\Desktop\metabolite_for_REMI.txt', delimiterIn);
metname = metdata.textdata(:,1);
mratio = metdata.data(:,1);

indUP=find(mratio>up_cutoff);
indDOWN=find(mratio<down_cutoff);

regMetRatio=[mratio(indUP);mratio(indDOWN)];
regMets=[metname(indUP);metname(indDOWN)];
%% We are going to add constraints for only relative metabolites (M model)
% now we add relative expression constarint 
% input argument : 1) model1 represents condition 1
%                  2) model2 represents condition2 
%                  3) regulated mets 
%                  4) regulated metratios
%                  5) false (for only metaboliet), true(if expression is alreay intgerated)
[M]=addRelMetabolite(netM,netM,regMets,regMetRatio,false)
sol_M=solveTFBAmodel(M) 

% you can force consistency by putting MCS in the lower bound. 
% coM2.var_lb(end)=MCS; % this will force maximum consistency can not be 


% find alternatives
path_save='/Users/vikash/Desktop/REMI/data/AltSolOnlyMet.mat'

findAltCombi(2,M,M.metB,path_save,time);

%% We want to add relative expression and metabolites together (GexM model)

[gexM]=addRelMetabolite(gex,gex,regMets,regMetRatio,true);
sol=solveTFBAmodel(gexM)
% find alternatives
path_save='C:\Users\Administrator\Desktop\ALTERNATIVE_for_REMI3h_7.12.mat'
comIdx=[gexM.relExp.forB;gexM.metB];
findAltCombi(49,gexM,comIdx,path_save);

%% force model for one alternative solution

modelTMP=gexM;
sol1=solveTFBAmodel(modelTMP)
index1=comIdx(find(sol.x(comIdx)>0.98));
index0=comIdx(find(sol.x(comIdx)<0.08));
modelTMP.var_lb(index1)=1;
modelTMP.var_ub(index0)=0;
sol2=solveTFBAmodel(modelTMP);




