% PARAFAC functionality
addpath("..\N-way-shell\Scripts\"); % own scripts
addpath("..\N-way-shell\N-way toolbox/"); % from Rasmus Bro

%%
% Load data
df = readmatrix("./Data/salivaCounts.csv");
taxonomy = readmatrix("./Data/taxonomy_fixed.csv", FileType="text", OutputType="string");
sampleInfo = readmatrix("./Data/salivaSampleMeta.csv", OutputType="string");

% Remove timepoint 9 months from the data
mask = sampleInfo(:,24) ~= "9";
df = df(mask,:);
sampleInfo = sampleInfo(mask,:);

% Separate TM and TW samples
salivaMaskTM = sampleInfo(:,7) == "TM";
salivaMaskTW = sampleInfo(:,7) == "TW";

salivaTM = df(salivaMaskTM,:);
salivaTW = df(salivaMaskTW,:);

%%
% Filter based on sparsity per group per niche
sparsity_saliva_TM = sum(salivaTM==0, 1) / size(salivaTM,1);
sparsity_saliva_TW = sum(salivaTW==0, 1) / size(salivaTW,1);

% subplot(2,1,1); histogram(sparsity_saliva_TM);
% subplot(2,1,2); histogram(sparsity_saliva_TW);

% Saliva
threshold = 0.5;
featureSelection1 = sparsity_saliva_TM <= threshold;
featureSelection2 = sparsity_saliva_TW <= threshold;
featureSelection = featureSelection1 | featureSelection2;

saliva_filtered = df(:,featureSelection);
taxa_saliva_filtered = taxonomy(featureSelection,:);

%%
% CLR
saliva_clr = transformCLR(saliva_filtered);

%%
% Make into cube
keepIndividuals = true;
[saliva_cube, salivaSubjectMeta, salivaConditionMeta] = rawDataToCube(saliva_clr, sampleInfo(:,6), str2double(sampleInfo(:,24)), keepIndividuals);

%%
% Center and scale
saliva_cnt = centerData(saliva_cube, 1);
saliva_cnt_scl = scaleData(saliva_cnt, 2);

%%
% Prepare metadata and remove outliers
saliva_subjectMeta_filtered = sortrows(unique(sampleInfo(:,[6 7]), "rows"), 1);
saliva_featureMeta_filtered = taxa_saliva_filtered;
saliva_timeMeta_filtered = salivaConditionMeta;

saliva_cnt_scl_filtered = saliva_cnt_scl;

%%
% Initialize options
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
const = [0 0 0];

%%
% Run PARAFACs
path_start = "./test_run/Figures/";
maxComponents=4;
numReps=25;
maxIterations = 20;
metaData = {saliva_subjectMeta_filtered, saliva_featureMeta_filtered, saliva_timeMeta_filtered};
resort = [true true false];
legendIndex = [2 4 0];
numColsPerLegend = [2 3 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
balancedJackKnifing = true;
subjectGroupCol = 2;

[salivaModels, salivaCons, salivaVarExps, salivaBoots, salivaBootVarExps, salivaTuckers] = quickReport(saliva_cnt_scl_filtered, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, "GOH TRANS bootstrapped saliva", path_start+"GOH_TRANS_saliva");

%%
% Dump
path_start = "./test_run/Dump/";

dump(salivaModels, salivaCons, salivaVarExps, salivaBoots, salivaBootVarExps, salivaTuckers, path_start, "GOHTRANS_saliva");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.
salivaNumFactors = 2;

saliva_choice = find(salivaCons{salivaNumFactors}==max(salivaCons{salivaNumFactors}));
salivaModel = pickModel(salivaModels, salivaNumFactors, saliva_choice);

%%
% Save models
model_path = "./test_run/PARAFAC models/";

salivaAnnotatedModel = annotateModel(saliva_cnt_scl_filtered, salivaModel, metaData);
savePARAFAC(saliva_cnt_scl_filtered, salivaModel, salivaAnnotatedModel, model_path + "GOHTRANS_saliva");

%%
% Plot PARAFAC model - saliva
resort = [true true false];
legendIndex = [4 6 0];
numColsPerLegend = [2 3 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
path_start = "./test_run/Figures/";

plotPARAFAC4(salivaAnnotatedModel, salivaNumFactors, salivaVarExps, saliva_choice, resort, legendIndex, numColsPerLegend, xlabels, titles, "PARAFAC GOHTRANS saliva", path_start + "PARAFAC_GOHTRANS_saliva.jpg");
