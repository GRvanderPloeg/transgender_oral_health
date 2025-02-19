% PARAFAC functionality
addpath("..\N-way-shell\Scripts\"); % own scripts
addpath("..\N-way-shell\N-way toolbox/"); % from Rasmus Bro

%%
% Load data
df = readmatrix("./Data/tongueCounts.csv");
taxonomy = readmatrix("./Data/taxonomy_fixed.csv", FileType="text", OutputType="string");
sampleInfo = readmatrix("./Data/tongueSampleMeta.csv", OutputType="string");

% Remove timepoint 9 months from the data
mask = sampleInfo(:,13) ~= "9";
df = df(mask,:);
sampleInfo = sampleInfo(mask,:);

% Separate TM and TW samples
tongueMaskTM = sampleInfo(:,7) == "TM";
tongueMaskTW = sampleInfo(:,7) == "TW";
tongueTM = df(tongueMaskTM,:);
tongueTW = df(tongueMaskTW,:);

%%
% Filter based on sparsity per group per niche
sparsity_tongue_TM = sum(tongueTM==0, 1) / size(tongueTM,1);
sparsity_tongue_TW = sum(tongueTW==0, 1) / size(tongueTW,1);

% subplot(2,1,1); histogram(sparsity_tongue_TM);
% subplot(2,1,2); histogram(sparsity_tongue_TW);

% Tongue
threshold = 0.5;
featureSelection1 = sparsity_tongue_TM <= threshold;
featureSelection2 = sparsity_tongue_TW <= threshold;
featureSelection = featureSelection1 | featureSelection2;

tongue_filtered = df(:,featureSelection);
taxa_tongue_filtered = taxonomy(featureSelection,:);

%%
% CLR
tongue_clr = transformCLR(tongue_filtered);

%%
% Make into cube
keepIndividuals = true;
[tongue_cube, tongueSubjectMeta, tongueConditionMeta] = rawDataToCube(tongue_clr, str2double(sampleInfo(:,6)), str2double(sampleInfo(:,13)), keepIndividuals);

%%
% Center and scale
tongue_cnt = centerData(tongue_cube, 1);
tongue_cnt_scl = scaleData(tongue_cnt, 2);

%%
% Prepare metadata and remove outliers
tongue_subjectMeta_filtered = sortrows(unique(sampleInfo(:,[6 7]), "rows"), 1);
tongue_featureMeta_filtered = taxa_tongue_filtered;
tongue_timeMeta_filtered = tongueConditionMeta;

tongue_cnt_scl_filtered = tongue_cnt_scl;

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
metaData = {tongue_subjectMeta_filtered, tongue_featureMeta_filtered, tongue_timeMeta_filtered};
resort = [true true false];
legendIndex = [2 4 0];
numColsPerLegend = [2 3 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
balancedJackKnifing = true;
subjectGroupCol = 2;

[tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers] = quickReport(tongue_cnt_scl_filtered, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, "GOH TRANS bootstrapped tongue", path_start+"GOH_TRANS_tongue");

%%
% Dump
path_start = "./test_run/Dump/";

dump(tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers, path_start, "GOHTRANS_tongue");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.
tongueNumFactors = 1;

tongue_choice = find(tongueCons{tongueNumFactors}==max(tongueCons{tongueNumFactors}));
tongueModel = pickModel(tongueModels, tongueNumFactors, tongue_choice);

%%
% Save models
model_path = "./test_run/PARAFAC models/";

tongueAnnotatedModel = annotateModel(tongue_cnt_scl_filtered, tongueModel, metaData);

savePARAFAC(tongue_cnt_scl_filtered, tongueModel, tongueAnnotatedModel, model_path + "GOHTRANS_tongue");

%%
% Plot PARAFAC model - tongue
resort = [true true false];
legendIndex = [3 5 0];
numColsPerLegend = [2 3 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
path_start = "./test_run/Figures/";

plotPARAFAC4(tongueAnnotatedModel, tongueNumFactors, tongueVarExps, tongue_choice, resort, legendIndex, numColsPerLegend, xlabels, titles, "PARAFAC GOHTRANS tongue", path_start + "PARAFAC_GOHTRANS_tongue.jpg");