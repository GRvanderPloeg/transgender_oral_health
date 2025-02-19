% PARAFAC functionality
addpath("..\N-way-shell\Scripts\"); % own scripts
addpath("..\N-way-shell\N-way toolbox/"); % from Rasmus Bro

%%
% Load data
df = readmatrix("./Data/20241209_cytokines.csv");
featureMeta = readmatrix("./Data/20241209_cytokines_featureMeta.csv", FileType="text", OutputType="string");
sampleInfo = readmatrix("./Data/20241209_cytokines_sampleMeta.csv", OutputType="string");

% Separate TM and TW samples
cytokineMaskTW = sampleInfo(:,2) == "TW";
df = df(cytokineMaskTW,:);
sampleInfo = sampleInfo(cytokineMaskTW,:);

% Prepare y
% Testosterone is 4
% Estradiol is 5
mask0 = (sampleInfo(:,3) == "0");
mask3 = (sampleInfo(:,3) == "3");
mask12 = (sampleInfo(:,3) == "12");

testosterone0 = sampleInfo(mask0, [1 5]);
testosterone0 = testosterone0(testosterone0(:,2) ~= "NA",:); % remove NAs
testosterone3 = sampleInfo(mask3, [1 5]);
testosterone3 = testosterone3(testosterone3(:,2) ~= "NA",:); % remove NAs
testosterone12 = sampleInfo(mask12, [1 5]);
testosterone12 = testosterone12(testosterone12(:,2) ~= "NA",:); % remove NAs

testosterone0 = sortrows(testosterone0,1);
testosterone3 = sortrows(testosterone3,1);
testosterone12 = sortrows(testosterone12,1);

subsetSubjects0 = ismember(sampleInfo(:,1), testosterone0(:,1));
subsetSubjects3 = ismember(sampleInfo(:,1), testosterone3(:,1));
subsetSubjects12 = ismember(sampleInfo(:,1), testosterone12(:,1));

% For 0m, 3m
usableSubjects = unique(sampleInfo(subsetSubjects0 & subsetSubjects3, 1));
testosterone = [usableSubjects testosterone0(ismember(testosterone0(:,1), usableSubjects),2) testosterone3(ismember(testosterone3(:,1), usableSubjects),2)];

% For 0m, 12m
%usableSubjects = unique(sampleInfo(subsetSubjects0 & subsetSubjects12, 23));
%testosterone = [usableSubjects testosterone0(ismember(testosterone0(:,1), usableSubjects),2) testosterone12(ismember(testosterone12(:,1), usableSubjects),2)];

% Finally
testosterone = sortrows(testosterone,1);

% Where to save output
path = "./20241212_Cytokines_TW/";

%%
% Log transform
cytokines_log = log(df);

%%
% Make into cube
keepIndividuals = true;
[cytokine_cube, cytokineSubjectMeta, cytokineConditionMeta] = rawDataToCube(cytokines_log, sampleInfo(:,1), str2double(sampleInfo(:,3)), keepIndividuals);

%%
% Remove subjects with missing testosterone measurements
mask = ismember(cytokineSubjectMeta, usableSubjects);
cytokine_cube = cytokine_cube(mask,:,:);
cytokineSubjectMeta = cytokineSubjectMeta(mask,:);

%%
% Center and scale
cytokine_cnt = centerData(cytokine_cube, 1);
cytokine_cnt_scl = scaleData(cytokine_cnt, 2);

%%
% Prepare metadata and remove outliers
cytokine_subjectMeta_filtered = sortrows(unique(sampleInfo(:,[1 2]), "rows"), 1);
cytokine_subjectMeta_filtered = cytokine_subjectMeta_filtered(mask,:);
cytokine_featureMeta_filtered = featureMeta;
cytokine_timeMeta_filtered = cytokineConditionMeta;

cytokine_cnt_scl_filtered = cytokine_cnt_scl;

%%
% Center testosterone
testosterone = str2double(testosterone(:,2:end));
testosterone = testosterone(:,2) - testosterone(:,1);
testosterone = testosterone - mean(testosterone);

%%
% Initialize options
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
const = [0 0 0];

%%
% Determine correct number of NPLS components per niche
maxComponents = 10;
[XValResult_cytokine,~] = ncrossreg(cytokine_cnt_scl_filtered, testosterone, maxComponents, 0);

%%
% Plot RMSEP of CV and save it
%set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
%plot(1:maxComponents, XValResult.RMSEP); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP"); yyaxis right; plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); ylim([0 100]); ylabel("Variance explained (%)");
%saveas(gcf, path+"RMSEP_CV.jpg");
%close();
XValResult = XValResult_cytokine;

% NEW APPROACH
df = [XValResult.RMSEP' XValResult.Percent.Xexp(:,1) XValResult.Percent.Yexp(:,1) XValResult.PRESS'];

%plot(1:maxComponents, XValResult.RMSEP); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP"); yyaxis right; plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); ylim([0 100]); ylabel("Variance explained (%)");
%plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("Variance explained (%)"); plot(1:maxComponents, XValResult.Percent.Yexp(:,1));

%colororder({'r', 'b', 'm'});
%plot(1:maxComponents, df(:,[2 3])); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("Variance explained (%)"); yyaxis right; plot(1:maxComponents, df(:,1)); ylabel("RMSEP"); legend("VarExpX", "VarExpY", "RMSEP");

% ALTERNATIVE
plot(1:maxComponents, df(:,[2 3])); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("Variance explained (%)"); legend("X", "Y");
saveas(gcf, path+"RMSEP_CV_new_varExps.jpg");
close();
plot(1:maxComponents, df(:,1)); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP");
saveas(gcf, path+"RMSEP_CV_new_RMSEP.jpg");
close();
plot(1:maxComponents, df(:,4)); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("PRESS");
saveas(gcf, path+"RMSEP_CV_new_PRESS.jpg");
close();


%%
% Find optimal number of components based on RMSEP CV
%[bestRMSEP_saliva, numComponents_saliva] = min(XValResult_saliva.RMSEP);

% OVERRIDE
numComponents_cytokines = 1;

%%
% Run NPLS
[Xfactors_cytokine,Yfactors_cytokine,Core_cytokine,B_cytokine,ypred_cytokine,ssx_cytokine,ssy_cytokine,reg_cytokine] = npls(cytokine_cnt_scl_filtered, testosterone, numComponents_cytokines);

%%
% Save NPLS models
metaData = {cytokine_subjectMeta_filtered, cytokine_featureMeta_filtered, cytokine_timeMeta_filtered};
cytokineAnnotatedModel = annotateModel(cytokine_cnt_scl_filtered, Xfactors_cytokine, metaData);
savePARAFAC(cytokine_cnt_scl_filtered, Xfactors_cytokine, cytokineAnnotatedModel, path + "GOHTRANS_cytokine");

%%
% Save ypred of NPLS models
%y_saliva = [ypred_saliva saliva_subjectMeta_filtered testosterone];
y_cytokine = [ypred_cytokine(:,:,1) cytokine_subjectMeta_filtered testosterone];
writematrix(y_cytokine, path + "Cytokine_" + numComponents_cytokines + "_ypred.csv");

%%
% Save coefficients of NPLS models
for i=1:size(reg_cytokine, 2)
    writematrix(reg_cytokine{i}, path + "Cytokine_" + numComponents_cytokines + "_coeff_" + i + ".csv");
end

%%
% Save crossvalidated coefficients of the NPLS models
[coeff_mean_cytokine, coeff_std_cytokine] = CV_coeff_NPLS(cytokine_cnt_scl_filtered, testosterone, numComponents_cytokines, 1);
writematrix(coeff_mean_cytokine, path + "Cytokine_" + numComponents_cytokines + "_CVcoeff_means.csv");
writematrix(coeff_mean_cytokine, path + "Cytokine_" + numComponents_cytokines + "_CVcoeff_stds.csv");