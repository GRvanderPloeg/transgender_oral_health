% PARAFAC functionality
addpath("..\N-way-shell\Scripts\"); % own scripts
addpath("..\N-way-shell\N-way toolbox/"); % from Rasmus Bro

%%
% Load data
df = readmatrix("./Data/20241209_cytokines.csv");
featureMeta = readmatrix("./Data/20241209_cytokines_featureMeta.csv", FileType="text", OutputType="string");
sampleInfo = readmatrix("./Data/20241209_cytokines_sampleMeta.csv", OutputType="string");

% Prepare y - binary for TM/TW
subjectInfo = unique(sampleInfo(:,[1 2]), "rows");
Y = zeros(size(subjectInfo,1), 1);
Y(subjectInfo(:,2) == "TW",:) = 1;
Y = Y - mean(Y);
Y = [subjectInfo(:,1) Y];

%testosterone = sortrows(testosterone,1);

% Where to save output
path = "./20241212_Cytokines_combined_binary/";

%%
% Log transform
[I,J] = size(df);
vectorizedX = reshape(df, I*J, 1);
pseudocount = min(vectorizedX(vectorizedX>0));

cytokines_log = log(df+pseudocount);

%%
% Make into cube
keepIndividuals = true;
[cytokine_cube, cytokineSubjectMeta, cytokineConditionMeta] = rawDataToCube(cytokines_log, sampleInfo(:,1), str2double(sampleInfo(:,3)), keepIndividuals);

%%
% Center and scale
cytokine_cnt = centerData(cytokine_cube, 1);
cytokine_cnt_scl = scaleData(cytokine_cnt, 2);

%%
% Prepare metadata and remove outliers
cytokine_subjectMeta_filtered = sortrows(unique(sampleInfo(:,[1 2]), "rows"), 1);
cytokine_subjectMeta_filtered([10 19],:) = []; % remove subjects 19 (TM) and 27 (TW) due to being outliers

cytokine_featureMeta_filtered = featureMeta;
cytokine_timeMeta_filtered = cytokineConditionMeta;

cytokine_cnt_scl_filtered = cytokine_cnt_scl;
cytokine_cnt_scl_filtered([10 19], :, :) = [];

Y = str2double(Y(:,2));
Y([10 19],:) = [];

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
[XValResult_cytokine,~] = ncrossreg(cytokine_cnt_scl_filtered, Y, maxComponents, 0);

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
[Xfactors_cytokine,Yfactors_cytokine,Core_cytokine,B_cytokine,ypred_cytokine,ssx_cytokine,ssy_cytokine,reg_cytokine] = npls(cytokine_cnt_scl_filtered, Y, numComponents_cytokines);

%%
% Save NPLS models
metaData = {cytokine_subjectMeta_filtered, cytokine_featureMeta_filtered, cytokine_timeMeta_filtered};
cytokineAnnotatedModel = annotateModel(cytokine_cnt_scl_filtered, Xfactors_cytokine, metaData);
savePARAFAC(cytokine_cnt_scl_filtered, Xfactors_cytokine, cytokineAnnotatedModel, path + "GOHTRANS_cytokine");

%%
% Save ypred of NPLS models
%y_saliva = [ypred_saliva saliva_subjectMeta_filtered testosterone];
y_cytokine = [ypred_cytokine(:,:,1) cytokine_subjectMeta_filtered Y];
writematrix(y_cytokine, path + "Cytokine_" + numComponents_cytokines + "_ypred.csv");

%%
% Save coefficients of NPLS models
for i=1:size(reg_cytokine, 2)
    writematrix(reg_cytokine{i}, path + "Cytokine_" + numComponents_cytokines + "_coeff_" + i + ".csv");
end

%%
% Save crossvalidated coefficients of the NPLS models
[coeff_mean_cytokine, coeff_std_cytokine] = CV_coeff_NPLS(cytokine_cnt_scl_filtered, Y, numComponents_cytokines, 1);
writematrix(coeff_mean_cytokine, path + "Cytokine_" + numComponents_cytokines + "_CVcoeff_means.csv");
writematrix(coeff_mean_cytokine, path + "Cytokine_" + numComponents_cytokines + "_CVcoeff_stds.csv");