% PARAFAC functionality
addpath("..\N-way-shell\Scripts\"); % own scripts
addpath("..\N-way-shell\N-way toolbox/"); % from Rasmus Bro

%%
% Load data
df = readmatrix("./Data/20240503_UNOISE_new/tongueCounts.csv");
taxonomy = readmatrix("./Data/20240503_UNOISE_new/taxonomyTongue_fixed.csv", FileType="text", OutputType="string");
sampleInfo = readmatrix("./Data/20240503_UNOISE_new/tongueSampleMeta.csv", OutputType="string");

% Remove timepoint 9 months from the data
% Yes, this is REQUIRED
mask = sampleInfo(:,24) ~= "9";
df = df(mask,:);
sampleInfo = sampleInfo(mask,:);

% Separate TM and TW samples
tongueMaskTM = sampleInfo(:,7) == "TM";
tongueMaskTW = sampleInfo(:,7) == "TW";
tongueTM = df(tongueMaskTM,:);
tongueTW = df(tongueMaskTW,:);

% Prepare y - binary for TM/TW
subjectInfo = unique(sampleInfo(:,[23, 7]), "rows");
Y = zeros(size(subjectInfo,1), 1);
Y(subjectInfo(:,2) == "TW",:) = 1;
Y = Y - mean(Y);
Y = [subjectInfo(:,1) Y];

% Finally
%testosterone = sortrows(testosterone,1);

% Where to save output
path = "./20241212_Tongue_UNOISE_combined_binary/";

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

%%
% CLR
tongue_clr = transformCLR(df);

%% Filter after CLR
tongue_filtered = tongue_clr(:,featureSelection);
taxa_tongue_filtered = taxonomy(featureSelection,:);

%%
% Make into cube
keepIndividuals = true;
[tongue_cube, tongueSubjectMeta, tongueConditionMeta] = rawDataToCube(tongue_filtered, sampleInfo(:,23), str2double(sampleInfo(:,24)), keepIndividuals);

%%
% Center and scale
tongue_cnt = centerData(tongue_cube, 1);
tongue_cnt_scl = scaleData(tongue_cnt, 2);

%%
% Prepare metadata and remove outliers
tongue_subjectMeta_filtered = sortrows(unique(sampleInfo(:,[23 7]), "rows"), 1);
tongue_featureMeta_filtered = taxa_tongue_filtered;
tongue_timeMeta_filtered = tongueConditionMeta;

tongue_cnt_scl_filtered = tongue_cnt_scl;

Y = str2double(Y(:,2));

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
[XValResult_tongue,~] = ncrossreg(tongue_cnt_scl_filtered, Y, maxComponents, 0);

%%
% Plot RMSEP of CV and save it
%set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
%plot(1:maxComponents, XValResult.RMSEP); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP"); yyaxis right; plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); ylim([0 100]); ylabel("Variance explained (%)");
%saveas(gcf, path+"RMSEP_CV.jpg");
%close();
XValResult = XValResult_tongue;

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
%[bestRMSEP_tongue, numComponents_tongue] = min(XValResult_tongue.RMSEP);

% OVERRIDE
numComponents_tongue = 1;
bestRMSEP_tongue = XValResult_tongue.RMSEP(numComponents_tongue);

%%
% Run NPLS
[Xfactors_tongue,Yfactors_tongue,Core_tongue,B_tongue,ypred_tongue,ssx_tongue,ssy_tongue,reg_tongue] = npls(tongue_cnt_scl_filtered, Y, numComponents_tongue);

%%
% Save NPLS models
metaData = {tongue_subjectMeta_filtered, tongue_featureMeta_filtered, tongue_timeMeta_filtered};
tongueAnnotatedModel = annotateModel(tongue_cnt_scl_filtered, Xfactors_tongue, metaData);
savePARAFAC(tongue_cnt_scl_filtered, Xfactors_tongue, tongueAnnotatedModel, path + "GOHTRANS_tongue");

%%
% Save ypred of NPLS models
y_tongue = [ypred_tongue tongue_subjectMeta_filtered Y];
writematrix(y_tongue, path + "Tongue_" + numComponents_tongue + "_ypred.csv");

%%
% Save coefficients of NPLS models
for i=1:size(reg_tongue, 2)
    writematrix(reg_tongue{i}, path + "Tongue_" + numComponents_tongue + "_coeff_" + i + ".csv");
end

%%
% Save crossvalidated coefficients of the NPLS models
[coeff_mean_tongue, coeff_std_tongue] = CV_coeff_NPLS(tongue_cnt_scl_filtered, Y, numComponents_tongue, 1);
writematrix(coeff_mean_tongue, path + "Tongue_" + numComponents_tongue + "_CVcoeff_means.csv");
writematrix(coeff_mean_tongue, path + "Tongue_" + numComponents_tongue + "_CVcoeff_stds.csv");