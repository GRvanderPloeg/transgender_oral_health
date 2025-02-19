% PARAFAC functionality
addpath("..\N-way-shell\Scripts\"); % own scripts
addpath("..\N-way-shell\N-way toolbox/"); % from Rasmus Bro

%%
% Load data
df = readmatrix("./Data/20240503_UNOISE_new/salivaCounts.csv");
taxonomy = readmatrix("./Data/20240503_UNOISE_new/taxonomySaliva_fixed.csv", FileType="text", OutputType="string");
sampleInfo = readmatrix("./Data/20240503_UNOISE_new/salivaSampleMeta.csv", OutputType="string");

% Remove timepoint 9 months from the data
mask = sampleInfo(:,24) ~= "9";
df = df(mask,:);
sampleInfo = sampleInfo(mask,:);

% Also use only TM samples
salivaMaskTW = sampleInfo(:,7) == "TW";
df = df(salivaMaskTW,:);
sampleInfo = sampleInfo(salivaMaskTW,:);

% Prepare y
% Testosterone is 21
% Estradiol is 17
mask0 = (sampleInfo(:,24) == "0");
mask3 = (sampleInfo(:,24) == "3");
mask6 = (sampleInfo(:,24) == "6");
mask12 = (sampleInfo(:,24) == "12");

testosterone0 = sampleInfo(mask0, [23 21]);
testosterone0 = testosterone0(testosterone0(:,2) ~= "NA",:); % remove NAs
testosterone3 = sampleInfo(mask3, [23 21]);
testosterone3 = testosterone3(testosterone3(:,2) ~= "NA",:); % remove NAs
testosterone6 = sampleInfo(mask6, [23 21]);
testosterone6 = testosterone6(testosterone6(:,2) ~= "NA",:); % remove NAs
testosterone12 = sampleInfo(mask12, [23 21]);
testosterone12 = testosterone12(testosterone12(:,2) ~= "NA",:); % remove NAs

subsetSubjects0 = ismember(sampleInfo(:,23), testosterone0(:,1));
subsetSubjects3 = ismember(sampleInfo(:,23), testosterone3(:,1));
subsetSubjects6 = ismember(sampleInfo(:,23), testosterone6(:,1));
subsetSubjects12 = ismember(sampleInfo(:,23), testosterone12(:,1));

% For 0m, 3m
%usableSubjects = unique(sampleInfo(subsetSubjects0 & subsetSubjects3, 23));
%testosterone = [usableSubjects testosterone0(ismember(testosterone0(:,1), usableSubjects),2) testosterone3(ismember(testosterone3(:,1), usableSubjects),2)];

% For 0m, 12m
usableSubjects = unique(sampleInfo(subsetSubjects0 & subsetSubjects12, 23));
testosterone = [usableSubjects testosterone0(ismember(testosterone0(:,1), usableSubjects),2) testosterone12(ismember(testosterone12(:,1), usableSubjects),2)];

% Finally
testosterone = sortrows(testosterone,1);

% Where to save output
path = "./temp/";

%%
% Filter based on sparsity per group per niche
sparsity_saliva = sum(df==0, 1) / size(df,1);

% subplot(2,1,1); histogram(sparsity_tongue_TM);
% subplot(2,1,2); histogram(sparsity_tongue_TW);

% Saliva
threshold = 0.5;
featureSelection = sparsity_saliva <= threshold;
%%
% CLR
saliva_clr = transformCLR(df);

%% Filtering after CLR
saliva_filtered = saliva_clr(:,featureSelection);
taxa_saliva_filtered = taxonomy(featureSelection,:);

%%
% Make into cube
keepIndividuals = true;
[saliva_cube, salivaSubjectMeta, salivaConditionMeta] = rawDataToCube(saliva_filtered, sampleInfo(:,23), str2double(sampleInfo(:,24)), keepIndividuals);

%%
% Remove subjects with missing testosterone measurements
mask = ismember(salivaSubjectMeta, usableSubjects);
saliva_cube = saliva_cube(mask,:,:);
salivaSubjectMeta = salivaSubjectMeta(mask,:);

%%
% Center and scale
saliva_cnt = centerData(saliva_cube, 1);
saliva_cnt_scl = scaleData(saliva_cnt, 2);

%%
% Prepare metadata and remove outliers
saliva_subjectMeta_filtered = sortrows(unique(sampleInfo(:,[23 7]), "rows"), 1);
saliva_subjectMeta_filtered = saliva_subjectMeta_filtered(mask,:);
saliva_featureMeta_filtered = taxa_saliva_filtered;
saliva_timeMeta_filtered = salivaConditionMeta;

saliva_cnt_scl_filtered = saliva_cnt_scl;

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
[XValResult_saliva,~] = ncrossreg(saliva_cnt_scl_filtered, testosterone, maxComponents, 0);

%%
% Plot RMSEP of CV and save it
%set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
%plot(1:maxComponents, XValResult.RMSEP); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP"); yyaxis right; plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); ylim([0 100]); ylabel("Variance explained (%)");
%saveas(gcf, path+"RMSEP_CV.jpg");
%close();
XValResult = XValResult_saliva;

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
numComponents_saliva = 2;
bestRMSEP_saliva = XValResult_saliva.RMSEP(numComponents_saliva);

%%
% Run NPLS
[Xfactors_saliva,Yfactors_saliva,Core_saliva,B_saliva,ypred_saliva,ssx_saliva,ssy_saliva,reg_saliva] = npls(saliva_cnt_scl_filtered, testosterone, numComponents_saliva);

%%
% Save NPLS models
metaData = {saliva_subjectMeta_filtered, saliva_featureMeta_filtered, saliva_timeMeta_filtered};
salivaAnnotatedModel = annotateModel(saliva_cnt_scl_filtered, Xfactors_saliva, metaData);
savePARAFAC(saliva_cnt_scl_filtered, Xfactors_saliva, salivaAnnotatedModel, path + "GOHTRANS_saliva");

%%
% Save ypred of NPLS models
% Save ypred of NPLS models
T = str2double(salivaAnnotatedModel{1}(:,1:numComponents_saliva));
b = inv(T'*T) * T' * testosterone;
ypred_saliva_combined = T * b;
y_saliva = [ypred_saliva_combined saliva_subjectMeta_filtered testosterone];
writematrix(y_saliva, path + "Saliva_" + numComponents_saliva + "_ypred.csv");


%y_saliva = [ypred_saliva saliva_subjectMeta_filtered testosterone];
%y_saliva = [ypred_saliva(:,:,1) saliva_subjectMeta_filtered testosterone];
%writematrix(y_saliva, path + "Saliva_" + numComponents_saliva + "_ypred.csv");

%%
% Save coefficients of NPLS models
for i=1:size(reg_saliva, 2)
    writematrix(reg_saliva{i}, path + "Saliva_" + numComponents_saliva + "_coeff_" + i + ".csv");
end

%%
% Save crossvalidated coefficients of the NPLS models
[coeff_mean_saliva, coeff_std_saliva] = CV_coeff_NPLS(saliva_cnt_scl_filtered, testosterone, numComponents_saliva, 1);
writematrix(coeff_mean_saliva, path + "Saliva_" + numComponents_saliva + "_CVcoeff_means.csv");
writematrix(coeff_mean_saliva, path + "Saliva_" + numComponents_saliva + "_CVcoeff_stds.csv");