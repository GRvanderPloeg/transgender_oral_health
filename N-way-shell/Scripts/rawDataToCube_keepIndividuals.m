function Xcube=rawDataToCube(Xlong, Xmeta_subjects, Xmeta_thirdMode, keepIndividuals)
%
% rawDataToCube_keepIndividuals: converts a (long) matricized array with
% every row corresponding to a sample to a data cube for PARAFAC modelling.
% This version of rawDataToCube conserves missing samples by adding NaNs.
% 
% See also:
% 'parafac', 'rawDataToCube'
%
% Xcube=rawDataToCube_keepIndividuals(Xlong, Xmeta_subjects, Xmeta_thirdMode, convertThirdModeToNumeric)
%
% Creates a data cube out of a matricized array using supplied sample and
% third mode metadata. Missing sample-timepoint combinations are added as
% NaN vectors.
%
% NOTE: this sorts the subjects alphabetically. Make sure to adapt your
% metadata accordingly!
%
% ---------------------- INPUT ---------------------
%
% Xlong             This is the (long) matricized input data array. This
%                   data is expected to be completely numeric with the
%                   metadata removed.
%
% Xmeta_subjects    Metadata of the subjects. Use only the sampleIDs (i.e.
%                   the names of your subjects). Expected to be strings.
%
% Xmeta_thirdMode   Metadata of the third mode. In most cases this will be
%                   a vector of timepoints. Can be string or numeric.
%
% keepIndividuals   Boolean specifying if missing samples should be
%                   conserved as NaN values.
%
% ---------------------- OUTPUT ---------------------
%
% Xcube             The dataset in cube format. The subjects will be in the
%                   first mode, the features in the second and the time in
%                   the third. Note that missing samples are added as NaN.

% Created by: G.R. van der Ploeg (g.r.ploeg@uva.nl)

% Create metadata overview
Xmeta = [Xmeta_subjects Xmeta_thirdMode];
uniqueThirdMode = unique(Xmeta_thirdMode);
numUniqueThirdMode = size(uniqueThirdMode,1);

% Identify how many samples are missing
subjects = unique(Xmeta_subjects);
expectedNumSamples = length(subjects) * numUniqueThirdMode;
realNumSamples = size(Xlong,1);

numMissingSamples = expectedNumSamples - realNumSamples;

% Identify which samples are missing
missingSamples = zeros(numMissingSamples, 2);
missingSamples = string(missingSamples);
numOccurences = groupcounts(Xmeta(:,1));
missingSampleIterator = 1;

% Identification of individuals missing one or more samples
for i=1:length(subjects)
    subject = subjects(i);
    discrepancy = numUniqueThirdMode - numOccurences(i);
    if (discrepancy > 0)
        for j=1:discrepancy
            missingSamples(missingSampleIterator, 1) = subject;
            missingSampleIterator = missingSampleIterator + 1;
        end
    end 
end

% Identification of the timepoints that are missing
missingSubjects = unique(missingSamples(:,1));
missingSampleIterator = 1;

for i=1:length(missingSubjects)
    subject = missingSubjects(i,1);
    timepoints = str2double(Xmeta(Xmeta(:,1) == subject,2));
    expectedTimepoints = uniqueTimepoints;

    ii = ~ismember(expectedTimepoints,timepoints);
    values = expectedTimepoints(ii);

    if (length(values) > 1)
        for j=1:length(values)
            missingSamples(missingSampleIterator,2) = string(values(j));
            missingSampleIterator = missingSampleIterator + 1;
        end
    else
        missingSamples(missingSampleIterator,2) = string(values);
        missingSampleIterator = missingSampleIterator + 1;
    end
end

% Append them to the dataset
missingData = nan(numMissingSamples, size(Xlong,2));
Xlong = [Xlong; missingData];

% Sort the data
Xmeta = [Xmeta; missingSamples];
Xmeta(:,end+1) = 1:size(Xmeta,1);

subjectIDs = Xmeta(:,1);
timepoints = str2double(Xmeta(:,2));
index = str2double(Xmeta(:,3));
T = table(subjectIDs, timepoints, index);
T_sorted = sortrows(T, [2 1]);

%Xmeta = sortrows(Xmeta, [2 1]); % sort on visit number (ascending), then subject name (alphabetical)
%keepRowIndices = str2double(Xmeta(:,end)); % new sorting of data
keepRowIndices = T_sorted.index;
Xlong = Xlong(keepRowIndices, :);

% Reshape into cube
I = size(Xlong, 1) / numTimepoints;
J = size(Xlong, 2);
K = numTimepoints;

Xcube = reshape(Xlong, I, K, J);
Xcube = permute(Xcube, [1 3 2]);

