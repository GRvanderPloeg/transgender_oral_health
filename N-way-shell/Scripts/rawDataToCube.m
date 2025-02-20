function [Xcube,newSubjectMeta,newConditionMeta]=rawDataToCube(Xlong, subjectMeta, conditionMeta, keepIndividuals)
%
% rawDataToCube: converts a (long) matricized array with every row 
% corresponding to a sample to a data cube for PARAFAC modelling.
% Can add missing samples as NaN values using the keepIndividuals
% parameter.
% 
% See also:
% 'parafac'
%
% [Xcube,newSubjectMeta,newConditionMeta]=rawDataToCube(Xlong, subjectMeta, conditionMeta, keepIndividuals)
%
% Creates a data cube out of a matricized array using supplied sample and
% third mode metadata. Missing samples can be added as NaN vectors.
%
% Expects matricization to be in the 'long' direction, i.e. the rows are
% repeated sets of samples 1 ... I, with each set corresponding to a third
% mode condition 1 ... K.
%
% If 'conditionMeta' contains multiple copies of the same condition for the
% same subject, then only the first occurence will be used for the cube.
%
% ---------------------- INPUT ---------------------
%
% Xlong             This is the (long) matricized input data array. This
%                   data is expected to be completely numeric with the
%                   metadata removed.
%
% subjectMeta       Metadata of the subjects. Use only the sampleIDs (i.e.
%                   the names of your subjects). Expected to be strings.
%
% conditionMeta     Metadata of the third mode. In most cases this will be
%                   a vector of timepoints. Can be string or numeric.
%
% keepIndividuals   Boolean specifying if missing samples should be
%                   conserved as NaN values.
%
% ---------------------- OUTPUT ---------------------
%
% Xcube             The dataset in cube format. The subjects will be in the
%                   first mode and the features in the second.
%
% newSubjectMeta    Gives back the way in which the subjects were ordered
%                   corresponding to the first mode of the cube.
%
% newConditionMeta  Gives back the way in which the third mode was ordered.

% Created by: G.R. van der Ploeg (g.r.ploeg@uva.nl)

% Check for duplicates in metadata
meta = [subjectMeta conditionMeta];
[u,I,~] = unique(meta, "rows", "first");

if size(u) < size(meta,1)
    disp('Duplicate sample metadata detected.');
    disp('Cube output will only use the first occurence of each unique subject-condition combination.');
    disp('This behavior can be avoided by making sure every subject-condition combination only occurs once.')
    disp('The following number of samples will be removed:')
    removed = size(Xlong,1)-length(I);
    disp(removed)
    Xlong = Xlong(I,:);
    subjectMeta = subjectMeta(I,:);
    conditionMeta = conditionMeta(I,:);
end

% Determine what the third mode looks like
uniqueConditions = unique(conditionMeta);
numConditions = size(uniqueConditions,1);

% Identify how many samples are missing
subjects = unique(subjectMeta);
expectedNumSamples = length(subjects) * numConditions;
realNumSamples = size(Xlong,1);

numMissingSamples = expectedNumSamples - realNumSamples;

% Identification of individuals missing one or more samples
if isa(subjectMeta, "string")
    missingSubjects = strings(numMissingSamples, 1);
else
    missingSubjects = zeros(numMissingSamples, 1);
end

numOccurences = groupcounts(subjectMeta);
missingSampleIterator = 1;

for i=1:length(subjects)
    subject = subjects(i);
    discrepancy = numConditions - numOccurences(i);
    if (discrepancy > 0)
        for j=1:discrepancy
            missingSubjects(missingSampleIterator) = subject;
            missingSampleIterator = missingSampleIterator + 1;
        end
    end 
end

if keepIndividuals == true    

    % Identification of the third mode conditions that are missing
    if isa(conditionMeta, "string")
        missingConditions = strings(numMissingSamples, 1);
    else
        missingConditions = zeros(numMissingSamples, 1);
    end
    
    uniqueMissingSubjects = unique(missingSubjects);
    missingSampleIterator = 1;
    
    for i=1:length(uniqueMissingSubjects)
        subject = uniqueMissingSubjects(i);
        conditions = conditionMeta(subjectMeta == subject);
        expectedConditions = uniqueConditions;
    
        ii = ~ismember(expectedConditions,conditions);
        values = expectedConditions(ii);
    
        if (length(values) > 1)
            for j=1:length(values)
                missingConditions(missingSampleIterator) = string(values(j));
                missingSampleIterator = missingSampleIterator + 1;
            end
        else
            missingConditions(missingSampleIterator) = string(values);
            missingSampleIterator = missingSampleIterator + 1;
        end
    end
    
    % Append missing data to the dataset
    missingData = nan(numMissingSamples, size(Xlong,2));
    Xlong = [Xlong; missingData];

    % Create metadata table
    Xmeta = table([subjectMeta; missingSubjects], [conditionMeta; missingConditions], 'VariableNames', {'subjectID', 'condition'});

else

    % Remove subjects that have missing samples from the data
    removeMask = ismember(subjectMeta, missingSubjects);
    subjectMeta(removeMask) = [];
    conditionMeta(removeMask) = [];
    Xlong(removeMask,:) = [];
    
    % Create metadata table
    Xmeta = table(subjectMeta, conditionMeta, 'VariableNames', {'subjectID', 'condition'});
    
end

% Sort the data by condition, then subjectID
Xmeta.index = (1:size(Xmeta,1))';
Xmeta_sorted = sortrows(Xmeta, [2 1]);
keepRowIndices = Xmeta_sorted.index;
Xlong = Xlong(keepRowIndices, :);

% Manually remove two cases of 5 replicates

% Reshape into cube
I = size(Xlong, 1) / numConditions;
J = size(Xlong, 2);
K = numConditions;

Xcube = reshape(Xlong, I, K, J);
Xcube = permute(Xcube, [1 3 2]);

% Supply metadata output
newSubjectMeta = unique(Xmeta_sorted.subjectID, "stable");
newConditionMeta = unique(Xmeta_sorted.condition, "stable");
