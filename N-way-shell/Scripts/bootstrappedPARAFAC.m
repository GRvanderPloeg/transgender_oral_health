function [allModels,allVarExps]=bootstrappedPARAFAC(X, subjectMeta, subjectGroupCol, numFactors, maxIterations, const, balanced)
%
% bootstrappedPARAFAC: creates jack-knifed PARAFAC models to determine
% loading stability.
% 
% See also:
% 'parafac'
%
% [allModels,allVarExps]=bootstrappedPARAFAC(X, subjectMeta, subjectGroupCol, numFactors, maxIterations, const, balancedJackKnifing)
%
% Creates many jack-knifed versions of a PARAFAC model while storing all
% model output for later inspection.
%
% Has two modes:
% 1. Jack-knife in a balanced way across subject groups.
% 2. Jack-knife random samples.
%
% ---------------------- INPUT ---------------------
%
% X                 This is the input data cube. This data is expected
%                   to be completely numeric with the metadata removed.
%
% subjectMeta       Metadata of the subjects.
%
% subjectGroupCol   Column number index to be used for subject groups. Set
%                   to 0 if not used.
%
% numFactors        Number of factors of the PARAFAC model.
%
% maxIterations     The maximum number of jack-knifed models that must be
%                   created. Note that this will throw an error if it exceeds
%                   the total number of samples or sample combinations.
%
% const             PARAFAC model constraints (see parafac.m documentation)
%
% balanced          Boolean setting the jack-knifing mode:
%                   True: jack-knife in a balanced way acrosss subject
%                   groups.
%                   False: jack-knife random samples.
%
% ---------------------- OUTPUT ---------------------
%
% allModels         Cell object containing the loadings of every generated
%                   model.
%
% allVarExps        List object containing the variance explained of every
%                   generated model.

% Created by: G.R. van der Ploeg (g.r.ploeg@uva.nl)

Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 1;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

numModes = size(size(X),2);
I = size(X,1);
allModels = {};
for i=1:numModes
    allModels{i} = [];
end

if balanced == true && subjectGroupCol == 0
    % Find subject groups
    subjectGroups = unique(subjectMeta(:,subjectGroupCol));
    numSubjectGroups = size(subjectGroups, 1);
    groupMembership = {};
    
    for i=1:numSubjectGroups
        groupMembership{i} = find(subjectMeta(:,subjectGroupCol) == subjectGroups(i));
    end

    % Establish which individuals to remove randomly
    combinations = nchoosek(1:size(X,1), numSubjectGroups);
    combinationGroups = {};
    
    for i=1:numSubjectGroups
        combinationGroups{i} = ismember(combinations, groupMembership{i});
    end
    
    % Check for balanced jack-knifing - one subject from each group must be
    % removed.
    valid = ones(size(combinations,1), 1);
    
    for i=1:numSubjectGroups
        valid = valid & (sum(combinationGroups{i},2) == 1);
    end
    
    combinations = combinations(valid,:);
    
else
    combinations = randi([1 size(X,1)],1,maxIterations)';
end

jackKnifedSamples = 1:maxIterations;
allVarExps = zeros(size(jackKnifedSamples,2),1);
textprogressbar('Creating bootstrapped PARAFAC models ');
for n=1:size(jackKnifedSamples,2)
    textprogressbar((n/size(jackKnifedSamples,2))*100);
    Xboot = X;
    %IDremoval = [RFlow(RFlow_remove(n)) RFmid(RFmid_remove(n)) RFhigh(RFhigh_remove(n))];
    selection = jackKnifedSamples(n);
    IDremoval = combinations(selection,:);

    % Remove row without hardcoding dimensions of the array
    S = [];
    S.subs = repmat({':'},1,ndims(Xboot));
    S.subs{1} = IDremoval;
    S.type = '()';
    Xboot = subsasgn(Xboot,S,[]);
    
    % Run PARAFAC algorithm
    evalc('[Factors, it, err, c] = silent_parafac(Xboot, numFactors, Options, const)');
    sortedModel = sortParafacComponents(Xboot, Factors);

    % Make ID loading vector with NaNs for the removed individuals
    Amodified = zeros(size(X,1),numFactors);
    IDiterator = 1;
    for i=1:I
        if(ismember(i,IDremoval))
            Amodified(i,:) = nan(1,numFactors);
        else
            Amodified(i,:) = sortedModel{1}(IDiterator,:);
            IDiterator = IDiterator + 1;
        end
    end

    % Store results
    for i=1:numModes
        if i==1
            allModels{i}(:,:,n) = Amodified;
        else
            allModels{i}(:,:,n) = sortedModel{i};
        end
    end

    allVarExps(n) = calcVarExplained(Xboot, sortedModel);
end

textprogressbar(' done');