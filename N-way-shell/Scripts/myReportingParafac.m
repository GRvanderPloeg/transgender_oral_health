function [allModels,allIterations,allErrors,allCorcondias,allVarExp,allTuckers, allFPIs, allEPIs]=myReportingParafac(X, numFactors, numReps, Options, const)

numModes = size(size(X),2);
allModels = cell(1,numModes);
allIterations = 1:numReps;
allErrors = 1:numReps;
allCorcondias = 1:numReps;
allVarExp = 1:numReps;
allTuckers = cell(1,numModes);
allFPIs = {};
allEPIs = {};

for i=1:numModes
    allModels{i} = [];
end

textprogressbar('Creating PARAFAC models ');
for i=1:numReps
    textprogressbar(i/numReps*100);
    evalc('[Factors, it, err, c, FPI, EPI] = silent_parafac(X, numFactors, Options, const)');

    for j=1:numModes
        allModels{j}(:,:,i) = Factors{j};
    end
    allIterations(i) = it;
    allErrors(i) = err;
    allCorcondias(i) = c;
    allVarExp(i) = calcVarExplained(X, Factors);
    allEPIs{i} = EPI;
    Tuckers = calcTuckerCongruence(Factors);
    for j=1:numModes
        allTuckers{j}(i,:) = Tuckers{j};
        allFPIs{j,i} = FPI{j};
    end
end

textprogressbar(' done');