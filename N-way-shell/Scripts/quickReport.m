function [allModels,allCons,allVarExps, bootstrappedModels, bootVarExps, allTuckers]=quickReport(X, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, overallTitle, path_start)
% [allModels,allCons,allVarExps, bootstrappedModels, bootVarExps, allTuckers]=quickReport(X, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, overallTitle, path_start)

numModes = size(size(X),2);
sizeVector = size(X);
allModels = cell(numModes,maxComponents);
allCons = cell(1,maxComponents);
allVarExps = cell(1,maxComponents);
bootstrappedModels = cell(numModes, maxComponents);
bootVarExps = cell(1,maxComponents);
allTuckers = cell(1, maxComponents);

for n=1:maxComponents
    [models, Its, Errs, Cons, varExps, Tuckers, allFPIs, allEPIs] = myReportingParafac(X, n, numReps, Options, const);
    plotReporterOutput(Its, Errs, Cons, varExps, overallTitle, path_start+"_model_report_comp"+n+".jpg");

    % Report the number of models that have one or more flipped component
    negativityMatrix = checkComponentNegativity(models);
    [~, numFixed] = fixComponentNegativity(models, negativityMatrix);
    disp(numFixed + " out of " + numReps + " randomly initialized " + n + "-component models had one or more flipped components.")

    % for i=1:numModes
    %     plotLossSpace(allFPIs, allEPIs, i, path_start + "_loss");
    % end

    if n >= 2
        plotTuckerCongruences(Tuckers, n, overallTitle, path_start+"_model_TCC_report_comp"+n+".jpg");
        allTuckers{n} = Tuckers;
    end

    for i=1:numModes
        allModels{i,n} = models{i};
    end
    allCons{n} = Cons;
    allVarExps{n} = varExps;

    [boots,varExps] = bootstrappedPARAFAC(X, metaData{1}, subjectGroupCol, n, maxIterations, const, balancedJackKnifing);

    for i=1:numModes
        bootstrappedModels{i,n} = boots{i};
    end
    bootVarExps{n} = varExps;

    % Report the number of models that have one or more flipped component
    negativityMatrix = checkComponentNegativity(boots);
    [~, numFixedBoots] = fixComponentNegativity(boots, negativityMatrix);
    disp(numFixedBoots + " out of " + maxIterations + " jack-knifed " + n + "-component models had one or more flipped components.")

    plotBootstrappedPARAFAC(boots, bootVarExps{n}, metaData, resort, false, legendIndex, numColsPerLegend, xlabels, titles, overallTitle, path_start+"_bootstrapped_model_comp"+n+".jpg");
    plotBootstrappedPARAFAC(boots, bootVarExps{n}, metaData, resort, true, legendIndex, numColsPerLegend, xlabels, titles, overallTitle, path_start+"_bootstrapped_model_comp"+n+"_fixed.jpg");

%     for i=1:numModes
%         plotBootstrappedPARAFAC_oneMode(boots, i, false, bootVarExps{n}, overallTitle, path_start+"_bootstrapped_model_comp"+n+"_mode"+i+".jpg");
%         plotBootstrappedPARAFAC_oneMode(boots, i, true, bootVarExps{n}, overallTitle, path_start+"_bootstrapped_model_comp"+n+"_mode"+i+"_fixed.jpg");
%     end
end