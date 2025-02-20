function plotTuckerCongruences(allTuckers, numComponents, overallTitle, path)

numModes = size(allTuckers, 2);

if ~exist('path','var')
    path="";
end

numModels = size(allTuckers{1},1);
numComparisons = size(allTuckers{1},2);

% Find the original number of components used by solving an equation
% nFactors1 = (1 + sqrt(1+4*1*2*nComparisons)) / 2; % Quadratic equation formula
% nFactors2 = (1 - sqrt(1+4*1*2*nComparisons)) / 2; % Quadratic equation formula
% 
% if nFactors1 > 0
%     nComponents = nFactors1;
% elseif nFactors2 > 0
%     nComponents = nFactors2;
% end

% Generate titles for all tucker comparisons made
comparisonTitles = strings(1,numComparisons);
comparisonIterator = 1;
for f1=1:(numComponents-1)
    for f2=(f1+1):numComponents
        comparisonTitles(comparisonIterator) = "Component " + f1 + " vs. " + f2;
        comparisonIterator = comparisonIterator + 1;
    end
end

% Make the plots
plotIterator = 1;
set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

for i=1:numComparisons

    % Plot TCC for every mode
    overallTCC = ones(size(allTuckers{1}(:,i)));
    for j=1:numModes
        subplot(numComparisons,numModes+1,plotIterator); hold on;
        bar(1:numModels, allTuckers{j}(:,i));
        title("Tuckers congruences, mode " + j + " (" + comparisonTitles(i) + ")");
        xlabel("Model number");
        ylabel("Tucker congruence coefficient");
        plotIterator = plotIterator + 1;
        overallTCC = overallTCC .* allTuckers{j}(:,i);
    end

    % Plot overall TCC by multiplying the TCCs
    subplot(numComparisons,numModes+1,plotIterator); hold on;
    bar(1:numModels, overallTCC);
    title("Overall tucker congruence (" + comparisonTitles(i) + ")");
    xlabel("Model number");
    ylabel("Tucker congruence coefficient");
    plotIterator = plotIterator + 1;
end

sgtitle(overallTitle);

if path ~= ""
    saveas(gcf, path)
    close()
end