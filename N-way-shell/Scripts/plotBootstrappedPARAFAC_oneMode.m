function plotBootstrappedPARAFAC_oneMode(models, mode, fixNegativity, varExps, overallTitle, path)

numComponents = size(models{1}, 2);

if fixNegativity == true 
    negativityMatrix = checkComponentNegativity(models);
    models = fixComponentNegativity(models, negativityMatrix);
end

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
plotIterator = 1;

for i=1:numComponents
    [I,~,~] = size(models{mode});
    df = reshape(models{mode}(:,i,:), I, []);
    plottedData = df';

    if isa(plottedData, "string") % convert to numeric if the annotation is a string
        plottedData = str2double(plottedData);
    end

    % Box plot
    subplot(numComponents, 1, plotIterator), hold on;
    boxplot(plottedData);

    ylabel("Comp. " + i);
    xlabel("Index");

    hold off;
    plotIterator = plotIterator + 1;
end

sgtitle(overallTitle + ", variance explained = " + round(mean(varExps),2) + " +/- " + round(std(varExps),2) + " %");
saveas(gcf, path)
close()
end