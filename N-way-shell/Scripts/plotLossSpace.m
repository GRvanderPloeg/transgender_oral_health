function plotLossSpace(FactorsPerIteration, ErrorPerIteration, modeToPlot, path_start)

MAXPLOTSPERIMAGE = 25;

colours =           [255 0 0;
                    0 255 0;
                    0 0 255;
                    0 0 0;
                    0 255 255;
                    255 0 255;
                    128 128 128;
                    128 0 0;
                    128 128 0;
                    0 128 0;
                    128 0 128;
                    0 128 128;
                    0 0 128;
                    191 191 191;
                    0 0 0]/255;

numModes = size(FactorsPerIteration, 1);
numReps = size(FactorsPerIteration, 2);
numComponents = size(FactorsPerIteration{1}, 2);

numItemsToPlot = size(FactorsPerIteration{modeToPlot,1}, 1);
if numItemsToPlot > MAXPLOTSPERIMAGE
    numRows = ceil(sqrt(MAXPLOTSPERIMAGE));
    numCols = numRows;
    splitUp = true;
else
    numRows = ceil(sqrt(numItemsToPlot));
    numCols = numRows;
    splitUp = false;
end

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
plotIterator = 1;
partCounter = 1;

for j=1:numItemsToPlot
    subplot(numRows,numCols,plotIterator); hold on;

    for k=1:numReps
        numIterations = size(ErrorPerIteration{k},2);
        for l=1:numComponents
            df = reshape(FactorsPerIteration{modeToPlot,k}(j,l,:), [], numIterations);
            plot(1:numIterations, df, "Color",colours(l,:));
        end
    end

    xlabel("Iteration");
    ylabel("Item " + j + " loading");
    hold off;
    plotIterator = plotIterator + 1;

    % When max plots per image is reached, start over.
    if plotIterator > MAXPLOTSPERIMAGE
        saveas(gcf, path_start+"_ncomp"+numComponents+"_mode"+modeToPlot+"_part"+partCounter+".jpg")
        close()
        set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
        plotIterator = 1;
        partCounter = partCounter + 1;
    end
end

if splitUp == true
    saveas(gcf, path_start+"_ncomp"+numComponents+"_mode"+modeToPlot+"_part"+partCounter+".jpg")
else
    saveas(gcf, path_start+"_ncomp"+numComponents+"_mode"+modeToPlot+".jpg")
end
close()

end
