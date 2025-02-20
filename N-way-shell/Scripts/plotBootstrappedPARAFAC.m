function plotBootstrappedPARAFAC(models, varExps, metaData, resort, fixNegativity, legendIndex, numColsPerLegend, xlabels, titles, overallTitle, path)

colours =    [122 62 118; % rgb(122, 62, 118)
                    61 105 123; % rgb(61, 105, 123)
                    49 155 49; % rgb(49, 155, 49)
                    183 107 1; % rgb(183, 107, 1)
                    184 3 0; % rgb(184, 3, 0)
                    180 106 175; % rgb(180, 106, 175)
                    91 151 174; % rgb(91, 151, 174)
                    84 201 84; % rgb(84, 201, 84)
                    254 153 11; % rgb(254, 153, 11)
                    255 34 31; % rgb(255, 34, 31)
                    204 153 201; % rgb(204, 153, 201)
                    158 193 207; % rgb(158, 193, 207)
                    158 224 158; % rgb(158, 224, 158)
                    254 177 68; % rgb(254, 177, 68)
                    255 102 99; % rgb(255, 102, 99)
                    255 0 0;
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

numModes = size(models,2);
numComponents = size(models{1}, 2);

if fixNegativity == true 
    negativityMatrix = checkComponentNegativity(models);
    models = fixComponentNegativity(models, negativityMatrix);
end

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
plotIterator = 1;

for i=1:numComponents
    for j=1:numModes
        [I,~,~] = size(models{j});
        df = reshape(models{j}(:,i,:), I, []);
        plottedData = [median(df, 2, "omitnan") quantile(df, 0.25, 2) quantile(df, 0.75, 2)];

        if legendIndex(j) ~= 0
            names = metaData{j}(:,legendIndex(j));
            grp = grp2idx(metaData{j}(:,legendIndex(j)));
            index = 1:length(grp);
            plottedData = [plottedData grp names index'];

            if resort(j) == true
                plottedData = sortrows(plottedData, [4 6]);
            end
            grp = plottedData(:,4);
            if isa(grp, "string")
                grp = str2double(grp);
            end
            names = unique(plottedData(:,5), "stable");
            plottedData = plottedData(:,[1 2 3]);
        end

        if isa(plottedData, "string") % convert to numeric if the annotation is a string
            plottedData = str2double(plottedData);
        end

        % Bar plot
        subplot(numComponents, numModes, plotIterator), hold on;
        plot = bar(plottedData(:,1), "FaceColor", "flat");

        % Add error bars
        er = errorbar(1:size(plottedData,1), plottedData(:,1), plottedData(:,1)-plottedData(:,2), plottedData(:,3)-plottedData(:,1));
        er.Color = [0 0 0];
        er.CapSize = 0;
        er.LineStyle = 'none';

        % Add axis text and subtitles
        if i==1
            title(titles(j));
        end
        if rem(plotIterator, numModes) == 1 % Only add y label to leftmost plot per row
            ylabel("Comp. " + i);
        end
        xlabel(xlabels(j));

        if legendIndex(j) ~= 0

            % Colouring of the bars
            for k=1:size(grp, 1)
                plot.CData(k,:) = colours(grp(k),:);
            end

            % Addition of legend
            h = zeros(max(grp), 1);
            for k=1:max(grp)
                h(k) = patch(NaN, NaN, colours(k,:));
            end
            legend(h, string(names), "Location", "southoutside", "NumColumns", numColsPerLegend(j));
        end

        hold off;
        plotIterator = plotIterator + 1;
    end
end

sgtitle(overallTitle + ", variance explained = " + round(mean(varExps),2) + " +/- " + round(std(varExps),2) + " %");
saveas(gcf, path)
close()
end