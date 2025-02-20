function plotPARAFAC4(annotatedModel, numFactors, varExps, choice, resort, legendIndex, numColsPerLegend, xlabels, titles, overallTitle, path)

colours =           [122 62 118; % rgb(122, 62, 118)
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

numModes = size(annotatedModel, 2);
set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
plotIterator = 1;

for i=1:numFactors
    for j=1:numModes
        plottedData = annotatedModel{j}(:,i);

        if legendIndex(j) ~= 0
            names = annotatedModel{j}(:,legendIndex(j));
            grp = grp2idx(annotatedModel{j}(:,legendIndex(j)));
            index = 1:length(grp);
            plottedData = [plottedData grp names index'];

            if resort(j) == true
                plottedData = sortrows(plottedData, [2 4]); % sort on group number and then index
            end

            grp = plottedData(:,2);
            if isa(grp, "string")
                grp = str2double(grp);
            end
            names = unique(plottedData(:,3), "stable");
            plottedData = plottedData(:,1);
        end

        if isa(plottedData, "string") % convert to numeric if the annotation is a string
            plottedData = str2double(plottedData);
        end

        % Bar plot
        subplot(numFactors, numModes, plotIterator); hold on;
        plot = bar(plottedData, "FaceColor", "flat");

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


varExp = varExps{numFactors}(choice);
sgtitle(overallTitle + ", variance explained =  " + round(varExp,2) + "%");
saveas(gcf, path)
close()
end