function flippedModes=checkComponentNegativity(models)

numModes = size(models, 2);
numComponents = size(models{1}, 2);
numReps = size(models{1}, 3);

correctModel = {};
flippedModes = {};

% Find the "correct model" by looking at the median loading
for i=1:numModes
    correctLoadings = [];
    for j=1:numComponents
        df = reshape(models{i}(:,j,:), [], numReps);
        correctLoadings(:,j) = median(df, 2, "omitnan");
    end
    correctModel{i} = correctLoadings;
end

% Check the loadings of every bootstrapped model  against the "correct"
% version.
% If flipping the sign gets the loadings closer to the correct version, 
% record it.
for i=1:numModes
    flippedModes{i} = zeros(numReps, numComponents);
    for j=1:numComponents
        for k=1:numReps

            df = models{i}(:,j,k);
            correctLoadings = correctModel{i}(:,j);
            flippedDf = -1 * df;

            diff = sumsqr(correctLoadings - df);
            flippedDiff = sumsqr(correctLoadings - flippedDf);

            if flippedDiff < diff
                flippedModes{i}(k,j) = 1;
            end
        end
    end
end