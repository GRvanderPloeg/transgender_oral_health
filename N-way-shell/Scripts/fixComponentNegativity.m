function [fixedModels, n]=fixComponentNegativity(models, negativityMatrix)

numModes = size(models, 2);
numComponents = size(models{1}, 2);
numReps = size(models{1}, 3);

fixedModels = models;
n = 0;

% Sum the negativity matrices together to obtain the number of modes that
% are flipped per repetition/component.
combinedNegativity = negativityMatrix{1};
for i=2:numModes
    combinedNegativity = combinedNegativity + negativityMatrix{i};
end

% If the number of flipped modes is even, flip them.
flippedComponents = (combinedNegativity > 0 & rem(combinedNegativity, 2) == 0);

for i=1:numReps
    for j=1:numComponents
        if flippedComponents(i,j) == true
            for k=1:numModes
                if negativityMatrix{k}(i,j) == true
                    fixedModels{k}(:,j,i) = -1 * fixedModels{k}(:,j,i);
                end
            end
            n = n + 1;
        end
    end
end