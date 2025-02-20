function sortedModel=sortParafacComponents(X, model)

numModes = size(model,2);
numComponents = size(model{1}, 2);
varExplained = 1:numComponents;

for i=1:numComponents
    fakeModel = {};
    for j=1:numModes
        fakeModel{j} = model{j}(:,i);
    end
    varExplained(i) = calcVarExplained(X, fakeModel);
end

[~, Index] = sort(varExplained, 'descend');

sortedModel = model;
for i=1:numModes
    sortedModel{i} = sortedModel{i}(:,Index);
end