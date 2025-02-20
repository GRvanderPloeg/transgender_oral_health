function Factors=pickModel(allModels,numComponents,choiceNumber)

numModes = size(allModels,1);
Factors = cell(1,numModes);

for i=1:numModes
    Factors{i} = allModels{i,numComponents}(:,:,choiceNumber);
end