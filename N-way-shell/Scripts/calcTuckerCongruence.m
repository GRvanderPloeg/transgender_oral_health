function Tuckers=calcTuckerCongruence(model)

numFactors = size(model{1},2);
numModes = size(model,2);
numCombinations = numFactors * (numFactors-1) / 2;
numCombo = 1;
Tuckers = {};
for i=1:numModes
    Tuckers{i} = zeros(1,numCombinations);
end

if numFactors == 1
    %disp('You cannot calculate similarity in a 1-component model.');
    return;
end


for f1=1:(numFactors-1)
    for f2=(f1+1):numFactors
        % tuckersA(numCombo) = sum(A(:,f1) .* A(:,f2)) / sqrt(sumsqr(A(:,f1)) * sumsqr(A(:,f2)));
        % tuckersB(numCombo) = sum(B(:,f1) .* B(:,f2)) / sqrt(sumsqr(B(:,f1)) * sumsqr(B(:,f2)));
        % tuckersC(numCombo) = sum(C(:,f1) .* C(:,f2)) / sqrt(sumsqr(C(:,f1)) * sumsqr(C(:,f2)));

        for i=1:numModes
            Tuckers{i}(numCombo) = sum(model{i}(:,f1) .* model{i}(:,f2)) / sqrt(sumsqr(model{i}(:,f1)) * sumsqr(model{i}(:,f2)));
        end
        numCombo = numCombo + 1;
    end
end