function savePARAFAC(X, numericModel, annotatedModel, path)

numModes = size(annotatedModel, 2);
numComponents = size(numericModel{1}, 2);
writematrix(X, path + "_input.csv");

for i=1:numModes
    writematrix(annotatedModel{i}, path + "_mode" + i + ".csv");
end

if numModes == 3
    [A,B,C] = fac2let(numericModel);
    M = A*krb(C,B)';
    for i=1:numComponents
        Mcomp = A(:,i)*krb(C(:,i),B(:,i))';
        %Mcomp = reshape(Mcomp, size(X));
        writematrix(Mcomp, path + "_component_" + i + ".csv");
    end
    
elseif numModes == 4
    [A,B,C,D] = fac2let(numericModel);
    M = (B*krb(krb(A,D), C)')';
    for i=1:numComponents
        Mcomp = (B(:,i)*krb(krb(A(:,i),D(:,i)), C(:,i))')';
        writematrix(Mcomp, path + "_component_" + i + ".csv");
        %Mcomp = reshape(Mcomp, size(X));
    end
end

%M = reshape(M, size(X));
writematrix(M, path + "_model.csv");