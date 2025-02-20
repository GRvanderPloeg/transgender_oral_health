function varExplained=calcVarExplained(X, model)

numModes = size(model,2);
numComponents = size(model{1}, 2);

if numModes == 3
    [A,B,C] = fac2let(model);
    M = A*krb(C,B)';
    %M = reshape(M, size(X));
elseif numModes == 4
    [A,B,C,D] = fac2let(model);
    M = (B*krb(krb(A,D), C)')';
end

Xflat = reshape(X, size(M));
notNaN = ~isnan(Xflat);
varExplained = (sumsqr(M(notNaN)) / sumsqr(Xflat(notNaN))) * 100;