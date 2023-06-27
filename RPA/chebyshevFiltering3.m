function filteredVectors = chebyshevFiltering3(X, Hfun, lambda_cutoff, eigmin, eigmax, chebyshevDegree)
    e = (eigmax - lambda_cutoff)/2;
    c = (lambda_cutoff + eigmax)/2;
    sigma = e/(c - eigmin);
    tau = 2/sigma;
    Y = (Hfun(X) - c*X) * (sigma/e);
    for time = 1:chebyshevDegree
        sigmaNew = 1 / (tau - sigma);
        Yt = (Hfun(Y) - c*Y) * (2*sigmaNew/e) - (sigma*sigmaNew)*X;
        X = Y;
        Y = Yt;
        sigma = sigmaNew;
    end
    filteredVectors = Y;
end