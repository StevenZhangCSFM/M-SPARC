function Erpa = RPA_integration(Neig, omegaNumber, eigenValues, mesh01, weights)
    traces = sum(log(ones(Neig, omegaNumber) - eigenValues) + eigenValues, 1);
    % Replace $ \int_0^\infty dx f(x) $ with $ \int_0^1 dz f(1/z - 1)/z^2 $.
    Erpa = (traces ./ (mesh01.^2))*weights';
    Erpa = Erpa / (2*pi);
end