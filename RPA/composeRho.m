function deltaRho = composeRho(psis, deltaPsis)
    deltaRho = 4*real(sum(conj(psis).*deltaPsis, 2));
end