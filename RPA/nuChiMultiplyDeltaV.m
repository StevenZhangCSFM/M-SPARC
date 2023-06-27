function nuDeltaRho = nuChiMultiplyDeltaV(DL11, DL22, DL33, DG1, DG2, DG3, Heff, ...
    S, spin, kpt, psis, weight, epsilons, Nd, Ns, gammas, omega, deltaV)
    if ispc % windows
        addpath('..\src\');
    else % max/linux
        addpath('../src/');
    end
    numberCol = size(deltaV, 2);
    deltaRho = zeros(Nd, numberCol);
    nuDeltaRho = zeros(Nd, numberCol);
    [Qws, Pws] = composeLinearOperator(psis, weight, epsilons, Nd, Ns, gammas, omega);
    phi_guess = [];
    for i = 1:numberCol
%         [deltaPsiws] = linearSolver(H, Qws, Pws, deltaV(:, i), psis, epsilons, Nd, Ns, omega); % to be modified
        [deltaPsiws] = sternheimerSolver(DL11, DL22, DL33, DG1, DG2, DG3, Heff, ...
            S, S.kptgrid(kpt, :), spin, Qws, Pws, deltaV(:, i), psis, epsilons, Nd, Ns, omega);
        deltaRho(:, i) = composeRho(psis, deltaPsiws); % the sum of \delta rho should be zero
%         nuDeltaRho(:, i) = AARSolver(-1/(4*pi)*L, deltaRho(:, i), 6, 3);
        nuDeltaRho(:, i) = aar(S.Lap_std, -4*pi*deltaRho(:, i), phi_guess, ...
            S.poisson_tol, S.MAXIT_POISSON, 0.6,0.6,7,6, S.LapPreconL, S.LapPreconU);
        nuDeltaRho(:, i) = nuDeltaRho(:, i) - sum(nuDeltaRho(:, i))/Nd;
    end
end