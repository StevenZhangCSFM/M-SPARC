function deltaPsis = sternheimerSolver(DL11, DL22, DL33, DG1, DG2, DG3, Heff,...
    S, kpt_vec, spin, Qs, Ps, deltaV, psis, epsilons, Nd, Ns, omega)
    if ispc % windows
        addpath('..\src\');
    else % max/linux
        addpath('../src/');
    end
    unit = eye(Nd);
    deltaVdiag = diag(deltaV);
    deltaPsis = zeros(Nd, Ns);
    for occ = 1:Ns
        Q = reshape(Qs(:, occ), Nd, Nd);
        P = reshape(Ps(:, occ), Nd, Nd);
        lhsfun = @(x) h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,x,S,kpt_vec,spin) ...
            + (Q - epsilons(occ)*unit + 1i*omega*unit)*x; % 1i*omega at here should be negative of 1i*omega in composeLinearOperator.m
        rhs = (P - unit)*deltaVdiag*psis(:, occ);
        deltaPsi = gmres(lhsfun, rhs, [], 1e-8, 3000, S.LapPreconL, S.LapPreconU);
        deltaPsis(:, occ) = deltaPsi;
    end
end

%     testlhs2 = Q*deltaPsi;
%     testrhs2 = (P)*deltaVdiag*psis(:, occ);
%     testlhs1 = h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,deltaPsi,S,kpt_vec,spin) + (- epsilons(occ)*unit + 1i*omega*unit)*deltaPsi;
%     testrhs1 = (- unit)*deltaVdiag*psis(:, occ);
%     error1 = testlhs1 - testrhs1;
%     error2 = testlhs2 - testrhs2;