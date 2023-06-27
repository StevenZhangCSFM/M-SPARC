function [Qs, Ps] = composeLinearOperator(psis, weight, epsilons, Nd, Ns, gammas, omega)
    Qs = zeros(Nd*Nd, Ns);
    Ps = zeros(Nd*Nd, Ns);
    for occ = 1:Ns % \delta\psi_n
        Q = zeros(Nd);
        P = zeros(Nd);
        for occ2 = 1:Ns % m
            outerProd = psis(:, occ2)*(psis(:, occ2))';
            Q = Q + gammas(occ2)*outerProd;
            if occ2 ~= occ
                P = P + gammas(occ2)/(epsilons(occ) - epsilons(occ2) - 1i*omega)*outerProd; % 1i*omega at here should be negative of 1i*omega in linearSolver.m
            else
                P = P + 1.0*weight * outerProd;
            end
        end
        Qs(:, occ) = Q(:);
        Ps(:, occ) = P(:);
    end
end