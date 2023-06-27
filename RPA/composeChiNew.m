function chi = composeChiNew(occPsis, occEpsilons, unoccPsis, unoccEpsilons, Nd, Ns, omega)
    chi = zeros(Nd);
    for occ = 1:Ns
        for unocc = 1:(Nd-Ns)
            vecr = conj(occPsis(:, occ)).*unoccPsis(:, unocc);
            vecrp = conj(unoccPsis(:, unocc)).*occPsis(:, occ);
            chi = chi + 2*real(vecr*vecrp' / (occEpsilons(occ) - unoccEpsilons(unocc) - omega*1i));
%             vecr = conj(allPsis(:, unocc)).*allPsis(:, occ);
%             vecrp = conj(allPsis(:, occ)).*allPsis(:, unocc);
%             chi = chi - vecr*vecrp' / (allEpsilons(unocc, unocc) - allEpsilons(occ, occ) - omega*1i); 
        end
    end
    chi = 2*chi;
end