%% Currently, the script cannot handle k-points and spin-polarization
%% The script cannot calculate metal, fractional occupations
load("testcase/S_gga_6state.mat");
if ispc % windows
    addpath('..\src\');
else % max/linux
    addpath('../src/');
end
%% step 1: extract wave functions, their eigenvalues ane occupations from S, 
Neig = 400; % the number of eigenvalues of chi0 to be computed by PDEP
% S.psi;S.occ;S.EigVal; 
psi = S.psi;
% psi = S.psi ./ vecnorm(S.psi); % unify wave functions
occ = S.occ;
EigVal = S.EigVal;
Nd = S.N;
Nev = S.Nev;
for i = 1:Nev
    if occ(i) < 1e-5
        occPos = i;
        break
    end
end
occPsis = psi(:, 1:occPos - 1);
occEigs = EigVal(1:occPos - 1);
unoccPsis = psi(:, occPos:Nev);
unoccEigs = EigVal(occPos:Nev);
Ns = occPos - 1;
weight = S.W(1); % <psi*|psi>=1
%% step 2: extract potential vector and laplacian operators from S
% use function h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,X,S,kptvec,spin) directly
spin = 1;
ks = 1; 
kpt = 1;% spin and ks, kpt can be modified. See eigSolver.m
Heff = S.Veff(:,spin);
rng('default'); % Initialize random number generator
rng(ks+1);
%opts = struct('maxit', 10000, 'tol', 1e-6, 'p', S.Nev+10, 'v0', rand(S.N,1), 'isreal', true);
opts = struct('maxit', 100, 'tol', S.TOL_LANCZOS, 'v0', rand(S.N*S.nspinor,1));
kpt_vec = S.kptgrid(kpt,:);
[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,kpt_vec);
% Hfun = @(x) h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,x,S,kpt_vec,spin);
gammas = ones(Ns, 1);
%% Large for loop over i\omega, at first only set a value
omegaNumber = 1;
[omega_mesh, mesh01, weights] = imaginary_omega_mesh(omegaNumber);
eigenValues = zeros(Neig, omegaNumber);
for omegaIndex = 1:omegaNumber
    omega = omega_mesh(omegaIndex);
%     % define the function \chi0(iw)*\delta V (to be removed)
%     chi0 = composeChiNew(occPsis, occEigs, unoccPsis, unoccEigs, Nd, Ns, omega);
%     % confirm the correctness of the function by using \chi0 generated before (to be removed)
%     trialDeltaV = rand(Nd, 1);
%     trialDeltaV = trialDeltaV / norm(trialDeltaV);
%     [Qws, Pws] = composeLinearOperator(occPsis, weight, occEigs, Nd, Ns, gammas, omega);
%     deltaPsis = sternheimerSolver(DL11, DL22, DL33, DG1, DG2, DG3, Heff, ...
%         S, kpt_vec, spin, Qws, Pws, trialDeltaV, occPsis, occEigs, Nd, Ns, omega);
%     trialDeltaRho = composeRho(occPsis, deltaPsis);
%     chi0DeltaV_error = trialDeltaRho/weight - chi0*trialDeltaV;
%     % brutally find the matrix nu_chi0, and compute its eigenpairs (to be removed)
    Hfun = @(x) nuChiMultiplyDeltaV(DL11, DL22, DL33, DG1, DG2, DG3, Heff, ...
        S, spin, kpt, occPsis, weight, occEigs, Nd, Ns, gammas, omega, x);
%     eyeNd = eye(Nd);
%     nu_chi0 = Hfun(eyeNd);
%     [nuChiPsis_direct, nuChiEigs_direct] = eig(nu_chi0);
%     sorted_nuChiEigs_direct = sort(diag(nuChiEigs_direct));
%% step 3: define the function \nu\chi0(iw)\delta V, \delta V = [v1, v2, ...]
%     % compose linear operator P and Q
%     % Sternheimer equation solver to get \delta\psi for every occupied orbital, iterative linear solver (call gmres)
%     % compose \delta\rho
%     % solve laplacian*\delta\phi = \delta\rho, (AAR solver)
%     trial_nuChi_DeltaV = nuChiMultiplyDeltaV(DL11, DL22, DL33, DG1, DG2, DG3, Heff, ...
%         S, spin, kpt, occPsis, weight, occEigs, Nd, Ns, gammas, omega, trialDeltaV);
%     verifyDeltaRho = -1/(4*pi)*(lapVec(DL11,DL22,DL33,DG1,DG2,DG3,trial_nuChi_DeltaV,S));
%     laplaSolver_residual = trialDeltaRho - verifyDeltaRho;
%% step 4: solve eigenvalues and eigenvectors of \nu\chi0
    % define chebyshev degree and time of loop
    Ncheb = 3;
    chebyshevDegree = 3;
%% Small chebyshev filtering loop 
    for ncheb = 1:Ncheb
        % if 1st loop, lanczos to get \lambda_min, \ecut can be -0.00001
        % else \ecut = \lambda(Neig)
        if ncheb == 1
            opts.maxit = 300; % WARNING: might need more accuracy
            [eigmin_guessVec, eigmin] = eigs(Hfun, Nd, 1, 'sr', opts);
            eigmax = 0.0;
            lambda_cutoff = -0.01;
            nuChiPsis = rand(Nd, Neig);
            nuChiPsis = nuChiPsis ./ vecnorm(nuChiPsis);
        else
            lambda_cutoff = colRRChiEigs(Neig) + 0.05;
        end
        % chebyshev filtering, if 1st loop input rand vectors; else input vectors from the last loop
        [filteredVectors] = chebyshevFiltering3(nuChiPsis, Hfun, lambda_cutoff, eigmin, eigmax, chebyshevDegree);
        % orthogonalization
        orthFilteredVecs = orth(filteredVectors);
        Nev1 = size(orthFilteredVecs,2);    % WARNING: ORTH(psi) might change the size of psi, that's why we update Nev
		assert(Nev1 == Neig,'Number of states have changed within Chebyshev filtering');
        % reileigh-ritz to get reduced eigenpairs
        Vtnuchi0V = orthFilteredVecs' * Hfun(orthFilteredVecs); % Vtchi0V is very small, and it is non-symmetric
        [RRChiPsis, RRChiEigs] = eig(Vtnuchi0V); % orthFilteredVecs'*orthFilteredVecs = I after orthonormal
        colRRChiEigs = sort(real(diag(RRChiEigs)));
        % rotate reduced eigenvectors back to get real eigenvectors
        nuChiPsis = orthFilteredVecs * RRChiPsis;
        nuChiPsis = nuChiPsis ./ vecnorm(nuChiPsis);
    end
%% end of Small chebyshev loop
    % use definition of eigenvalues |Hx-\lambda x| to verify the correctness (to be removed)
    verifyPsis = Hfun(nuChiPsis);
    error = nuChiPsis*RRChiEigs - verifyPsis;
    eigenValues(:, omegaIndex) = colRRChiEigs;
%% end of Large for loop
end
%% step 5: Gauss-Lerengre integration
Erpa = RPA_integration(Neig, omegaNumber, eigenValues, mesh01, weights);