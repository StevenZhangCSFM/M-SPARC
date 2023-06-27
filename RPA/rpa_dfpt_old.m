load("S_grid12_n300.mat");
%% step 1: extract wave functions, their eigenvalues ane occupations from S, 
% S.psi;S.occ;S.EigVal; 

%% step 2: extract potential vector and laplacian operators from S
% use function h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,X,S,kptvec,spin) directly

%% Large for loop over i\omega, at first only set a value
% and generate \chi0(iw) by definition. w1=0; w2=1.6 (to be removed in the real code)

% define the function \chi0(iw)*\delta V (to be removed)

% confirm the correctness of the function by using \chi0 generated before (to be removed)

%% step 3: define the function \nu\chi0(iw)\delta V, \delta V = [v1, v2, ...]
    % compose linear operator P and Q

    % Sternheimer equation solver to get \delta\psi for every occupied orbital, iterative linear solver (call gmres?)

    % compose \delta\rho

    % solve laplacian*\delta\phi = \delta\rho, AAR solver

%% step 4: solve eigenvalues and eigenvectors of \nu\chi0
% define chebyshev degree and time of loop
%% Small chebyshev filtering loop 
    % if 1st loop, lanczos to get \lambda_min, \ecut can be -0.00001
    % else \ecut = \lambda(Neig)

    % chebyshev filtering, if 1st loop input rand vectors; else input vectors from the last loop

    % orthogonalization

    % reileigh-ritz to get reduced eigenpairs
    
    % rotate reduced eigenvectors back to get real eigenvectors

%% end of Small chebyshev loop

% use definition of eigenvalues |Hx-\lambda x| to verify the correctness (to be removed)

%% end of Large for loop

%% step 5: Gauss-Lerengre integration