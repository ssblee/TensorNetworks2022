function [M,E0,Eiter,Sv,dw] = DMRG_GS_CBE_Ex (M,Hs,Nkeep,Nsweep,delta,varargin)
% < Description >
%
% [M,E0,Eiter,Sv,dw] = DMRG_GS_CBE_Ex (M,Hs,Nkeep,Nsweep,delta [,'Krylov',nKrylov] [,'tol',tol])
%
% Single-site density-matrix renormalization group (DMRG) complemented with
% the controlled bond expansion (CBE), for ground state search. For the
% details of the CBE-DMRG algorithm, refer to Gleis2022a [A. Gleis, J.-W.
% Li, and J. von Delft, arXiv:2207.14712 (2022)], especially the section
% "Controlled bond expansion" and Fig. 2.
%
% < Input >
% M : [1 x N cell array] Matrix product state (MPS) as the initial guess of
%       the variational search. Each M{n} is a rank-3 tensor acting on site
%       n. Their legs are ordered as left-right-bottom(=physical). The
%       length M, i.e., numel(Hs), defines the system length N.
% Hs : [1 x N cell array] Matrix product operator (MPO) of the Hamiltonian.
%       Each Hs{n} is a rank-4 tensor acting on site n. The order of legs
%       of Hs{n} is bottom-top-left-right, where the bottom (top) leg
%       contracts to the physical leg of bra (ket) tensor. 
% Nkeep : [numeric] Maximum bond dimension of the matrix product state
%       (MPS) to keep.
% Nsweep : [numeric] Number of round trips of sweeping. There will be
%       Nsweep pairs of right-to-left sweep and left-to-right sweep. That
%       is, the first sweep is right-to-left and the last is left-to-right.
% delta : [numeric] Ratio of (dimension of the relevant subspace of the
%       discarded space, to be appended to the kept bond space)/Nkeep,
%       which serves as a knob of the CBE scheme. Note that in the original
%       CBE scheme, the ratio's denominator is defined differently; it
%       involves another parameter \alpha and differs from Nkeep. In this
%       tutorial code we simply choose Nkeep as the denominator, for
%       simplicity. For details, see the description of parameter \delta in
%       the section "Sweeping" of Gleis2022a.
%
% < Option >
% 'Krylov', .. : [numeric] The maximum dimension of the Krylov subspace to
%       be considered, within the Lanczos method used for updating each MPS
%       tensor. That is, a tridiagonal matrix of size n-by-n (at maximal)
%       is diagonalized in each tensor update, with n set by 'Krylov',.. .
%       (Default: 5)
% 'tol', .. : [numeric] The tolerance for elements on the +-1 diagonals
%       (those next to the main diagonal) within the Lanczos method. If an
%       element is smaller than the tolerance, then the Lanczos
%       tridiagonalization procedure stops and only the part of the
%       tridiagonal matrix constructed so far is diagonalized.
%       (Default: 1e-8)
%
% < Output >
% M : [1 x N cell array] The result MPS which is obtained variationally to
%       have the minimum expectation value of the Hamiltonian Hs. It is in
%       a left-canonical form, since the last sweep is left-to-right.
% E0 : [numeric] The energy of M.
% Eiter : [(N-1) x (2*Nsweep) numeric array] Each element Eiter(m,n) means
%       the variational energy in the m-th iteration within the n-th sweep.
%       Odd n is for right-to-left sweep and even n for left-to-right
%       sweep. Note that the iteration index m matches with the index of
%       the site just left to the orthogonality center for left-to-right
%       sweep; the iteration m corresponds to the site (N-m) for
%       right-to-left sweep.
% Sv : [1 x (N+1) cell array] Sv{n} contains a vector of singular values on
%       the bond between sites n-1 and n, for 1 < n < N. Sv{1} and Sv{end}
%       are the norms of the MPS, sitting on the left- and rightmost legs
%       (that are dummy), respectively.
% dw : [(N-1) x (2*Nsweep) array] dw(n,m) is the discarded weight, i.e.,
%       the sum of the squares of the discarded singular values, for the
%       bond space sitting between M{n-1} and M{n}, in the m-th sweep. The
%       sum along the MPS during the last sweep, sum(dw(:,end)), indicates
%       the discarded weight for the output MPS. See Eq. (58) of
%       Schollwoeck2011 [U. Schollwoeck, Ann. Phys. 326, 96 (2011)] for how
%       to use the discarded weights as an error estimate.
%
% Written by S.Lee (Nov.08,2022): for the lecture course at SNU.

tobj = tic2;

% default vales of optional input parameters
nKrylov = 5;
tol = 1e-8;

% parse options
while ~isempty(varargin)
    if numel(varargin) < 2
        error('ERR: Option should be set by a pair of option name and value.');
    end
    switch varargin{1}
        case 'Krylov'
            nKrylov = varargin{2};
            varargin(1:2) = [];
        case 'tol'
            tol = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: Unknown input.');
    end
end

% % sanity check for input and option
N = numel(M);
if N < 2
    error('ERR: chain is too short.');
elseif N ~= numel(Hs)
    error('ERR: The lengths of ''M'' and ''Hs'' should be equal.');
end

for itN = (1:N)
    if ~all(size(M{itN},3) == [size(Hs{itN},1), size(Hs{itN},2)])
        error(['ERR: The first or second leg of Hs{',sprintf('%i',itN), ...
            '} has dimension inconsistent with the physical leg of an MPS tensor.']);
    end
end
% % %

% show message
disptime('CBE-DMRG: ground state search');
disptime(['# of sites = ',sprintf('%i',numel(Hs)), ...
    ', Nkeep = ',sprintf('%i',Nkeep),', # of sweeps = ',sprintf('%i',Nsweep),' x 2']);

% % result arrays
% ground-state energy for each iteration
Eiter = zeros(N-1,2*Nsweep);
% later, Eiter(end,end) will be taken as the final result E0
Sv = cell(1,N+1); % collection of singular value vectors
dw = zeros(N-1,2*Nsweep); % discarded weights

% bring into a site-canonical form without extra truncation
[M,Sv{N}] = canonForm(M,N-1,[],-1);

% Contractions of MPO and MPS tensors that represent the effective
% Hamiltonian for the left/right parts of the chain: When the orthogonality
% center (= central tensor in a site-canonical form) is at site n, then
% Hlr{m+1} for m < n (m > n) is the left (right) part of the effective
% Hamiltonian made of M(1:m) and Hs(1:m) [M(m:end) and Hs(m:end)].
Hlr = cell(1,N+2);
Hlr{1} = 1; % to initialize "zipper" contraction
Hlr{end} = 1; % to represent the Hamiltonian for the dummy space

% Since M is in left-canonical form by now, we build Hlr{..} by contracting
% left-to-right.
for itN = (1:N)
    Hlr{itN+1} = updateLeft(Hlr{itN},3,M{itN},Hs{itN},4,M{itN});
end

for itS = (1:Nsweep)
    % % % % TODO (start) % % % %
    % right -> left
    for itN = (N:-1:2)
        % shrewed selection
        
        % Use the subfunction eigs_1site_GS, put at the end of this file,
        % to obtain the ground state via the Lanczos method
        
        % SVD
        
        % update the next tensors
        
        % update the Hamiltonian in effective basis
        
        % % % % TODO (end) % % % %
    end
    
    % display information of the sweep
    disptime(['Sweep #',sprintf('%i/%i',[2*itS-1,2*Nsweep]),' (right -> left) : Energy = ', ...
        sprintf('%.7g',Eiter(N-1,2*itS-1))]);
    
    % left -> right
    for itN = (1:(N-1))
        % % % % TODO (start) % % % %
        % shrewed selection: exchange left<->right to recycle CBE_shrewd
        % implemented for the right-to-left sweep
        

        % SVD


        % update the next tensors


        % update the Hamiltonian in effective basis

    end

    % % % % TODO (end) % % % %

    % display informaiton of the sweep
    disptime(['Sweep #',sprintf('%i/%i',[2*itS,2*Nsweep]),' (left -> right) : Energy = ', ...
        sprintf('%.7g',Eiter(N-1,2*itS))]);
end

E0 = Eiter(N-1,2*itS); % take the last value
    
toc2(tobj,'-v');

end

function [Anew,Enew] = eigs_1site_GS (Hleft,Hcen,Hright,Aold,nKrylov,tol)
% < Description >
%
% Anew = eigs_1site_GS (Hleft,Hcen,Hright,Aold,nKrylov,tol)
%
% Update an MPS tensor acting on a single site, by solving the effective
% Hamiltonian via the Lanczos method. This is the same as the Lanczos
% sub-routine of the 1-site DMRG function DMRG/DMRG_GS_1site.m.
%
% < Input >
% Hleft : [rank-3 tensor] Left part of the effective Hamiltonian. Its legs
%       are ordered as bottom-top-right.
% Hcen : [rank-4 tensor] Center part of the effective Hamiltonian. It is
%       indeed an MPO tensor for the current site.
% Hright : [rank-3 tensor] Right part of the effective Hamiltonian. Its
%       legs are ordered as bottom-top-left.
% Aold : [rank-3 tensor] Current ket tensor.
% nKrylov, tol : Options for the Lacnzos method. See the documentation for
%       the parent function for details.
%
% < Output >
% Anew : [rank-3 tensor] The approximate ground state of the effective
%       Hamiltonian obtained by the Lanczos method, using the Krylov
%       subspace with dimension up to "nKrylov".
% Enew : [numeric] Expectation value of the effective Hamiltonian with
%       respect to "Anew".

% The collection of rank-3 tensors as Krylov basis; the 4th dimension is
% for indexing different rank-3 tensors
As = zeros([size(Aold,1),size(Aold,2),size(Aold,3),nKrylov]);

% Define the first Krylov vector as the input "Aold"
As(:,:,:,1) = Aold/sqrt(abs(contract(conj(Aold),3,(1:3),Aold,3,(1:3))));
% normalize; insert "abs" to avoid numerical noise

alphas = zeros(nKrylov,1); % main diagonal elements
betas  = zeros(nKrylov-1,1); % +-1 diagonal elements
cnt = 0; % counter for Krylov subspace dimension

for itn = (1:nKrylov)
    % "matrix-vector" multiplication
    Amul = contract(Hleft,3,2,As(:,:,:,itn),3,1);
    Amul = contract(Amul,4,[2 4],Hcen,4,[3 2]);
    Amul = contract(Amul,4,[2 4],Hright,3,[2 3],[1 3 2]);

    alphas(itn) = real(contract(conj(As(:,:,:,itn)),3,(1:3),Amul,3,(1:3)));
    % insert "real" to avoid numerical noise

    cnt = cnt+1;

    if itn < nKrylov
        % orthogonalize, to get the next Krylov vector
        for it2 = (1:2) % do twice to further reduce numerical noise
            T = contract(conj(As(:,:,:,1:itn)),4,(1:3),Amul,3,(1:3));
            T = contract(As(:,:,:,1:itn),4,4,T,2,1);
            Amul = Amul - T;
        end
    
        % norm
        Anorm = sqrt(abs(contract(conj(Amul),3,(1:3),Amul,3,(1:3)))); % insert "abs" to avoid numerical noise
        
        if Anorm < tol % for numerical stability
            break;
        end
    
        As(:,:,:,itn+1) = Amul/Anorm; % normalize
        betas(itn) = Anorm;
    end
end

Hkrylov = diag(betas(1:cnt-1),-1);
Hkrylov = Hkrylov + Hkrylov' + diag(alphas(1:cnt));

[V,D] = eig(Hkrylov);
[~,minid] = min(diag(D));
Anew = contract(As(:,:,:,1:cnt),4,4,V(:,minid),2,1);

% compute the epectation value of the effective Hamiltonian with respect to "Anew"
Amul = contract(Hleft,3,2,Anew,3,1);
Amul = contract(Amul,4,[2 4],Hcen,4,[3 2]);
Amul = contract(Amul,4,[2 4],Hright,3,[2 3],[1 3 2]);
Enew = real(contract(conj(Anew),3,(1:3),Amul,3,(1:3)));

end


function Atr = CBE_shrewd (HL,HR,Al,Ar,Hl,Hr,Sb,Nkeep,delta)
% < Description >
%
% Atr = CBE_shrewd (HL,HR,Al,Ar,Hl,Hr,Sb,delta)
%
% Perform the "shrewd selection" of identifying the relevant subspace of
% the discarded bond space, following the recipe described in Fig. 2 of
% Gleis2022a. The input tensors are assumed to be contracted as follows.
% 
%         1    2  1          2  1    2
%    /------ Al ---- diag(Sb) ---- Ar ------\
%    |       |3                    3|       |
%    |       |                      |       |
%   2|  3  3 |2                    2|  4  3 |2 
%    HL ---- Hl ----          ---- Hr ----- HR
%   1|      1|   4              3   |1      |1
%    |       |                      |       |
%
% This diagram corresponds to the diagram on the upper-left corner of Fig.
% 2 of Gleis2022a. As explained in the paper, this setup is to expand the
% isometry Al, for the update of Ar during a right-to-left sweep. For
% left-to-right sweeps, one can use this function with left<->right
% inversions.
%
% < Input >
% HL, HR : [rank-3 tensors] Effective Hamiltonian of the left and right
%       parts of the system, respectively. Their leg orders are bottom-top
%       -side (right for HL, left for HR).
% Al, Ar : [rank-3 tensors] Ket tensors of the MPS. Al and Ar are left- and
%       right-normalized, respectively. Their leg orders are left-right-
%       bottom(physical).
% Hl, Hr : [rank-4 tensors] Tensors from the MPO Hamiltonian. Their leg
%       orders are bottom-top-left-right.
% Sb : [column vector] Singular values sitting on the bond between Al and
%       Ar.
% Nkeep, delta : [numeric] Numerical parameters. See the description of the
%       parameters given in the main function.


D = numel(Sb); % MPS bond dimension
w = size(Hl,4); % MPO bond dimension
Dp = ceil(D/w); % bond dimension D' = D/w arising at an indermediate step
Dt = ceil(Nkeep*delta); % bond dimension \tilde{D} as the dimension for the truncated complement

% % % % TODO (start) % % % %


% % % % TODO (end) % % % %

end

