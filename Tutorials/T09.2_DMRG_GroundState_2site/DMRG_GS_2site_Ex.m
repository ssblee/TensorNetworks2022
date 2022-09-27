function [M,E0,Eiter,Sv] = DMRG_GS_2site_Ex (M,Hs,Nkeep,Nsweep,varargin)
% < Description >
%
% [M,E0,Eiter,Sv] = DMRG_GS_2site_Ex (M,Hs,Nkeep,Nsweep [,'Krylov',nKrylov] [,'tol',tol])
%
% Two-site density-matrix renormalization group (DMRG) calculation to
% search for the ground state and its energy of a one-dimensional system,
% whose Hamiltonian is given by the matrix product operator. Here we
% discard singluar values smaller than 1e-8 and their corresponding
% singular vectors.
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
%
% < Option >
% 'Krylov', .. : [numeric] The maximum dimension of the Krylov subspace to
%       be considered, within the Lanczos method used for updating MPS
%       tensors. That is, a tridiagonal matrix of size n-by-n (at maximal)
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
% Eiter : [N x (2*Nsweep) numeric array] Each element Eiter(m,n) means the
%       variational energy in the m-th iteration within the n-th sweep.
%       Odd n is for right-to-left sweep and even n for left-to-right
%       sweep. Note that the iteration index m matches with the site index
%       for left-to-right sweep; the iteration m corresponds to the site
%       (N+1-m) for right-to-left sweep.
% Sv : [1 x (N+1) cell array] Sv{n} contains a vector of singular values on
%       the bond between sites n-1 and n, for 1 < n < N. Sv{1} and Sv{end}
%       are empty (= []), since they are for the dummy legs and the
%       orthogonality center in the two-site update algorithm does not move
%       to the left- and rightmost legs.
%
% Written by S.Lee (May 28,2019)
% Updated by S.Lee (May 23,2020): Revised for SoSe2020.
% Rewritten by S.Lee (Sep.23,2022) for the 2022 fall semester at SNU.

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
disptime('Two-site DMRG: ground state search');
disptime(['# of sites = ',sprintf('%i',numel(Hs)), ...
    ', Nkeep = ',sprintf('%i',Nkeep),', # of sweeps = ',sprintf('%i',Nsweep),' x 2']);

% bring into a left-canonical form without extra truncation
M = canonForm(M,N,[],0);

% ground-state energy for each iteration
Eiter = zeros(N-1,2*Nsweep);
% later, Eiter(end,end) will be taken as the final result E0

Sv = cell(1,N+1); % collection of singular value vectors

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
    % right -> left
    for itN = (N-1:-1:1)
        % % % % TODO (start) % % % %
        % Use the subfunction eigs_2site_GS, put at the end of this file,
        % to obtain the ground state via the Lanczos method
        

        % SVD; now we can safely truncate small singular values;


        % update the next tensor

        % update the Hamiltonian in effective basis

        % permute the left and right legs, to use updateLeft
        % % % % TODO (end) % % % %
    end
    
    % display information of the sweep
    disptime(['Sweep #',sprintf('%i/%i',[2*itS-1,2*Nsweep]),' (right -> left) : Energy = ', ...
        sprintf('%.7g',Eiter(N-1,2*itS-1))]);
    
    % left -> right
    for itN = (1:N-1)
        % % % % TODO (start) % % % %


        % SVD; now we can safely truncate small singular values;


        % update the next tensor

        % update the Hamiltonian in effective basis

        % % % % TODO (end) % % % %
    end

    % display informaiton of the sweep
    disptime(['Sweep #',sprintf('%i/%i',[2*itS,2*Nsweep]),' (left -> right) : Energy = ', ...
        sprintf('%.7g',Eiter(N-1,2*itS))]);
end

E0 = Eiter(N-1,2*itS); % take the last value
    
toc2(tobj,'-v');

end

function [Anew,Enew] = eigs_2site_GS (Hleft,Hcen1,Hcen2,Hright,Aold,nKrylov,tol)
% < Description >
%
% Anew = eigs_1site_GS (Hleft,Hcen1,Hcen2,Hright,Aold,nKrylov,tol)
%
% Update an MPS tensor acting on two neighboring sites, by solving the
% effective Hamiltonian via the Lanczos method.
%
% < Input >
% Hleft : [rank-3 tensor] Left part of the effective Hamiltonian. Its legs
%       are ordered as bottom-top-right.
% Hcen1, Hcen2 : [rank-4 tensor] Center, local parts of the effective
%       Hamiltonian. They are MPO tensors for the current two sites.
% Hright : [rank-3 tensor] Right part of the effective Hamiltonian. Its
%       legs are ordered as bottom-top-left.
% Aold : [rank-4 tensor] Product of two rank-3 ket tensors. Its legs are
%       ordered as (left bond)-(right bond)-(left physical)-(right
%       physical).
% nKrylov, tol : Options for the Lacnzos method. See the documentation for
%       the parent function for details.
%
% < Output >
% Anew : [rank-3 tensor] The approximate ground state of the effective
%       Hamiltonian obtained by the Lanczos method, using the Krylov
%       subspace with dimension up to "nKrylov".
% Enew : [numeric] Expectation value of the effective Hamiltonian with
%       respect to "Anew".

% % % % TODO (start) % % % %

% The collection of rank-4 tensors as Krylov basis; the 5th dimension is
% for indexing different rank-4 tensors


% Define the first Krylov vector as the input "Aold"



for itn = (1:nKrylov)
    % "matrix-vector" multiplication

    % insert "real" to avoid numerical noise
end


% % % % TODO (end) % % % %

end