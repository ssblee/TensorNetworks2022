function [ts,M,Ovals,EE] = TDVP_1site_Ex (M,Hs,O,Nkeep,dt,tmax,varargin)
% < Description >
%
% [ts,M,Ovals,EE] = TDVP_1site_Ex (M,Hs,O,Nkeep,dt,tmax [,'Krylov',nKrylov] [,'tol',tol])
%
% Time-dependent variational principle (TDVP) method for simulating
% real-time evolution of matrix product state (MPS). The expectation
% values of local operator O for individual sites are evaluated for
% discrete time instances. The 1D chain system is described by the matrix
% product operator (MPO) Hamiltonian H.
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
% O : [matrix] Rank-2 tensor as a local operator acting on a site. The
%       expectation value of this operator at each chain site is to be
%       computed; see the description of the output 'Ovals' for detail.
% Nkeep : [integer] Maximum bond dimension.
% dt : [numeric] Real time step size. Each real-time evolution by step dt
%       is achieved by one pair of sweeps (left-to-right and right-to-left).
% tmax : [numeric] Maximum time range.
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
% ts : [numeric] Row vector of discrete time values.
% M : [cell] The final MPS after real-time evolution.
% Ovals : [matrix] Ovals(m,n) indicates the expectation value of local
%       operator O (input) at the site n and time ts(m).
% EE : [matrix] EE(m,n) indicates the entanglement entropy (with base 2) of
%       the MPS with respect to the bipartition at the bond between the
%       sites n and n+1, after acting the m-th full sweep. Note that the
%       time evolution from time ts(k-1) to ts(k) consists of two sweeps
%       associated with the rows EE(2*k+(-1:0),:). Since the base 2 
%       is chosen, the value 1 of the entanglement entropy means one
%       "ebit".
%
% Written by S.Lee (Jun.13,2019): Written for SoSe 2019.
% Updated by S.Lee (Jun.13,2019): Revised for SoSe 2020.
% Updated by S.Lee (Oct.29,2022): Revised for the course at SNU.

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

Nstep = ceil(tmax/dt);
N = numel(M);
ts = dt*(1:Nstep).';

% % sanity check for input
if Nstep < 1
    error('ERR: No time evolution to be done? Check dt and tmax.');
end

if N < 2
    error('ERR: chain is too short.');
elseif numel(M) ~= numel(Hs)
    error('ERR: M has different lengths from that of Hs.');
elseif ~ismatrix(O)
    error('ERR: local operator O should be rank 2.');
end

for itN = (1:N)
    if size(Hs{itN},1) ~= size(Hs{itN},2)
        error(['ERR: The first and second legs of Hs{', ...
            sprintf('%i',itN),'} have different dimensions.']);
    elseif size(Hs{itN},2) ~= size(M{itN},3)
        error(['ERR: The second leg of Hs{', ...
            sprintf('%i',itN),'} and the third leg of M{', ...
            sprintf('%i',itN),'} have different dimensions.']);
    end
end
% % % 

% show message
disp('TDVP : Real-time evolution with local measurements');
fprintf(['N = ',sprintf('%i',numel(M)),', Nkeep = ',sprintf('%i',Nkeep), ...
    ', dt = ',sprintf('%.4g',dt),', tmax = ',sprintf('%g',ts(end)), ...
    ' (',sprintf('%.4g',Nstep),' steps)\n']);

% results
Ovals = zeros(Nstep,N);
EE = zeros(Nstep*2,N);

% bring into a left- and then a right-canonical form, to normalize the MPS
% and to truncate the bond dimensions at both ends
M = canonForm(M,N,Nkeep,-1); 
M = canonForm(M,0,Nkeep,-1);
% set Skeep as -1, not to truncate zero singular values, since the 1-site
% update cannot increase bond dimensions

% Contractions of MPO and MPS tensors that represent the effective
% Hamiltonian for the left/right parts of the chain: When the orthogonality
% center (= central tensor in a site-canonical form) is at site n, then
% Hlr{m+1} for m < n (m > n) is the left (right) part of the effective
% Hamiltonian made of M(1:m) and Hs(1:m) [M(m:end) and Hs(m:end)].
Hlr = cell(1,N+2);
Hlr{1} = 1; % to initialize "zipper" contraction
Hlr{end} = 1; % to represent the Hamiltonian for the dummy space

% % % % TODO (start) % % % %

% Since M is in right-canonical form by now (except for M{1}), Hlr{..} are
% the right parts of the Hamiltonian. That is, Hlr{n+1} is the right part
% of Hamiltonian which is obtained by contracting M(n:end) with Hs(n:end).
% (Note the index for Hlr is n+1, not n, since Hlr{1} is dummy.)
for itN = (N:-1:1)
    
end

disptime('Time evolution: start.');

for itt = (1:Nstep)
    % Hint: complete and use the subfunctions TDVP_1site_expHA and
    % TDVP_1site_expHC!

    % left -> right
    for itN = (1:(N-1))
        
    end
    
    itN = N; % right end
    

    % right -> left
    for itN = ((N-1):-1:1)
        
    end

    % Measurement of local operators O; currently M is in site-canonical
    % with respect to site 1
    
    
    if (mod(itt,round(Nstep/10)) == 0) || (itt == Nstep)
        disptime(['#',sprintf('%i/%i',[itt,Nstep]), ...
            ' : t = ',sprintf('%g/%g',[ts(itt),ts(end)])]);
    end
end

% % % % TODO (end) % % % %

toc2(tobj,'-v');
chkmem;


end


function Anew = TDVP_1site_expHA (Hleft,Hcen,Hright,Aold,dt,nKrylov,tol)
% < Description >
%
% Anew = TDVP_1site_expHA (Hleft,Hcen,Hright,Aold,dt,nKrylov,tol)
%
% Time-evolve a rank-3 tensor by time step dt with respect to the effective
% Hamiltonian. This subfunction is adapted from the Lanczos subfunction of
% 'DMRG/DMRG_GS_1site.m'. Here the difference is that it considers the time
% evolution, not the ground state.
%
% < Input >
% Hleft : [rank-3 tensor] Left part of the effective Hamiltonian. Its legs
%       are ordered as bottom-top-right.
% Hcen : [rank-4 tensor] Center part of the effective Hamiltonian. It is
%       indeed an MPO tensor for the current site.
% Hright : [rank-3 tensor] Right part of the effective Hamiltonian. Its
%       legs are ordered as bottom-top-left.
% Aold : [rank-3 tensor] Current ket tensor.
% dt : [numeric] Time step size.
% nKrylov, tol : Options for the Lacnzos method. See the documentation for
%       the parent function for details.
%
% < Output >
% Anew : [rank-3 tensor] The state after the time evolution.

% % % % TODO (start) % % % %


% % % % TODO (end) % % % %

end


function Anew = TDVP_1site_expHC (Hleft,Hright,Aold,dt,nKrylov,tol)
% Time-evolve a rank-2 tensor on a bond by time step dt with respect to the
% effective Hamiltonian. This subfunction is adapted from another
% subfunction "TDVP_1site_expHA" of this function. The differences here
% are: (i) the ket tensor is rank-2 and (ii) there is no the local
% Hamiltonian 'Hcen'.

% % % % TODO (start) % % % %


% % % % TODO (end) % % % %

end