function [ts,M,Ovals,EE] = TDVP_2site (M,Hs,O,Nkeep,dt,tmax,varargin)
% < Description >
%
% [ts,M,Ovals,EE] = TDVP_2site (M,Hs,O,Nkeep,dt,tmax [,'Krylov',nKrylov] [,'tol',tol])
%
% Two-site time-dependent variational principle (TDVP) method for
% simulating real-time evolution of matrix product state (MPS). The
% expectation values of local operator O for individual sites are evaluated
% for discrete time instances. The 1D chain system is described by the
% matrix product operator (MPO) Hamiltonian H.
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
% Written by S.Lee (Nov.04,2022).

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
disp('2-site TDVP : Real-time evolution with local measurements');
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
    % permute left<->right, to make use of updateLeft
    T = permute(M{itN},[2 1 3]);
    Hlr{itN+1} = updateLeft(Hlr{itN+2},3,T,permute(Hs{itN},[1 2 4 3]),4,T);
end

disptime('Time evolution: start.');

for itt = (1:Nstep)
    % left -> right
    for itN = (1:(N-2))
        % time evolution of rank-4 tensor M{itN}*M{itN+1}, via DVP_2site_expHA2
        Aold = contract(M{itN},3,2,M{itN+1},3,1,[1 3 2 4]);
        Anew = TDVP_2site_expHA2 (Hlr{itN},Hs{itN},Hs{itN+1},Hlr{itN+3},Aold,dt/2,nKrylov,tol); % half time step dt/2
        
        % SVD; now we can safely truncate small singular values;
        [M{itN},S2,M{itN+1}] = svdTr(Anew,4,[1 3],Nkeep,1e-8);
        S2 = S2/norm(S2); % normalize
        M{itN} = permute(M{itN},[1 3 2]);
        M{itN+1} = contract(diag(S2),2,2,M{itN+1},3,1);
        
        % entanglement entropy
        Spart = -(S2.^2).*log(S2.^2)/log(2);
        EE(itt*2-1,itN) = sum(Spart(~isnan(Spart)));
        
        % update the Hamiltonian in effective basis
        Hlr{itN+1} = updateLeft(Hlr{itN},3,M{itN},Hs{itN},4,M{itN});
        
        % inverse time evolution of rank-3 tensor M{itN+1}, via TDVP_2site_expHA
        M{itN+1} = TDVP_2site_expHA (Hlr{itN+1},Hs{itN+1},Hlr{itN+3},M{itN+1},-dt/2,nKrylov,tol);  % half time step -dt/2
    end
    
    % right -> left
    for itN = ((N-1):-1:1)
        % time evolution of rank-4 tensor M{itN}*M{itN+1}, via TDVP_2site_expHA2
        Aold = contract(M{itN},3,2,M{itN+1},3,1,[1 3 2 4]);
        Anew = TDVP_2site_expHA2 (Hlr{itN},Hs{itN},Hs{itN+1},Hlr{itN+3},Aold, ...
            dt/2*((itN==(N-1))+1),nKrylov,tol); % full time step dt at the right end, otherwise dt/2
        
        % SVD; now we can safely truncate small singular values;
        [M{itN},S2,M{itN+1}] = svdTr(Anew,4,[1 3],Nkeep,1e-8);
        S2 = S2/norm(S2); % normalize
        M{itN} = contract(M{itN},3,3,diag(S2),2,1,[1 3 2]);
    
        % entanglement entropy
        Spart = -(S2.^2).*log(S2.^2)/log(2);
        EE(itt*2,itN) = sum(Spart(~isnan(Spart)));

        % update Hlr{itN+2} in effective basis
        % permute left<->right, to make use of updateLeft
        T = permute(M{itN+1},[2 1 3]);
        Hlr{itN+2} = updateLeft(Hlr{itN+3},3,T,permute(Hs{itN+1},[1 2 4 3]),4,T);
        
        if itN > 1
            % inverse time evolution of rank-3 tensor M{itN}, via TDVP_2site_expHA
            M{itN} = TDVP_2site_expHA (Hlr{itN},Hs{itN},Hlr{itN+2},M{itN},-dt/2,nKrylov,tol);  % half time step -dt/2
        end
    end

% % % % TODO (end) % % % %

    % Measurement of local operators O; currently M is in site-canonical
    % with respect to site 1
    MM = []; % contraction of bra/ket tensors from left
    for itN = (1:N)
        Ovals(itt,itN) = trace(updateLeft(MM,2,M{itN},O,2,M{itN})); % the right part is the identity
        MM = updateLeft(MM,2,M{itN},[],[],M{itN});
    end
    
    if (mod(itt,round(Nstep/10)) == 0) || (itt == Nstep)
        disptime(['#',sprintf('%i/%i',[itt,Nstep]), ...
            ' : t = ',sprintf('%g/%g',[ts(itt),ts(end)])]);
    end
end

toc2(tobj,'-v');
chkmem;

end


function Anew = TDVP_2site_expHA (Hleft,Hcen,Hright,Aold,dt,nKrylov,tol)
% < Description >
%
% Anew = TDVP_2site_expHA (Hleft,Hcen,Hright,Aold,dt,nKrylov,tol)
%
% Time-evolve a rank-3 tensor by time step dt with respect to the effective
% Hamiltonian. This subfunction is taken from the Lanczos subfunction of
% the 1-site TDVP, 'DMRG/TDVP_1site.m'.
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
[Vkr,Ekr] = eig(Hkrylov);

Anew = contract(As(:,:,:,1:cnt),4,4,Vkr*(diag(exp((-1i*dt)*diag(Ekr)))*Vkr(1,:)'),2,1);

end


function Anew = TDVP_2site_expHA2 (Hleft,Hcen1,Hcen2,Hright,Aold,dt,nKrylov,tol)
% < Description >
%
% Anew = TDVP_2site_expHA2 (Hleft,Hcen1,Hcen2,Hright,Aold,dt,nKrylov,tol)
%
% Time-evolve a rank-4 two-site tensor by time step dt with respect to the
% effective Hamiltonian. This subfunction is based on the Lanczos
% subfunctions of 'DMRG/DMRG_GS_1site.m' and 'DMRG/TDVP_1site.m'.
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
% Anew : [rank-4 tensor] The state after the time evolution.
% Enew : [numeric] Expectation value of the effective Hamiltonian with
%       respect to "Anew".

% % % % TODO (start) % % % %

% The collection of rank-4 tensors as Krylov basis; the 5th dimension is
% for indexing different rank-4 tensors
As = zeros([size(Aold,1),size(Aold,2),size(Aold,3),size(Aold,4),nKrylov]);

% Define the first Krylov vector as the input "Aold"
As(:,:,:,:,1) = Aold/sqrt(abs(contract(conj(Aold),4,(1:4),Aold,4,(1:4))));
% normalize; insert "abs" to avoid numerical noise

alphas = zeros(nKrylov,1); % main diagonal elements
betas  = zeros(nKrylov-1,1); % +-1 diagonal elements
cnt = 0; % counter for Krylov subspace dimension

for itn = (1:nKrylov)
    % "matrix-vector" multiplication
    Amul = contract(Hleft,3,2,As(:,:,:,:,itn),4,1);
    Amul = contract(Amul,5,[2 4],Hcen1,4,[3 2]);
    Amul = contract(Amul,5,[3 5],Hcen2,4,[2 3]);
    Amul = contract(Amul,5,[2 5],Hright,3,[2 3],[1 4 2 3]);

    alphas(itn) = real(contract(conj(As(:,:,:,:,itn)),4,(1:4),Amul,4,(1:4)));
    % insert "real" to avoid numerical noise

    cnt = cnt+1;

    if itn < nKrylov
        % orthogonalize, to get the next Krylov vector
        for it2 = (1:2) % do twice to further reduce numerical noise
            T = contract(conj(As(:,:,:,:,1:itn)),5,(1:4),Amul,4,(1:4));
            T = contract(As(:,:,:,:,1:itn),5,5,T,2,1);
            Amul = Amul - T;
        end
    
        % norm
        Anorm = sqrt(abs(contract(conj(Amul),4,(1:4),Amul,4,(1:4)))); % insert "abs" to avoid numerical noise
        
        if Anorm < tol % for numerical stability
            break;
        end
    
        As(:,:,:,:,itn+1) = Amul/Anorm; % normalize
        betas(itn) = Anorm;
    end
end

Hkrylov = diag(betas(1:cnt-1),-1);
Hkrylov = Hkrylov + Hkrylov' + diag(alphas(1:cnt));
[Vkr,Ekr] = eig(Hkrylov);

Anew = contract(As(:,:,:,:,1:cnt),5,5,Vkr*(diag(exp((-1i*dt)*diag(Ekr)))*Vkr(1,:)'),2,1);

% % % % TODO (end) % % % %

end