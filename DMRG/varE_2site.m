function [varE,varE1,varE2] = varE_2site (M,Hs)
% < Description >
%
% [varE,varE1,varE2] = varE_2site (M,Hs)
%
% Compute the two-site energy variance of a given MPS and and an MPO
% Hamiltonian, following the recipe givin in Sec. IV of [A. Gleis, J.-W.
% Li, and J. von Delft, arXiv:2207.13161 (2022)].
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
%
% < Output >
% varE : [numeric] Sum of the all 1- and 2-site terms of the energy
%       variance.
% varE1 : [1 x N numeric array] varE1(n) is the 1-site term of the energy
%       variance that involves the MPO tensor Hs{n}.
% varE2 : [1 x (N-1) numeric array] varE2(n) is the 2-site term of the
%       energy variance that involves the MPO tensors Hs{n} and Hs{n+1}.
%
% Written by S.Lee (Nov.04,2022): the course at SNU.

tobj = tic2;

N = numel(M);

% % sanity check for input
if N < 2
    error('ERR: chain is too short.');
elseif numel(M) ~= numel(Hs)
    error('ERR: M has different lengths from that of Hs.');
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

% bring into site-canonical form at site 1
[M,S] = canonForm(M,1,[],0); % no truncation

% Contractions of MPO and MPS tensors that represent the effective
% Hamiltonian for the left/right parts of the chain: When the orthogonality
% center (= central tensor in a site-canonical form) is at site n, then
% Hlr{m+1} for m < n (m > n) is the left (right) part of the effective
% Hamiltonian made of M(1:m) and Hs(1:m) [M(m:end) and Hs(m:end)].
Hlr = cell(1,N+2);
Hlr{1} = 1; % to initialize "zipper" contraction
Hlr{end} = 1; % to represent the Hamiltonian for the dummy space

% Since M is in right-canonical form by now (except for M{1}), Hlr{..} are
% the right parts of the Hamiltonian. That is, Hlr{n+1} is the right part
% of Hamiltonian which is obtained by contracting M(n:end) with Hs(n:end).
% (Note the index for Hlr is n+1, not n, since Hlr{1} is dummy.)
for itN = (N:-1:2)
    % permute left<->right, to make use of updateLeft
    T = permute(M{itN},[2 1 3]);
    Hlr{itN+1} = updateLeft(Hlr{itN+2},3,T,permute(Hs{itN},[1 2 4 3]),4,T);
end

% result arrays
varE1 = zeros(1,N); % 1-site terms
varE2 = zeros(1,N-1); % 2-site terms

% % % % TODO (start) % % % %

for itN = (1:N)
    % % 1-site terms
    % contract left Hamiltonian, left-normalized ket tensor, and rank-4 Hamiltonian tensor
    Tl = contract(Hlr{itN},3,2,M{itN},3,1);
    Tl = contract(Tl,4,[2 4],Hs{itN},4,[3 2]); 
    % leg order: bottom of Hlr{itN}, right of M{itN}, bottom of Hs{itN}, right of Hs{itN}

    % contract Tl with bond singular value tensor
    Tls = contract(Tl,4,2,diag(S),2,1,[1 4 2 3]);
    % leg order: bottom of Hlr{itN}, right of diag(S), bottom of Hs{itN}, right of Hs{itN}

    % contract Tl with right Hamiltonian
    Tlr1 = contract(Tls,4,[2 4],Hlr{itN+2},3,[2 3]);
    % leg order: bottom of Hlr{itN}, bottom of Hs{itN}, bottom of Hlr{itN+2}

    % contract Tlr1 with the kept projector
    Tlr1k = contract(conj(M{itN}),3,[1 3],Tlr1,3,[1 2]);
    Tlr1k = contract(M{itN},3,2,Tlr1k,2,1);
    % leg order: left of M{itN}, bottom of M{itN}, bottom of Hs{itN}, bottom of Hlr{itN+2}

    varE1(itN) = sum(abs(Tlr1k(:) - Tlr1(:)).^2);

    if itN < N
        % % 2-site terms
        % contract Tls with right-normalized ket tensor and local and right
        % Hamiltonians
        Tlr2 = contract(Tls,4,2,M{itN+1},3,1);
        Tlr2 = contract(Tlr2,5,[3 5],Hs{itN+1},4,[3 2]);
        Tlr2 = contract(Tlr2,5,[3 5],Hlr{itN+3},3,[2 3]);
        % leg order: bottom of Hlr{itN}, bottom of Hs{itN}, bottom of Hs{itN+1}, bottom of Hlr{itN+3}
    
        % contract Tlr2 with the left kept projector
        Tlr2k0 = contract(conj(M{itN}),3,[1 3],Tlr2,4,[1 2]);
        Tlr2k0 = contract(M{itN},3,2,Tlr2k0,3,1);
        % leg order: left of M{itN}, bottom of M{itN}, bottom of Hs{itN+1}, bottom of Hlr{itN+3}
    
        % contract Tlr2 with the right kept projector
        Tlr20k = contract(Tlr2,4,[3 4],conj(M{itN+1}),3,[3 2]);
        Tlr20k = contract(Tlr20k,3,3,M{itN+1},3,1);
        % leg order: bottom of Hlr{itN}, bottom of Hs{itN}, bottom of M{itN+1}, right of M{itN+1}
    
        % contract Tlr2k0 with the right kept projector
        Tlr2kk = contract(Tlr2k0,4,[3 4],conj(M{itN+1}),3,[3 2]);
        Tlr2kk = contract(Tlr2kk,3,3,M{itN+1},3,1);
        % leg order: left of M{itN}, bottom of M{itN}, bottom of M{itN+1}, right of M{itN+1}
    
        varE2(itN) = sum(abs(Tlr2(:)- Tlr2k0(:) - Tlr20k(:) + Tlr2kk(:)).^2);

        % shift the orthogonality center from a bond between M{itN} and
        % M{itN+1} to the next bond between M{itN+1} and M{itN+2}
        M{itN+1} = contract(diag(S),2,2,M{itN+1},3,1);
        [M{itN+1},S,V] = svdTr(M{itN+1},3,[1 3],[],0); % no truncation
        M{itN+1} = permute(M{itN+1},[1 3 2]);
        if itN < N-1
            M{itN+2} = contract(V,2,2,M{itN+2},3,1);

            % compute Hlr{itN+3} again, since M{itN+2} is updated
            T = permute(M{itN+2},[2 1 3]);
            Hlr{itN+3} = updateLeft(Hlr{itN+4},3,T,permute(Hs{itN+2},[1 2 4 3]),4,T);
        else
            % absorb V to M{itN+1}, since V is scalar
            M{itN+1} = contract(M{itN+1},3,2,V,2,1,[1 3 2]);
        end
        
        Hlr{itN+1} = contract(conj(M{itN}),3,[1 3],Tl,4,[1 3]);
    end
end

% % % % TODO (end) % % % %

varE = sum(varE1) + sum(varE2);

toc2(tobj,'-v');
chkmem;

end