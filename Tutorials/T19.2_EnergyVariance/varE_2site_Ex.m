function [varE,varE1,varE2] = varE_2site_Ex (M,Hs)
% < Description >
%
% [varE,varE1,varE2] = varE_2site_Ex (M,Hs)
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
[M,S] = canonForm(M,1,[],[]); % no truncation

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

end

% % % % TODO (end) % % % %

varE = sum(varE1) + sum(varE2);

toc2(tobj,'-v');
chkmem;

end