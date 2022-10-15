function C = mtimes_MPO (B,A,Nkeep,Nsweep)
% < Description >
%
% C = mtimes_MPO (B,A,Nkeep,Nsweep)
%
% Variational multiplication of two matrix product operators (MPOs) B and
% A, where A multiplies to B from above (i.e., A is the top row, and B is
% the row below A). In the variational scheme, the multiplication result C
% is initialized as A, and then undergoes a two-site update, as described
% in App. D 1 of [B.-B. Chen et al., Phys. Rev. X 8, 031082 (2018)].
%
% < Input >
% B, A : [1 x L cell array] MPOs to be multiplied together. Each tensor
%       A{n} or B{n} is a rank-4 tensor acting on site n. Its legs are
%       ordered as bottom-top-left-right, where the bottom (top) leg
%       contracts to the physical leg of bra (ket) tensor. 
% Nkeep : [numeric] The maximum bond dimension for the result MPO C.
% Nsweep : [numeric] Number of round trips of sweeping. There will be
%       Nsweep pairs of left-to-right sweep and right-to-left sweep. That
%       is, the first sweep is left-to-right and the last is right-to-left.
%
% < Output >
% C : [1 x L cell array] Multiplication of A and B. Each tensor C{n}
%       follows the same leg order convention as thoes in A and B.
%
% Written by S.Lee (Oct.10,2022)

N = numel(A);

% sanity check
if N ~= numel(B)
    error('ERR: Length of two input MPOs do not match.');
end
for itN = (1:N)
    if (itN == 1) && ~all([size(A{itN},3),size(B{itN},3)] == 1)
        error('ERR: The leftmost leg of an MPO should be dummy.');
    elseif (itN == N) && ~all([size(A{itN},4),size(B{itN},4)] == 1)
        error('ERR: The rightmost leg of an MPO should be dummy.');
    elseif (itN < N) && (size(A{itN},4) ~= size(A{itN+1},3))
        error(['ERR: The fourth (= right) leg of A{',sprintf('%i',itN), ...
            '} and the third (= left) leg of A{',sprintf('%i',itN+1), ...
            '} do not have the same dimensions.']);
    elseif (itN < N) && (size(B{itN},4) ~= size(B{itN+1},3))
        error(['ERR: The fourth (= right) leg of B{',sprintf('%i',itN), ...
            '} and the third (= left) leg of B{',sprintf('%i',itN+1), ...
            '} do not have the same dimensions.']);
    elseif size(A{itN},1) ~= size(B{itN},2)
        error(['ERR: The first (= bottom) leg of A{',sprintf('%i',itN), ...
            '} and the second (= top) leg of B{',sprintf('%i',itN), ...
            '} do not have the same dimensions.']);
    end
end
% % % 

% Initialize C with A
C = A;

% Bring C into right-canonical form
% First, convert rank-4 tensors into rank-3, by merging physical legs, to
% use the canonForm function that canonicalize MPSs
Aloc = cell(1,N); % isometries for merging the bottom and top legs of MPO tensors
for itN = (1:N)
    Aloc{itN} = getIdentity(C{itN},1,C{itN},2);
    C{itN} = contract(C{itN},4,[1 2],Aloc{itN},3,[1 2]);
end
% Use canonForm for MPS
C = canonForm(C,0,Nkeep,[]);
% Bring back to rank-4
for itN = (1:N)
    C{itN} = contract(C{itN},3,3,conj(Aloc{itN}),3,3,[3 4 1 2]);
end

% Contractions of A, B, and C^\dagger. They correspond to the "effective
% Hamitonian" for the left and right parts in the ground and excited states
% search within the DMRG. As in the DMRG, we set ABC{1} and ABC{N+2} as
% 1's, to simplify the code. 
ABC = cell(1,N+2);
ABC{1} = 1;
ABC{end} = 1;
% % % % TODO (start) % % % %
% Feel free to define sub-functions for tensor contractions that appear
% multiple times in the code. By doing so, you will be able to write a
% simpler code with less probability of encountering bugs!

% contract from right, since the first sweep is left-to-right
for itN = (N:-1:1)
    ABC{itN+1} = contract_CBA(ABC{itN+2},permute(C{itN},[1 2 4 3]), ...
        permute(B{itN},[1 2 4 3]),permute(A{itN},[1 2 4 3]));
end

for itS = (1:Nsweep)
    for itN = (1:(N-2)) % left-to-right sweep
        T = contract_CBA2(ABC{itN},B{itN},A{itN},ABC{itN+3},B{itN+1},A{itN+1});
        [U,S,Vd] = svdTr(T,6,(1:3),Nkeep,[]);
        C{itN} = permute(U,[2 3 1 4]);
        C{itN+1} = contract(diag(S),2,2,Vd,4,1,[2 3 1 4]);
        ABC{itN+1} = contract_CBA(ABC{itN},C{itN},B{itN},A{itN});
    end

    for itN = ((N-1):-1:1) % right-to-left sweep
        T = contract_CBA2(ABC{itN},B{itN},A{itN},ABC{itN+3},B{itN+1},A{itN+1});
        [U,S,Vd] = svdTr(T,6,(1:3),Nkeep,[]);
        C{itN+1} = permute(Vd,[2 3 1 4]);
        C{itN} = contract(U,4,4,diag(S),2,1,[2 3 1 4]);
        ABC{itN+2} = contract_CBA(ABC{itN+3},permute(C{itN+1},[1 2 4 3]),...
            permute(B{itN+1},[1 2 4 3]),permute(A{itN+1},[1 2 4 3]));
    end
end
% % % % TODO (end) % % % %

end

function T = contract_BA (T,B,A)
% This function contracts the following:
% (numbers next the legs indicate the orders)
%
%        2 |
%          |
%       3  |  4
%   /----- A ----
%   |      | 
%   | 3    | 1
%   |      | 2
%   |  2 3 |   4
%   T ---- B ----
%   |      |
%   |      | 1
%   | 1
%   |
%   \----
%
% The resulting rank-5 tensor has its legs in the order of (1st of T)-(1st
% of B)-(4th of B)-(2nd of A)-(4th of A)

    T = contract(T,3,3,A,4,3);
    T = contract(T,5,[2 3],B,4,[3 2],[1 4 5 2 3]);
end

function T = contract_CBA (T,C,B,A)
% This function contracts the following:
% (numbers next the legs indicate the orders)
%
%        2 |
%          |
%       3  |  4
%   /----- A ----
%   |      | 
%   | 3    | 1
%   |      | 2
%   |  2 3 |   4
%   T ---- B ----
%   |      |
%   |      | 1
%   | 1    | 2
%   |    1 |   4
%   \----- C' ---
%          |
%          | 1
%
% Here C is Hermitian conjugated; that is, it is complex conjugated and its
% first and second legs are permuted. And the bottom-most leg (1st of C')
% and the top-most leg (2nd of A) are contracted.
%
% The resulting rank-3 tensor has its legs in the order of (4th of C')-(4th
% of B)-(4th of A).

    T = contract_BA(T,B,A);
    T = contract(conj(C),4,[3 1 2],T,5,[1 2 4]);
end

function T = contract_CBA2 (Tl,Bl,Al,Tr,Br,Ar)
% This function contracts the following:
% (numbers next the legs indicate the orders)
%
%        2 |        | 2
%          |        |
%       3  |   4 3  |   4
%   /----- Al ----- Ar ------\
%   |      |        |        |
%   | 3    | 1      | 1      | 3
%   |      | 2      | 2      |
%   |  2 3 |   4 3  |   4 2  |
%  Tl ---- Bl ----- Br ----- Tr
%   |      |        |        |
%   |      | 1      | 1      |
%   | 1                      | 1
%   |                        |
%   \----                ----/
%
% The resulting rank-6 tensor has its legs in the order of (1st of Tl)-(1st
% of Bl)-(2nd of Al)-(1st of Br)-(2nd of Ar)-(1st of Tr).

    Tl = contract_BA(Tl,Bl,Al);
    Tr = contract_BA(Tr,permute(Br,[1 2 4 3]),permute(Ar,[1 2 4 3]));
    T = contract(Tl,5,[3 5],Tr,5,[3 5],[1 2 3 5 6 4]);
end