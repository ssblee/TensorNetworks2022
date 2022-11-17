function [C, T, sqnorm, Cvs] = symCTMRG (M,Nkeep,varargin)
% < Description >
%
% [C, T, sqnorm, Cvs] = symCTMRG (M,Nkeep [, option])
%
% Perform the symmetric corner transfer matrix renormalization group
% (CTMRG) for an infinite projected entangled-pair state (iPEPS) having a 1
% x 1 unit cell of tensor M. This function considers the following network:
%
%   C -------(  T  )----- C
%   |        |     |      |
%   ^ --- conj(M)--|----- ^
%   T        |   \ |      T
%   v ------------ M ---- v
%   |        |     |      |
%   C ------ (  T  )----- C
%
% Here C is the corner tensor (rank-2) and T is the transfer tensor
% (rank-4, revealing the bond legs that contract with the legs of M and its
% complex conjugate). C and T are symmetric under exchanging their legs
% according to the inversion with respect to symmetry axes, accompanied
% with complex conjugation.
%
% < Input >
% M : [rank-5 tensor] A tensor that represents the 1 x 1 unit cell of the
%       iPEPS. Its legs are ordered as left-up-physical-down-right.
% Nkeep : [numeric] Maximum dimension of the bonds connecting C and T.
% 
% < Option >
% 'Nctm', .. : [numeric] Maximum number of CTM coarse-graining moves. If
%       the vector of the eigenvalues of the corner contraction (i.e., the
%       diagonal of the updated corner tensor) does not change (up to
%       tolerance 'tol'; see below) after an iteration, then the iteration
%       of CTM moves stops.
%       (Default: 1e5)
% 'tol', .. : If the eigenvalues (sorted in descending order) of the corner
%       contraction change less than 'tol', then the iteration of CTM moves
%       stops.
%       (Default: 1e-14)
%
% < Output >
% C : [rank-2 tensor] Corner tensor. For the corner tensor on the upper
%       left, its legs are ordered as down-right. All the corner tensors
%       are related via symmetries.
% T : [rank-4 tensor] Transfer tensor. For the transfer tensor on the left
%       boundary edge, its legs are orderd as up-down-right (for bra,
%       conj(M) tensor)-right (for ket, M tensor). All the transfer tensors
%       are related vis symmetries.
% sqnorm : [numeric] The squared norm of the iPEPS per lattice site.
% Cvs : [matrix] After each CTM coarse-graining move, the corner tensor C
%       is diagonal. Cvs(:,n) contains the diagonal elements of C after the
%       n-th move.
%
% Written by S.Lee (Nov.15,2022): It's based on Jheng-Wei Li (LMU)'s
%       tutorial code.

tobj = tic2;

% default vales of optional input parameters
Nctm = 1e5;
tol = 1e-14;

% parse options
while ~isempty(varargin)
    if numel(varargin) < 2
        error('ERR: Option should be set by a pair of option name and value.');
    end
    switch varargin{1}
        case 'Nctm'
            Nctm = varargin{2};
            varargin(1:2) = [];
        case 'tol'
            tol = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: Unknown input.');
    end
end

disptime(['Symmetric CTMRG with Nkeep = ',sprintf('%i',Nkeep)]);

% % % % TODO (start) % % % %

% % initialize the corner tensor at the top-left corner
% contract the left, up, and physical legs
C = contract(conj(M),5,(1:3),M,5,(1:3),[1 3 2 4]); 
% merge bond legs
for itl = (1:2)
    C = contract(C,4-(itl-1),[1 2],getIdentity(C,1,C,2),3,[1 2]);
end
% leg order: down - right

% % initialize the transfer tensor at the left edge
% contract the left and physical legs
T = contract(conj(M),5,[1 3],M,5,[1 3],[1 4 2 5 3 6]);
% merge bond legs, except for the right legs
for itl = (1:2)
    T = contract(T,6-(itl-1),[1 2],getIdentity(T,1,T,2),3,[1 2]);
end
T = permute(T,[3 4 1 2]);
% leg order: up - down - right (for bra) - right (for ket)

% symmetrize the corner and transfer tensors
C = (C + permute(conj(C),[2 1]))/2;
T = (T + permute(conj(T),[2 1 3 4]))/2;

% result array for eigenvalues
Cvs = nan(Nkeep,Nctm);

% counter for convergence
iscvg = false; % set true if converged

for itc = (1:Nctm)
    % corner contraction that involves one C, two T's, conj(M) and M
    CTM = symCTMRG_corner (C, T, M);
    
    % bring the corner contraction into a matrix form
    sz = [size(CTM,1),size(CTM,2),size(CTM,3)];
    Cmat = reshape(CTM,prod(sz)*[1 1]);
    
    % Hermitianize and diagonalize
    [VC,DC] = eig((Cmat + Cmat')/2);
    [DC,ids] = sort(diag(DC),'descend');
    VC = VC(:,ids);

    % remove negative eigenvalues (that are from numerical noises) and the
    % corresponding eigenvectors
    oks = (DC > 0);
    DC = DC(oks);
    VC = VC(:,oks);

    % truncate
    Ntr = min(Nkeep,numel(DC));
    DC = DC(1:Ntr);
    VC = VC(:,1:Ntr);

    % reshape into rank-4 isometry
    VC = reshape(VC,[sz Ntr]);
    
    % compute the new corner tensor
    DC = DC/norm(DC); % first normalize the vector of eigenvalues
    Cnew = diag(DC);
    Cvs(1:numel(DC),itc) = DC;

    % obtain the new transfer tensor
    Tnew = contract(conj(VC),4,1,T,4,2);
    Tnew = contract(Tnew,6,[1 5],conj(M),5,[4 1]);
    Tnew = contract(Tnew,7,[4 6 1],M,5,[1 3 4]);
    Tnew = contract(Tnew,6,[2 3 5],VC,4,(1:3),[1 4 2 3]);
    % symmetrize & normalize Tnew
    Tnew = (Tnew + permute(conj(Tnew),[2 1 3 4]))/2;
    Tnew = Tnew/norm(Tnew(:));

    % check whether the corner tensor is converged; here we don't consider
    % the tranfer tensor
    if all(size(Cnew) == size(C)) && (max(abs(Cnew(:)-C(:))) < tol)
        iscvg = true;
        break;
    end

    % update C and T
    C = Cnew;
    T = Tnew;
end

if iscvg
    disptime(['Converged after #',sprintf('%02i',itc)]);
else
    disptime('Not converged');
end

sqnorm = symCTMRG_all(C,T,M)*symCTMRG_C4(C)/(symCTMRG_CTC2(C,T)^2);
Cvs(:,itc+1:end) = [];

% % % % TODO (end) % % % %

toc2(tobj,'-v');
chkmem;

end

function [res, CT] = symCTMRG_corner (C, T, M)
% < Description >
%
% [res, CT] = symCTMRG_corner (C, T, M)
%
% With the input of C, T, M (see the description of the main function for
% their meaning), it performs the following corner contraction:
%
% res =  C -------(  T  )----- 4
%        |        |     |
%        ^ --- conj(M)--|----- 5
%        T        |   \ |
%        v ------------ M ---- 6
%        |        |     |
%        |1       |2    |3   
%
% Here the numbers next to the legs indicate the leg order. It also
% computes:
%
% CT =  C ---- 1
%       |     
%       ^ ---- 3
%       T     
%       v ---- 4
%       |     
%       |2    
%
% This contraction result can be used as the input to another subfunction.
    
CT = contract(C,2,1,T,4,1); % contract C and T on the left edge
res = contract(CT,4,1,T,4,1); % contract with T on the top edge
res = contract(res,6,[2 5],conj(M),5,[1 2]);
res = contract(res,7,[2 4 5],M,5,[1 2 3],[1 3 5 2 4 6]);

end

function res = symCTMRG_all (C, T, M)
% < Description >
%
% res = symCTMRG_all (C, T, M)
%
% With the input of C, T, M (see the description of the main function for
% their meaning), it performs the following contraction of all the relevant
% tensors:
%
% res =  C -------(  T  )----- C
%        |        |     |      |
%        ^ --- conj(M)--|----- ^
%        T        |   \ |      T
%        v ------------ M ---- v
%        |        |     |      |
%        C ------ (  T  )----- C

[res, CT] = symCTMRG_corner (C, T, M);
res = contract(CT,4,[1 3 4],res,6,(1:3));
res = contract(res,4,(2:4),CT,4,[1 3 4]);
res = contract(res,2,[1 2],C,2,[1 2]);

end

function res = symCTMRG_CTC2 (C, T)
% < Description >
%
% res = symCTMRG_CTC2 (C, T)
%
% With the input of C, T (see the description of the main function for
% their meaning), it performs the following contraction:
%
% res =  C ---- C
%        |      |
%        ^ ---- ^
%        T      T
%        v ---- v
%        |      |
%        C ---- C

CT = contract(C,2,1,T,4,1); % contract C and T on the left edge
CTC = contract(CT,4,2,C,2,1);
res = contract(CTC,4,(1:4),CTC,4,(1:4));

end

function res = symCTMRG_C4 (C)
% < Description >
%
% res = symCTMRG_C4 (C)
%
% With the input of C (see the description of the main function for its
% meaning), it performs the following contraction:
%
% res =  C ---- C
%        |      |
%        C ---- C

res = trace(C*C*C*C);

end