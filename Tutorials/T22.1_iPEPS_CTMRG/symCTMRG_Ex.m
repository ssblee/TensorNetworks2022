function [C, T, sqnorm, Cvs] = symCTMRG_Ex (M,Nkeep,varargin)
% < Description >
%
% [C, T, sqnorm, Cvs] = symCTMRG_Ex (M,Nkeep [, option])
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
% contract the left and up legs

% merge bond legs

% leg order: down - right

% % initialize the transfer ternsor at the left edge
% contract the left legs

% merge bond legs, except for the right legs

% leg order: up - down - right (for bra) - right (for ket)

% symmetrize the corner and transfer tensors
C = (C + permute(conj(C),[2 1]))/2;
T = (T + permute(conj(T),[2 1 3 4]))/2;

% result arry for eigenvalues
Cvs = nan(Nkeep,Nctm);

% counter for convergence
iscvg = false; % set true if converged

for itc = (1:Nctm)
    % corner contraction that involves one C, two T's, conj(M) and M
    

    % bring the corner contraction into a matrix form
    

    % Hermitianize and diagonalize


    % remove negative eigenvalues (that are from numerical noises) and the
    % corresponding eigenvectors


    % truncate


    % reshape into rank-4 isometry

    
    % compute the new corner tensor
    DC = DC/norm(DC); % first normalize the vector of eigenvalues
    Cnew = diag(DC);
    Cvs(1:numel(DC),itc) = DC;

    % obtain the new transfer tensor
    
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

sqnorm = 
Cvs(:,itc+1:end) = [];

% % % % TODO (end) % % % %

toc2(tobj,'-v');
chkmem;

endend