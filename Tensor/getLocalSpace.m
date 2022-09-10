function varargout = getLocalSpace (varargin)
% < Description >
%
% [S,I] = getLocalSpace ('Spin',s)         % spin
% [F,Z,I] = getLocalSpace ('Fermion')      % spinless fermion
% [F,Z,S,I] = getLocalSpace ('FermionS')   % spinful (spin-1/2) fermion 
%
% Generates the local operators as tensors. The result operators F and S
% are rank-3, whose 1st and 2nd legs are to be contracted with bra and ket
% tensors, respectively. The 3rd legs of F and S encode the flavors of the
% operators, such as spin raising/lowering/z or particle flavor.
% Basis of the output tensors depend on the input as follows:
%   * 'Spin',s: +s, +s-1, ..., -s
%   * 'Fermion': |vac>, c'|vac>
%   * 'FermionS': |vac>, c'_up|vac>, c'_down|vac>, c'_down c'_up|vac>
% Here c' means fermion creation operator.
%
% < Input >
% s : [integer or half-integer] The value of spin (e.g., 1/2, 1, 3/2, ...).
%
% < Output >
% S : [rank-3 tensor] Spin operators.
%       S(:,:,1) : spin raising operator S_+ multiplied with 1/sqrt(2)
%       S(:,:,2) : spin-z operator S_z
%       S(:,:,3) : spin lowering operator S_- multiplied with 1/sqrt(2)
%       Then we can construct the Heisenberg interaction ($\vec{S} \cdot
%       \vec{S}$) by: contract(S,3,3,conj(S),3,3) that results in
%       (S^+ * S^-)/2 + (S^- * S^+)/2 + (S^z * S^z) = (S^x * S^x) + (S^y *
%       S^y) + (S^z * S^z).
%       There are two advantages of using S^+ and S^- rather than S^x and
%       S^y: (1) more compact. For spin-1/2 case for example, S^+ and S^-
%       have only one non-zero elements while S^x and S^y have two. (2) We
%       can avoid complex number which can induce numerical error and cost
%       larger memory; a complex number is treated as two double numbers.
% I : [rank-2 tensor] Identity operator.
% F : [rank-3 tensor] Fermion annihilation operators. For spinless fermions
%       ('Fermion'), the 2nd dimension of F is singleton, and F(:,:,1) is
%       the annihilation operator. For spinful fermions ('FermionS'),
%       F(:,:,1) and F(:,:,2) are the annihilation operators for spin-up
%       and spin-down particles, respectively.
% Z : [rank-2 tensor] Jordan-Wigner string operator for anticommutation
%        sign of fermions.
%
% Written by S.Lee (May 04,2017)
% Revised by S.Lee (Sep.08,2022)

% % parsing input
if (numel(varargin) == 0) || ~any(strcmp(varargin{1},{'Spin','Fermion','FermionS'}))
    error('ERR: Input #1 should be either ''Spin'', ''Fermion'', or ''FermionS''.');
end

switch varargin{1}
    case 'Spin'
        if nargout > 2
            error('ERR: Too many outputs are requested.');
        end
        if numel(varargin) < 2
            error('ERR: For ''Spin'', input #2 is required.');
        end
        s = varargin{2};
        if (abs(2*s - round(2*s)) > eps) || (s <= 0)
            error('ERR: Input #2 for ''Spin'' should be positive (half-)integer.');
        end
        s = round(2*s)/2;
        isFermion = false;
        isSpin = true; % create S tensor
        I = eye(2*s+1);
    case 'Fermion'
        if nargout > 3
            error('ERR: Too many outputs are requested.');
        end
        isFermion = true; % create F and Z tensors
        isSpin = false;
        I = eye(2);
    case 'FermionS'
        if nargout > 4
            error('ERR: Too many outputs are requested.');
        end
        isFermion = true;
        isSpin = true;
        I = eye(4);
end
% % %

if isFermion
    if isSpin % spinful fermion
        % basis: empty, up, down, two (= c_down^+ c_up^+ |vac>)
        F = zeros(4,4,2);
        % spin-up annihilation
        F(1,2,1) = 1; 
        F(3,4,1) = -1; % -1 sign due to anticommutation
        % spin-down annihilation
        F(1,3,2) = 1; 
        F(2,4,2) = 1;

        Z = diag([1 -1 -1 1]);

        S = zeros(4,4,3);
        S(2,3,1) = 1/sqrt(2); % spin-raising operator (/sqrt(2))
        S(:,:,3) = S(:,:,1)'; % spin-lowering operator (/sqrt(2))
        % spin-z operator
        S(2,2,2) = +1/2;
        S(3,3,2) = -1/2;
    else % spinless fermion
        % basis: empty, occupied
        F = zeros(2,2,1);
        F(1,2,1) = 1;

        Z = diag([1 -1]);
    end
else % spin
    % basis: +s, +s-1, ..., -s
    Sp = (s-1:-1:-s);
    Sp = diag(sqrt((s-Sp).*(s+Sp+1)), 1); % spin raising operator

    Sz = diag(s:-1:-s); % spin-z operator

    S = cat(3,Sp/sqrt(2),Sz,Sp'/sqrt(2));
end

% assign the tensors into varargout
if isFermion
    if isSpin % spinful fermion
        varargout = {F,Z,S,I};
    else % spinless fermion
        varargout = {F,Z,I};
    end
else % spin
    varargout = {S,I};
end

end