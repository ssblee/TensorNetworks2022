function A = getIdentity(B,idB,varargin)
% < Description >
%
% # Usage 1
%
% A = getIdentity(B,idB, [,p])
%
% Obtain the identity tensor in the space of the idB-th leg of B. For
% example, consider a ket tensor B. Then A = getIdentity(B,2) results in:
%
%   1      2    1       2
%  -->- B ->--*-->- A ->--
%       |
%     3 ^
%       |
%
% Here the numbers next to the legs mean the order of legs, and * indicates
% the location where the legs will be contracted.
%
% # Usage 2:
% A = getIdentity(B,idB,C,idC [,p])
%
% Obtain the identity tensor in the direct product space of the Hilbert
% space of the idB-th leg of B and the space of the idC-th leg of C. For
% example, consider a ket tensor B and the identity operator C at local
% site. Then A = getIdentity(B,2,C,2) results in another ket tensor A:
%
%   1      2    1       2
%  -->- B ->--*-->- A ->--
%       |           |
%     3 ^         3 ^
%       |           |
%                   *
%                 2 ^           
%                   |
%                   C
%                   |
%                 1 ^
%                   |
%
% < Input >
% B, C : [numeric array] Tensors.
% idB, idC : [integer] Indices for B and C, respectively.
%
% < Option >
% p : [interger array] If the option is given, the identity tensor is
%       permuted according to the permutation indices p. 
%       (Default: not given, i.e., no permutation)
%
% < Output >
% A : [numeric array] Identity tensor. If p option is not given, the 1st
%       and 2nd legs of A correspond to the idB-th leg of B and the idC-th
%       leg of C, respectively. If the 'p' option is given, the legs are
%       permuted accordingly.
%
% Written by S.Lee (May 2, 2017)
% Documentation edited by S.Lee (Sep.08,2022)

% parsing input
if nargin > 3 % combine the spaces of two tensors
    C = varargin{1};
    idC = varargin{2};
    varargin(1:2) = [];
else % consider only one space
    C = [];
    idC = [];
end

% default of options
p = []; % permutation of the contracted tensor (default: no permutation)

% % % parsing options
while ~isempty(varargin)
    if isnumeric(varargin{1})
        p = varargin{1};
        varargin(1) = [];
    else
        error('ERR: Unkown option.')
    end
end
% % %

DB = size(B,idB);
if ~isempty(C)
    DC = size(C,idC);
    A = reshape(eye(DB*DC),[DB DC DB*DC]);
else
    A = eye(DB);
end
    
if ~isempty(p)
    if numel(p) < numel(size(A))
        error('ERR: # of elements of permutation option ''p'' is smaller that the rank of ''A''.');
    end
    A = permute(A,p);
end


end