function T = contPlaq (varargin)
% < Description >
%
% T = contPlaq (A,B,C,D, ...)
%
% Contract rank-3 tensors A, B, C, ... placed around a plaquette. The
% tensors are placed around the plaquette in counter-clockwise order, and
% the legs of individual tensor are ordered in counter-clockwise order. For
% example, if there are 4 input tensors A, B, C, and D, the output T is
%
%     \2          2/
%      \   1  3   /  
%       A ------ D                  1\   /4
%      3|        |1                   \ /
%       |        |        ==>          T
%      1|        |3                   / \
%       B ------ C                  2/   \3
%      /   3  1   \
%     /2          2\
%
% < Input >
% A, B, .. : [tensors] Tensors aroung a plaquette. The order of tensors A,
%        B, .. is (counter-)clockwise and the order of the legs of each
%        tensor is also (counter-)clockwise.
%
% < Output >
% T : [tensor] Contraction result. The legs follows the same order of
%        input.
%
% Written by S.Lee (Jul.12,2017)
% Updated by S.Lee (Jul.08,2019): Documentation updated

n = nargin; % number of tensors
T = varargin{1};

for it = (2:n-1)
    T = contract(T,it+1,it+1,varargin{it},3,1);
end

T = contract(T,n+1,[1 n+1],varargin{n},3,[3 1]);

end