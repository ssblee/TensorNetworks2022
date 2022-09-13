function Cleft = updateLeft(Cleft,rankC,B,X,rankX,A)
% < Description >
%
% Cleft = updateLeft(Cleft,rankC,B,X,rankX,A)
%
% Contract the operator Cleft that act on the Hilbert space of the left
% part of the MPS (i.e., left of a given site) with the tensors B, X, and
% A, acting on the given site.
%
% < Input >
% Cleft : [tensor] Rank-2 or 3 tensor from the left part of the system. If
%       given as empty (i.e., []), then Cleft is considered as the identity
%       tensor of rank 2 (for rank(X) < 4) or rank 3 (for rank(X) == 4).
% rankC : [integer] Rank of Cleft. If Cleft is given as [], also set
%       rankC as [].
% B, A : [tensors] Ket tensors, whose legs are ordered as left-right-bottom
%       (local physical). In the contraction, the Hermitian conjugate
%       (i.e., bra form) of B is used, while A is contracted as it is. This
%       convention of inputting B as a ket tensor reduces extra
%       computational cost of taking the Hermitian conjugate of B.
% X : [tensor] Local operator with rank 2, 3, or 4. If given as [], then X
%       is considered as the identity. If rank-2, its legs are ordered as
%       bottom-top. If rank-3, bottom-top-right. If rank-4, bottom-top-left
%       -right.
% rankX : [integer] Rank of X. If X is given as [], also set X as [].
%
% < Output >
% Cleft : [tensor] Contracted tensor. The tensor network diagrams
%       describing the contraction are as follows.
%       * When Cleft is rank-3 and X is rank-2:
%                    1     2
%          /--------->- A ->--            /---->-- 2
%          |            | 3               |
%        2 ^            ^                 |
%          |    3       | 2               |      
%        Cleft---       X         =>    Cleft ---- 3
%          |            | 1               |
%        1 ^            ^                 |
%          |            | 3               |
%          \---------<- B'-<--            \----<-- 1
%                    1     2
%       * When Cleft is rank-2 and X is rank-3:
%                    1     2
%          /--------->- A ->--            /---->-- 2
%          |            | 3               |
%        2 ^            ^                 |
%          |          2 |   3             |      
%        Cleft          X ----    =>    Cleft ---- 3
%          |          1 |                 |
%        1 ^            ^                 |
%          |            | 3               |
%          \---------<- B'-<--            \----<-- 1
%                    1     2
%       * When both Cleft and X are rank-3:
%                    1     2
%          /--------->- A ->--            /---->-- 2
%          |            | 3               |
%        2 ^            ^                 |
%          |   3     3  | 2               |      
%        Cleft--------- X         =>    Cleft
%          |            | 1               |
%        1 ^            ^                 |
%          |            | 3               |
%          \---------<- B'-<--            \----<-- 1
%                    1     2
%       * When Cleft is rank-3 and X is rank-4:
%                    1     2
%          /--------->- A ->--            /---->-- 2
%          |            | 3               |
%        2 ^            ^                 |
%          |   3    3   | 2               |      
%        Cleft--------- X ---- 4   =>   Cleft ---- 3
%          |            | 1               |
%        1 ^            ^                 |
%          |            | 3               |
%          \---------<- B'-<--            \----<-- 1
%                    1     2
%       Here B' denotes the Hermitian conjugate (i.e., complex conjugate
%       and permute legs by [3 2 1]) of B.
%
% Written by H.Tu (May 3,2017); edited by S.Lee (May 19,2017)
% Rewritten by S.Lee (May 5,2019)
% Updated by S.Lee (May 27,2019): Case of rank-3 Cleft and rank-4 X is
%       added.
% Updated by S.Lee (Jul.28,2020): Minor fix for the case when Cleft == []
%       and rank(X) == 4.
% Updated by S.Lee (Sep.09,2022): Revised for the leg order convention used
%       for the course at SNU.

% sanity check
if isempty(Cleft)
    rankC = 2; % regard as the case of the rank-2 identity
end
if isempty(X)
    rankX = 2; % regard as the case of the rank-2 identity
end
if (isempty(rankC) || (rankC < ndims(Cleft))) || ...
        (isempty(rankX) || (rankX < ndims(X))) || ...
        ~ismember([rankC rankX],[2 2; 3 2; 2 3; 3 3; 3 4],'rows')
    error('ERR: Invalid ranks of C and X.');
end

B = conj(B); % take complex conjugate to B, without permuting legs

if ~isempty(X)
    T = contract(X,rankX,2,A,3,3);

    if ~isempty(Cleft)
        if (rankC > 2) && (rankX > 2)
            if rankX == 4 % (rankC,rankX) = (3,4)
                % contract the 3rd leg of Cleft and the 1st leg of X
                T = contract(Cleft,rankC,[rankC 2],T,rankX+1,[2 rankX]);
                Cleft = contract(B,3,[1 3],T,rankC+rankX-3,[1 2],[1 3 2]);
            else % (rankC,rankX) = (3,3)
                % contract the operator-flavor legs of Cleft and X
                T = contract(Cleft,rankC,[rankC 2],T,rankX+1,[2 rankX]);
                Cleft = contract(B,3,[1 3],T,rankC+rankX-3,[1 2]);
            end
        else % (rankC,rankX) = (2,3), (3,2)
            T = contract(Cleft,rankC,2,T,rankX+1,rankX);
            Cleft = contract(B,3,[1 3],T,rankC+rankX-1,[1 rankC],[1 3 2]);
        end
    elseif (rankX == 4) && (size(X,1) == 1) % no Cleft, rankX = 4
        % if X is rank-4 and its left leg is dummy, and Cleft is empty (i.e., identity)
        Cleft = contract(B,3,[1 3],T,rankX+1,[rankX 2],[1 4 3 2]);
        % permute the left dummy leg of X to be at the end; to be ignored
        % as singleton dimension
    else % no Cleft, rankX = 2,3
        Cleft = contract(B,3,[1 3],T,rankX+1,[rankX 1],[1 (3:rankX) 2]);
    end
elseif ~isempty(Cleft) % no X, rankC = 2, 3
    T = contract(Cleft,rankC,2,A,3,1);
    Cleft = contract(B,3,[1 3],T,rankC+1,[1 rankC+1],[1 (3:rankC) 2]);
else % no Cleft, no X
    Cleft = contract(B,3,[1 3],A,3,[1 3]);
end

end