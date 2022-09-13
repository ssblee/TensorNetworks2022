function [M,S,dw] = canonForm_Ex (M,id,Nkeep,Skeep)
% < Description >
%
% [M,S,dw] = canonForm_Ex (M,id,Nkeep,Skeep);
%
% Obtain the canonical forms of MPS, depending on the index id of the
% target bond. The left part of the MPS, M{1}, ..., M{id}, is brought into
% the left-canonical form and the right part of the MPS; M{id+1}, ...,
% M{end}, into the right-canonical form. Thus, if id is 0, the result is
% purely right-canonical form; if id is numel(M), the result is purely
% left-canonical form.  
%
% < Input >
% M : [cell array] MPS of length numel(M). Each cell element is a rank-3
%       tensor, where the first, second, and third dimensions are
%       associated with left, right, and bottom (i.e., physical) legs,
%       respectively.
% id : [integer] Index for the bond connecting the tensors M{id} and
%       M{id+1}. With respect to the bond, the tensors to the left
%       (right) are brought into the left-(right-)canonical form.
% Nkeep : [number] Maximal number of singular values to keep at each SVD.
%       If set empty ([]), it is interpreted as Inf, meaning no truncation
%       by the number of singular values.
% Skeep : [number] Minimum magnitude of the singluar value to keep at each
%       SVD. If set empty ([]), it is interpreted as 10*eps(S(1)), where
%       S(1) means the largest singular value. Not to truncate by the
%       magnitude of singular values, set Skeep = 0.
%
% < Output >
% M : [cell array] Left-, right-, or bond-canonical form from input M,
%       depending on id, as follows:
%       * id == 0: right-canonical form
%       * id == numel(M): left-canonical form
%       * otherwise: bond-canonical form
% S : [column vector] Singular values at the bond between M{id} and M{id+1}
%       if 0 < id < numel(M); the norm of the MPS if id = 1 or numel(M). 
% dw : [column vector] Vector of length numel(M)-1. dw(n) means the
%       discarded weight (i.e., the sum of the square of the singular  
%       values that are discarded) at the bond between M{n} and M{n+1}.
%
% Written by S.Lee (Apr.30,2019)
% Rewritten by S.Lee (Sep.12,2022)

try

% % check the integrity of input
if (numel(id) ~= 1) || (round(id) ~= id)
    error('ERR: 2nd input ''id'' needs to be a single integer.');
elseif (id < 0) || (id > numel(M))
    error('ERR: the 2nd input ''id'' needs to be in a range (0:numel(M))');
elseif size(M{1},1) ~= 1
    error('ERR: the first dimension (= left leg) of M{1} should be of size 1.');
elseif size(M{end},2) ~= 1
    error('ERR: the second dimension (= right leg) of M{end} should be of size 1.');
end
% % % %

dw = zeros(numel(M)-1,1); % discarded weights

% % % TODO (start) % % %

% % Bring the left part of MPS into the left-canonical form
for it = (1:id-1)

    % reshape M{it} and SVD


    % contract S and Vd with M{it+1}
    
    
end
    
% % Bring the right part into the right-canonical form
for it = (numel(M):-1:id+2)
    % % % TODO (start) % % %
    
    % reshape M{it} and SVD


    % contract U and S with M{it-1}

    
end

]
if id == 0 % purely right-canonical form
    % no trucation in SVD

    % U is a single number which serves as the overall phase factor to the
    % total many-site state. So we can pass over U to M{1}.
    
elseif id == numel(M) % purely left-canonical form
    % no trucation in SVD

    % V' is a single number which serves as the overall phase factor to the
    % total many-site state. So we can pass over V' to M{end}.

else % bond-canonical form
    % SVD with truncation

end

% % % TODO (end) % % %

catch e
    disp(getReport(e));
    disp('Something got wrong. Debug!')
    keyboard
end

end