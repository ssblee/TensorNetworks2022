function [M,S,dw] = canonForm (M,id,Nkeep,Skeep)
% < Description >
%
% [M,S,dw] = canonForm (M,id,Nkeep,Skeep);
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

% % Bring the left part of MPS into the left-canonical form
for it = (1:id-1)
    % % % TODO (start) % % %

    % reshape M{it} and SVD
    [M{it},S,Vd,dw(it)] = svdTr(M{it},3,[1 3],Nkeep,Skeep);
    M{it} = permute(M{it},[1 3 2]);

    % contract S and Vd with M{it+1}
    M{it+1} = contract(diag(S)*Vd,2,2,M{it+1},3,1);
    
    % % % TODO (end) % % %
end
    
% % Bring the right part into the right-canonical form
for it = (numel(M):-1:id+2)
    % % % TODO (start) % % %
    
    % reshape M{it} and SVD
    [U,S,M{it},dw(it-1)] = svdTr(M{it},3,1,Nkeep,Skeep);

    % contract U and S with M{it-1}
    M{it-1} = contract(M{it-1},3,2,U*diag(S),2,1,[1 3 2]);
    
    % % % TODO (end) % % %
end

if id == 0 % purely right-canonical form
    [U,S,M{1}] = svdTr(M{1},3,1,[],0); % no trucation
    % U is a single number which serves as the overall phase factor to the
    % total many-site state. So we can pass over U to M{1}.
    M{1} = contract(U,2,2,M{1},3,1);

elseif id == numel(M) % purely left-canonical form
    [M{end},S,Vd] = svdTr(M{end},3,[1 3],[],0); % no trucation
    % Vd is a single number which serves as the overall phase factor to the
    % total many-site state. So we can pass over Vd to M{end}.
    M{end} = contract(M{end},3,3,Vd,2,1,[1 3 2]);
    
else % bond-canonical form
    T = contract(M{id},3,2,M{id+1},3,1);
    [M{id},S,M{id+1},dw(id)] = svdTr(T,4,[1 2],Nkeep,Skeep);
    M{id} = permute(M{id},[1 3 2]);
end

catch e
    disp(getReport(e));
    disp('Something got wrong. Debug!')
    keyboard
end

end