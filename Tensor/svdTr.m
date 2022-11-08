function [U,S,Vd,dw] = svdTr (T,rankT,idU,Nkeep,Skeep)
% < Description >
%
% [U,S,Vd,dw] = svdTr (T,rankT,idU,Nkeep,Skeep)
%
% Singular value decomposition of tensor such that T = U*diag(S)*Vd. (Note
% that it is not U*S*V' as in the MATLAB built-in function 'svd'.) If
% truncation criterion Nkeep or Skeep is given, the tensors are truncated
% to discard the smallest singular values and the corresponding singular
% vectors.
%
% < Input >
% T : [tensor] Tensor.
% rankT : [number] Rank of T.
% idU : [integer vector] Indices of T to be associated with U. For example,
%       if rankT == 4 and idU == [1 3], the result U is rank-3 tensor whose
%       1st and 2nd legs correspond to the 1st and 3rd legs of T. The 3rd
%       leg of U is associated with the 1st leg of diag(S). And Vd is
%       rank-3 tensor whose 2nd and 3rd legs correspond to the 2nd and 4th
%       legs of T. Its 1st leg is associated with the 2nd leg of diag(S).
% Nkeep : [number] Maximal number of singular values to keep. If set empty
%       ([]), it is interpreted as Inf, meaning no truncation by the number
%       of singular values.
% Skeep : [number] Minimum magnitude of the singluar value to keep. If set
%       empty ([]), it is interpreted as 10*eps(S(1)), where S(1) means the
%       largest singular value. Not to truncate by the magnitude of
%       singular values, set Skeep = 0.
%
% < Output >
% U : [tensor] Tensor describing the left singular vectors. Its last leg
%       contracts with diag(S). The earlier legs are specified by input
%       idU; their order is determined by the ordering of idU.
% S : [vector] The column vector of singular values. If there were no
%       truncation, norm(S) indicates the norm of the tensor.
% Vd : [tensor] Tensor describing the right singular vectors. Its 1st leg
%       contracts with diag(S). The later legs conserve the order of the
%       legs of input T.
% dw : [tensor] Discarded weight (= sum of the square of the singular
%       values truncated).
%
% Written by S.Lee (May 22,2017)
% Rewritten by S.Lee (Sep.12,2022)

% sanity check
if rankT < ndims(T)
    error('ERR: Input ''rankT'' is smaller than the rank of other input ''T''.');
elseif ~all(any(idU(:) == (1:rankT),2))
    error('ERR: Invalid index for tensor U (e.g. out of bound, non-integer).');
end

% dimensions of T
Tsz = ones(1,rankT);
Tsz(1:ndims(T)) = size(T);

% leg indices to be associated with Vd
idV = (1:rankT);
idV(idU) = [];

% reshape to matrix form
T2 = reshape(permute(T,[idU,idV]),[prod(Tsz(idU)) prod(Tsz(idV))]);
[U2,S2,V2] = svd(T2,'econ'); % SVD
S2 = diag(S2);

% default truncation parameters
if isempty(Nkeep)
    Nkeep = Inf;
end
if isempty(Skeep) && ~isempty(S2)
    Skeep = 10*eps(S2(1));
end

oks = (S2 >= Skeep);
oks(min(numel(S2),Nkeep)+1:end) = false;

dw = sum(S2(~oks).^2);

S = S2(oks);
U = reshape(U2(:,oks),[Tsz(idU) numel(S)]);
Vd = reshape(V2(:,oks)',[numel(S) Tsz(idV)]);

end