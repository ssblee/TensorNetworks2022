function [Lambda,Gamma] = ortho_Orus_Ex (Lambda,Gamma)
% < Description >
%
% [Lambda,Gamma] = ortho_Orus_Ex (Lambda,Gamma)
%
% Orthonormalize an infinite MPS in Vidal's Gamma-Lambda representation,
% by using the method given in Orus2008 [R. Orus & G. Vidal, Phys. Rev. B
% 78, 155117 (2008)]. After the orthonormalization, the rank-3 tensors
% Lambda*Gamma and Gamma*Lambda become left- and right-normalized,
% respectively.
%
% < Input >
% Lambda : [vector] A column vector contanining the singular values.
% Gamma : [rank-3 tensor] Its first and second legs contract to the same
%       Lambda, given as the first input to this function.
%       When the unit cell contains n sites with n > 1, there will be
%       multiple Gamma and Lambda, say Gammas{1}, ..., Gammas{n}, and
%       Lambdas{1}, ..., Lambdas{n}. Here Lambdas{k} sits between Gammas{k}
%       and Gammas{k+1}. To orthonormalize the tensors in such a case, we
%       first contract the tensors ("coarse-graining" called in Orus2008)
%       to make a single Gamma tensor as the second input to this function:
%       Gammas{1}*Lambdas{1}*Gammas{2}* ... *Gammas{n-1}*Lambdas{n-1}*Gammas{n}
%       And Lambdas{n} becomes the first input to this function.
%
% < Output >
% Lambda : [vector] A column vector contanining the singular values, after
%       the orthonormalization.
% Gamma : [rank-3 tensor] Rank-3 Gamma tensor after the orthonormalization.
%       If the unit cell is larger than one site, we contract the outputs
%       as Lambda*Gamma*Lambda, and sequentially SVD it to obtain
%       individual Gamma's and Lambda's for individual sites.
%
% Written by S.Lee (Jun.16,2017)
% Rewritten by S.Lee (Oct.03,2022): Revised for the Tensor Networks course
%       at SNU.

% sanity check of input
if ~isvector(Lambda)
    error('ERR: ''Lambda'' should be a vector.');
elseif any(numel(Lambda) ~= [size(Gamma,1) size(Gamma,2)])
    error('ERR: The dimensions of the 1st and 2nd legs of ''Gamma'' should be equal to the length of ''Lambda''.');
end

% % % % TODO (start) % % % %

% do SVD in Fig. 2(ii) of Orus2008

% contraction in Fig. 2(ii) of Orus2008

% % % % TODO (end) % % % %

end



function X = ortho_Orus_vec (M,isright)
% < Description >
%
% X = ortho_Orus_vec (M,isright)
%
% This function obtains the eigenvector V with the eigenvalue of the
% largest absolute value (called "dominant eigenvector" in Orus2008), for a
% transfer operator made of a ket tensor M. And it decomposes V = X*X',
% where V is brought into a matrix form of size D-by-D, where D is the bond
% dimension for the left and right legs of M.
% Here we follow the recipe described in Sec. II and Fig. 2 of Orus2008,
% except that we simply use the MATLAB "eigs" function to get the dominant
% eigenvector.
%
% < Input >
% M : [rank-3 tensor] A ket tensor, to be used to construct a transfer
%       operator. Its legs are ordered as left-right-bottom(physical).
% isright : [logical] If true, this function finds the right eigenvector.
%       If false, it finds the left eigenvector.
%
% < Output >
% X : [matrix] The dominant eigenvector, brought in its matrix form, is
%       decomposed as X * X'. 
%
% Written by S.Lee (Jun.16,2017)
% Rewritten by S.Lee (Oct.03,2022): Revised for the Tensor Networks course
%       at SNU.

D = size(M,1);

% transfer operator
T = contract(conj(M),3,3,M,3,3,[3 1 4 2]);

% Convert T into a matrix form
T = reshape(T,D^2*[1 1]);

if ~isright
    T = T'; % Hermitian conjugate
end

% % % % TODO (start) % % % %

% Get the dominant eigenvector by using "eigs"

% Convert the dominant eigenvector into a matrix form, and diagonalize the
% matrix
% Note: The dominant eigenvector is equivalent up to overall sign;
% therefore, multiply -1 if the trace of the eigenvector in the matrix form
% is negative

% Once the overall sign is fixed, then the eigenvalues (for the dominant
% eigenvector in the matrix form) should be non-negative; negative values
% are numerical errors. Remove the negative eigenvalues and the
% corresponding eigenvectors.

% % % % TODO (end) % % % %

end





