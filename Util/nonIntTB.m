function [E_G, d_G, e_1p] = nonIntTB (t)
% < Description >
%
% [E_G, d_G, e_1p] = nonIntTB (t)
%
% Compute the energy and degeneracy of the many-body ground state of
% non-interacting spinless fermions on a tight-binding chain. The input
% parameterizes the hopping amplitudes.
% 
% t : [numeric vector] Each element t(n) indicates a hopping amplitude from
%       site n to site n+1. The length of the vector numel(t) plut 1
%       defines the number of chain sites.
%
% < Output >
% E_G : [numeric scalar] Ground-state energy, given by the sum of negative
%       elements of the single-particle eigen-energies "e_1p" (see below).
% d_G : [numeric scalar] Ground-state degeneracy, given by 2^(number of
%       zero elements of "e_1p"). Here the zero elements are identified up
%       to numerical precision noise.
% e_1p : [numeric vector] Eigenvalues of the single-particle Hamiltonian in
%       an ascending order.
%
% Written by S.Lee (Sep.02,2022)

% single-particle Hamiltonian
H_1p = diag(-t,-1); % sign factor -1 in -t as convention; the 1st diagonal below the main diagonal describes the forward hoppings
H_1p = H_1p + H_1p'; % Hermitianize
e_1p = sort(eig(H_1p),'ascend');
e_1p(abs(e_1p) < 10*size(H_1p,1)*max(eps(H_1p(:)))) = 0; % elements smaller than numerical precision noise are de facto zeros

E_G = sum(e_1p(e_1p < 0)); % ground-state energy
d_G = 2^(sum(e_1p == 0)); % degeneracy

end