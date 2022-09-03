%% Non-interacting fermions on a tight-binding chain
% Author: <https://cqm.snu.ac.kr/ Seung-Sup Lee>
%% 
% Here we consider non-interacting spinless fermions (or spin-polarized fermions, 
% equivalently) on a tight-binding chain. Its Hamiltonian is given by
% 
% $$\hat{H} = \sum_{\ell = 1}^{L-1} ( -t_\ell \hat{c}_{\ell+1}^\dagger \hat{c}_{\ell} 
% - t_\ell^* \hat{c}_{\ell}^\dagger \hat{c}_{\ell+1} ),$$
% 
% where the chain has $L$ sites, $t_\ell$ indicates the hopping amplitute between 
% sites $\ell$ and $\ell+1$, and $\hat{c}_\ell^\dagger$ creates a particle at 
% a site $\ell \in [1, L]$.
%% Exercise (a): Compute the energy and degeneracy of the many-body ground states
% Write a script or function that computes the ground-state energy and degeneracy 
% of this non-interacting tight-binding chain. The script or function takes the 
% following input and output:
% 
% *< Input >*
%% 
% * |t| : [numeric vector] Each element |t(l)| indicates a hopping amplitude 
% $t_\ell$. The length of the vector |numel(t)| equals to the number of chain 
% sites minus 1.
%% 
% *< Output >*
%% 
% * |E_G| : [numeric scalar] Ground-state energy.
% * |d_G| : [numeric scalar] Ground-state degeneracy.
%% 
% Once you complete a script or function, compute the ground-state energies 
% and degeneracies for the following three cases:
% 
% (1) $L = 10$, $t_\ell = 1$ for all $\ell$'s.
% 
% (2) $L = 11$, $t_\ell = 1$ for all $\ell$'s.
% 
% (3) $L = 11$, $t_\ell = e^{\mathrm{i} \ell}$.