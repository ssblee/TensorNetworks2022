%% Iterative diagonalization
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% 
% The Hilbert space dimension of a many-body system increases exponentially 
% with the number of single-particle bases (e.g., number of lattice/chain sites). 
% To keep the tensor size manageable, tensor-networks-based methods use various 
% ways to truncate the Hilbert space.
% 
% Iterative diagonalization is such an approach. Let's consider a one-dimensional 
% system having $N$ sites. At the $n$-th iteration of iterative diagonalization, 
% we diagonalize the Hamiltonian of the subsystem ranging from site 1 to site 
% $n$. The $N_\mathrm{keep}$ (and a few more) lowest-lying energy eigenvalues 
% and the corresponding eigenstates are kept to span the Hilbert space for the 
% enlarged subsystem, ranging from site 1 to site $n+1$, to be treated at the 
% next iteration.
% 
% Here we need to *keep all degenerate states whose energies are close to the 
% truncation threshold.* So the actual number of kept states at an iteration can 
% be larger than $N_\mathrm{keep}$, hence "a few more". The closeness to the threshold 
% is determined by the tolerance parameter |tol|, chosen in terms of numerical 
% precision; the states separated within this tolerance are regarded as degenerate. 
% The degeneracy often comes from physical symmetries, such as spin and particle-hole 
% symmetries. If we keep only the part of the degenerate states, then the Hilbert 
% space will not respect the symmetries anymore. This artificial symmetry breaking 
% would lead to qualitatively wrong result. (Of course, we can discard the degenerate 
% states altogether, rather than keeping them, which also preserves the symmetries.)
%% Exercise (a): Non-interacting tight-binding chain
% In this Exercise, we implement iterative diagonalization codes for computing 
% the ground-state energy of non-interacting spinless fermions on a tight-binding 
% chain. Since it's non-interacting, we can compare the iterative diagonalization 
% results with numerically exact single-particle calculations (namely, via |nonIntTB|). 
% The Hamiltonian is given by
% 
% $$\hat{H} = \sum_{\ell = 1}^{N-1} ( -t_\ell \hat{c}_{\ell+1}^\dagger \hat{c}_{\ell} 
% - t_\ell^* \hat{c}_{\ell}^\dagger \hat{c}_{\ell+1} ),$$
% 
% where the chain has $N$ sites and $\hat{c}_\ell^\dagger$ creates a particle 
% at a site $\ell \in [1, N]$. Here we consider two different chain types (i) 
% with uniform hopping amplitudes, $t_\ell = 1$, and (ii) with logarithmic hopping 
% amplitudes, $t_\ell = 2^{-(\ell-1)/2}$. Compute the ground-state energies of 
% both types for various chain lengths, $N = 2, 3, \cdots, 50$, and compare them 
% the single-particle calculations. Use $N_\mathrm{keep} = 300$.