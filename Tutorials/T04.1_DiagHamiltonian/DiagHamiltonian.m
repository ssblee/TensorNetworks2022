%% Diagonalize many-body Hamiltonians
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% 
% In this tutorial, we will construct many-body Hamiltonians by using tensor 
% networks and diagonalize them. We introduce two new functions that perform the 
% tasks needed in this tutorial.
%% 
% * |Tensor/getIdentity.m|: It creates a rank-2 identity tensor that spans the 
% space of a specified leg of an input tensor, or a rank-3 identity tensor that 
% merge two specified legs of input tensors.
% * |Tensor/getLocalSpace.m|: It generates local operators for spins, spinless 
% fermions, and spinful fermions.
%% 
% Check out the documentation of these functions for details.
%% Exercise (a): Spin-1/2 Heisenberg triangle (pen-and-paper)
% In Exercises (a) and (b), we study the system of three spin-1/2's, where every 
% pair of spins interact via the Heisenberg exchange of equal strength. The system's 
% Hamiltonian is specified as
% 
% $$\hat{H} = J (\hat{\vec{S}}_1 \cdot \hat{\vec{S}}_2 + \hat{\vec{S}}_1 \cdot 
% \hat{\vec{S}}_3 + \hat{\vec{S}}_2 \cdot \hat{\vec{S}}_3).$$
% 
% Identify the eigenvalues and their degeneracies of the Hamiltonian by hand.
% 
% (_Hint_: By exploting the SU(2) spin symmetry, the Hamiltonian can be much 
% simplified.)
%% Exercise (b): Spin-1/2 Heisenberg triangle (coding)
% In this Exercise, we solve the Hamiltonian introduced in Exercise (a) above, 
% with $J = 1$. Obtain the matrix elements of the Hamiltonian in the many-body 
% state basis, by constructing the identity tensors (without trunctating them) 
% and contracting them with spin operators. Then diagonalize the Hamiltonian to 
% obtain the energy eigenvalues. Check whether your numerical result is consistent 
% with pen-and-paper calculations done for Exercise (a) above.
%% Exercise (c): Non-interacting tight-binding chain
% In this Exercise, we consider non-interacting spinless fermions (or spin-polarized 
% fermions, equivalently) on a tight-binding chain. Its Hamiltonian is given by
% 
% $$\hat{H} = \sum_{\ell = 1}^{N-1} ( -t_\ell \hat{c}_{\ell+1}^\dagger \hat{c}_{\ell} 
% - t_\ell^* \hat{c}_{\ell}^\dagger \hat{c}_{\ell+1} ),$$
% 
% where the chain has $N$ sites, $t_\ell$ indicates the hopping amplitute between 
% sites $\ell$ and $\ell+1$, and $\hat{c}_\ell^\dagger$ creates a particle at 
% a site $\ell \in [1, N]$. Consider the case of $N = 11$ and $t_\ell = e^{\mathrm{i} 
% \ell}$. Obtain the matrix elements of the Hamiltonian in the many-body state 
% basis (without truncation), and diagonalize it to identify the ground-state 
% and the lowest-excited-state energies and their degeneracies. Compare your "many-body 
% calculation" result with one from |Util/nonIntTB.m|.