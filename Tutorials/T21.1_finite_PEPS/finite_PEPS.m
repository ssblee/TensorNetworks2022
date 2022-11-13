%% Contraction of finite PEPS
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% 
% In this tutorial, we will implement a scheme to measure correlation functions 
% with respect to a projected entangled-pair state (PEPS) on a square lattice 
% with open boundary conditions.
%% Resonating valence bond (RVB) state
% Let's first identify which kind of tensor networks need to be contracted to 
% measure correlation functions. For this, we take an example of the RVB state, 
% which allows compact PEPS representation, and the correlation function of spin-$z$ 
% operators at nearest neighbors.
% 
% Each PEPS tensor at a lattice site consists of (i) two valence bond (or dimer) 
% tensors defined on the right and down edges from the site and (ii) the projector 
% onto a physical space.
% 
% The valence bond tensor represents a state $| {\uparrow \downarrow} \rangle 
% - | {\downarrow \uparrow} \rangle + |0 0 \rangle$, where the first labels inside 
% $| \cdot \rangle$ indicate the states of the left (upper) "virtual particle" 
% for a horizontal (vertical) dimer. The leg order and directions are chosen as 
% left (incoming)-right (outgoing) for a horizontal dimer and up (incoming)-down 
% (outgoing) for a vertical dimer. Moreover, we use the ordered basis $\{ |{\uparrow}\rangle, 
% |{\downarrow}\rangle, |{0}\rangle \}$ for the bond spaces.

clear
VB = blkdiag([0,1;-1,0],1);
%% 
% The projector is rank-5, with legs sorted as left (incoming)-up (incoming)-physical 
% (incoming)-down (outgoing)-right(outgoing). The projector has an element 1 if 
% and only if only one bond carrys spin $| {\uparrow} \rangle$ or $| {\downarrow} 
% \rangle$ while the other three bonds are empty $| 0 \rangle$, and the bond spin 
% equals to the physical spin. Accordingly, the projector is constructed as

[S,I] = getLocalSpace('Spin',1/2);
P = zeros(3,3,2,3,3);
for it1 = (1:size(I,2))
    P(it1,3,it1,3,3) = 1;
    P(3,it1,it1,3,3) = 1;
    P(3,3,it1,it1,3) = 1;
    P(3,3,it1,3,it1) = 1;
end
%% 
% We obtain a rank-5 tensor by contracting the dimer tensors and the projector.

M = contract(P,5,4,VB,2,1);
M = contract(M,5,4,VB,2,1);
%% 
% For contractions, we derive two types of reduced tensors, one made of |M| 
% and its complex conjugate only, and the other made of |M|, its complext conjugate, 
% and the spin-$z$ operator.

MM = contract(conj(M),5,3,M,5,3);
for itl = (1:4)
    MM = contract(MM,8-(itl-1),[1 5-(itl-1)], ...
        getIdentity(MM,1,MM,5-(itl-1)),3,[1 2]);
end
MSM = contract(conj(M),5,3,S(:,:,2),2,1);
MSM = contract(MSM,5,5,M,5,3);
for itl = (1:4)
    MSM = contract(MSM,8-(itl-1),[1 5-(itl-1)], ...
        getIdentity(MSM,1,MSM,5-(itl-1)),3,[1 2]);
end
%% 
% To the spin-spin correlation function, we consider two different square-lattice 
% networks of the reduced tensors. To be more concrete, let's consider a lattice 
% of 11 rows and 12 columns, and consider putting the spin-$z$ operators at the 
% "center," i.e., sites (6,6) and (6,7). The correlation function is written as 
% $\langle \psi | \hat{S}^z_{(6,6)} \hat{S}^z_{(6,7)} | \psi \rangle / \langle 
% \psi | \psi \rangle$. First, the network for the denominator is given by

Nrow = 11;
Ncol = 12;
Tden = cell(Nrow,Ncol);
Tden(:) = {MM};
Tden(:,1) = cellfun(@(x) x(end,:,:,:), Tden(:,1), 'UniformOutput',false);
Tden(1,:) = cellfun(@(x) x(:,end,:,:), Tden(1,:), 'UniformOutput',false);
Tden(end,:) = cellfun(@(x) x(:,:,end,:), Tden(end,:), 'UniformOutput',false);
Tden(:,end) = cellfun(@(x) x(:,:,:,end), Tden(:,end), 'UniformOutput',false);
%% 
% Here we have projected the legs along the lattice boundary onto the last basis 
% of $|0 \rangle$ (for the bond space sfor kets) $\otimes | 0 \rangle$ (for the 
% bonds spaces for bras). The network for the numerator is defined in a similar 
% way,

Tnum = Tden;
Tnum(6,[6 7]) = {MSM};
%% Exercise (a): Complete the function for the MPO-MPS method of PEPS contraction
% There is a function |contract_finPEPS_Ex.m|, which is in the same sub-directory 
% with this script. It is supposed to contract tensor networks on a square lattice, 
% such as |Tden| and |Tnum| above. Complete the parts enclosed by the comments 
% |TODO (start)| and |TODO (end)|.
%% 
% Once you complete the function, you can verify your implementation by running 
% the following script.

Nkeep = 30;
Nsweep = 5;
Tnum_res = contract_finPEPS_Ex (Tnum,Nkeep,Nsweep);
disp(Tnum_res);
Tden_res = contract_finPEPS_Ex (Tden,Nkeep,Nsweep);
disp(Tden_res);
disp(Tnum_res/Tden_res);
%% 
% As we see, the contraction results can be very small or large, since each 
% tensor is not normalized and one contracts many of such tensors. Thus it is 
% crucial to properly treat the normalization within the contraction scheme.
%% Exercise (b): Spin-spin correlation function of the RVB state with different lattice sizes
% Consider the RVB state on a square lattice of "almost square" sizes, i.e., 
% $N_\mathrm{row} = N_\mathrm{col} \pm 1$, for different odd numbers $N_\mathrm{row}$. 
% And put a pair of spin-$z$ operators at the "center" of the lattice, i.e., sites 
% $((N_\mathrm{row} + 1)/2, N_\mathrm{col}/2)$, $((N_\mathrm{row} + 1)/2, N_\mathrm{col}/2 
% + 1)$. Compute the spin-spin correlation functions for different lattice sizes.
%% Exercise (c): Ground state of Kitaev's toric code
% Construct the ground state of Kitaev's toric code on a square lattice with 
% open boundary conditions, as a PEPS. Compute the correlation function of two 
% spin-$z$ operators. What do you get? Can you explain the reason?