%% W state
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% 
% In quantum information theory, there are several special quantum states. One 
% among them is the W state. The $N$-qubit W state is defined by the equal superposition 
% of all possible pure states in which exactly one of the qubits is in the $|1\rangle$ 
% state and the rest are in the $|0\rangle$ state: 
% 
% $$|\mathrm{W}\rangle = \frac{1}{\sqrt{N}} ( |1 0 0 \cdots 0 0 \rangle + |0 
% 1 0 \cdots 0 0\rangle + \cdots + | 0 0 0 \cdots 0 1\rangle).$$
% 
% The W state is named after *W*olfgang Dür, an Austrian physicist who proposed 
% the W state for three qubits, together with Guifré Vidal and J. Ignacio Cirac 
% in [<https://journals.aps.org/pra/abstract/10.1103/PhysRevA.62.062314 W. Dür, 
% G. Vidal, and J. I. Cirac, Phys. Rev. A *62*, 062314 (2000)>]. You will see 
% the latter two authors' names often during this lecture course, as they are 
% early founders of tensor networks!
% 
% The $N$-qubit W state can be represented in terms of a rank-$N$ tensor $A$,
% 
% $$| \mathrm{W} \rangle = | n_1 n_2 \cdots n_N \rangle A^{n_1+1, n_2+1, \cdots, 
% n_N + 1},$$
% 
% where $n_i = 0,1$ indicate the state of the $i$-th qubit, associated with 
% indices 1 and 2, respectively, along the $i$-th dimension of the tensor $A$. 
% The repeated indices $n_1, \cdots, n_N$ are assumed to be summed over.
%% Exercise (a): Tensor representation of the W state
% Write a script or function that generate the rank-$n$ tensor $A$, taking a 
% general input of $N$. Try to compose it in the most compact way, while keeping 
% its computational efficiency.