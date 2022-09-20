%% Expectation values in the AKLT state
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% Generate the AKLT state
% We will generate the AKLT state on a finite chain of length $L$. The tensor 
% at each bulk site (i.e., any site except the left- and rightmost sites) is defined 
% as 3-dimensional array, |AKLT|:

clear

AKLT = zeros(2,2,3);
% local spin S_z = +1
AKLT(1,2,1) = sqrt(2/3);
% local spin S_z = 0
AKLT(1,1,2) = -1/sqrt(3);
AKLT(2,2,2) = +1/sqrt(3);
% local spin S_z = -1
AKLT(2,1,3) = -sqrt(2/3);

L = 50; % number of sites
M = cell(1,L); % MPS
M(:) = {AKLT};
%% 
% As we are considering a finite system with open boundary condition, the left- 
% and rightmost legs of the MPS should have dimension 1, to represent a single 
% global quantum state. Therefore, we project the space of the left leg of |M{1}| 
% and the space of the right leg of |M{end}| onto the subspaces of size 1 for 
% each leg. Since the left and right legs have dimension 2, there are $4 = 2 \times 
% 2$ different global states, which are linearly independent. We'll denote them 
% by $| \psi (\alpha, \beta) \rangle$ with $\alpha = 1,2$ is chosen for the left 
% leg of |M{1}| and $\beta = 1,2$ for the right leg of |M{end}|. For example, 
% we obtain $|\psi (1,1) \rangle$ by:

M{1} = M{1}(1,:,:);
M{end} = M{end}(:,1,:);
%% 
% Let us check that the bulk tensors are *both left- and right-normalized at 
% the same time*.

% check whether left-normalized
T = contract(conj(AKLT),3,[1 3],AKLT,3,[1 3]);
disp(T - eye(size(T))); % all zeros
% check whether right-normalized
T = contract(AKLT,3,[2 3],conj(AKLT),3,[2 3]);
disp(T - eye(size(T))); % all zeros
%% 
% Note that the bra tensor is obtained as the complex conjugate of the ket tensor, 
% via |conj|. Of course, |M{1}| is only right-normalized and |M{end}| is only 
% left-normalized, since their left/right legs are projected by boundary condition.

T = contract(conj(M{1}),3,[1 3],M{1},3,[1 3]);
disp(T); % not left-normalized
T = contract(M{1},3,[2 3],conj(M{1}),3,[2 3]);
disp(T); % right-normalized
T = contract(conj(M{end}),3,[1 3],M{end},3,[1 3]);
disp(T); % left-normalized
T = contract(M{end},3,[2 3],conj(M{end}),3,[2 3]);
disp(T); % not right-normalized
%% 
% For computing the expectation values, it is convenient to normalize the MPS. 
% (Otherwise, one needs to divide the expectation value with the square of the 
% norm of the MPS.) Let's bring the MPS into its left-canonical form, without 
% loss of generality.

% transform into left-canonical form
[M,Sv] = canonForm(M,numel(M),[],0);
fprintf('Norm of MPS = %.4g\n',norm(Sv));
%% 
% Note that the norm was not unity. Now the MPS represented by the cell array 
% |M| is normalized.
%% Magnetization
% We can compute the local magnetization as the MPS expectation value of the 
% spin-$z$ operator $S_{[n]z}$ acting onto site $n = 10$. The spin-$z$ operator 
% can be obtained by using |getLocalSpace|:

[S,I] = getLocalSpace('Spin',1);
Sz = S(:,:,2); % spin-z
Sz
%% 
% Then we contract tensors in a way of "closing the zipper", from left to right. 
% We initialize the contraction result |T| as |1|. Here the spin-$z$ operator 
% is regarded as rank-2. Note that the |updateLeft| function has been slightly 
% updated to treat the case |rankC| = |rankX| = 2; see its documentation for details.

T = 1;
n = 10; % index of site on which Sz operator acts
for itN = (1:numel(M))
    if itN == n
        T = updateLeft(T,2,M{itN},Sz,2,M{itN});
    else
        T = updateLeft(T,2,M{itN},[],[],M{itN});
    end
end
fprintf('Magnetization at site %i = %.4g\n',n,T);
%% Exercise (a): Magnetization
% In the demonstration above, the code computes the magnetization at only one 
% site, by contracting all the tensors from the left end to the right end. However, 
% as the MPS was already brought into the left-canonical form, one can start from 
% the site $n$, not from the left end. Keeping this in mind, *compute the magnetization 
% for all chain sites*. And compare the results with the exact analytic results, 
% for all different boundary conditions $\alpha, \beta = 1, 2$:
% 
% $$\langle \psi(\alpha,\beta) | \hat{S}_{[n]z} | \psi(\alpha,\beta) \rangle 
% = 2 (-1)^{\alpha} \frac{ (-1/3)^{n} - (-1)^{\alpha+\beta} (-1/3)^{L - n + 1} 
% }{  1 + (-1)^{\alpha+\beta}(-1/3)^{L}  }.$$
% 
% 
%% Exercise (b): Spin-spin correlation
% Compute the correlation function between the spin-$z$ operators at nearest-neighbor 
% sites $n$ and $n+1$. Compare this with the exact analytic results:
% 
% $$\langle \psi (\alpha,\beta) | \hat{S}_{[n+1]z} \hat{S}_{[n]z} | \psi (\alpha,\beta) 
% \rangle = \frac{(-4/9) - 4 (-1)^{\alpha+\beta} (-1/3)^L}{1 + (-1)^{\alpha+\beta} 
% (-1/3)^L}.$$
% 
% Note that the correlation function does not depend on $n$.