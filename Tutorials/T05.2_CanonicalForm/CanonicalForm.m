%% Canonical forms of MPS
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% Exercise (a): Complete the function that transforms MPSs into canonical forms
% There is a function |*canonForm_Ex.m*| which is zipped together with this 
% document. The function is designed to transform the input MPS |M| into left-, 
% right-, and bond-canonical forms, depending on the index |id| of the target 
% bond.
% 
% Moreover, one can also set criteria |Nkeep| and |Skeep|, for truncating small 
% singular values and the corresponding singular vectors at each SVD. (Note that 
% we do not consider degneracies near the truncation threshold in this application, 
% in contrast to iterative diagonalization.)
% 
% *Read the documentation* of the function for the details of its input and 
% output properties.
% 
% The main computational parts of |canonForm_Ex| are missing; search for the 
% keyword "TODO" to identify missing parts. *Complete the function.* You will 
% notice that the function is wrapped by |try| ... |catch e| ... |end|, which 
% is usefull for debugging.
% 
% To test whether the function is correctly implemented, let's consider the 
% following example. First, define a random MPS |M|.

clear

N = 50; % number of sites
d = 3; % local space dimension
D = 30; % bond dimension

M = cell(1,N); % MPS; M{n} is the tensor at site n

for itN = (1:N)
    % assign individual tensors
    % leg order: left, right, bottom
    if itN == 1
    % left end; left leg is of size 1
        M{itN} = rand(1,D,d);
    elseif itN == N
    % right end; right leg is of size 1
        M{itN} = rand(D,1,d);
    else
        M{itN} = rand(D,D,d);
    end
end
%% 
% Obtain the left-canonical form of |M|. We do not truncate singular values 
% and vectors, by setting |Nkeep| = |[]| and |Skeep| = |0|.

[M_L,S_L] = canonForm_Ex(M,numel(M),[],0);
S_L
%% 
% Since |M| is not normalized from the outset, the norm |S| is not 1; it's actually 
% very huge!
% 
% We see that the sizes of the leftmost tensors have decreased after bringing 
% into the left-canonical form.

M(1:5)
M_L(1:5)
%% 
% Indeed, in this example, |M(1:3)| were redundantly large, so they could be 
% compactified just via the thin SVD.
% 
% Compute the overlap between the transformed MPS (with the norm |S| pulled 
% out) and the original one. It should be equal to the norm of the MPS.

Tovl = 1; % overlap
for itN = (1:N)
    Tovl = updateLeft(Tovl,2,M{itN},[],[],M_L{itN});
end
Tovl/S_L
%% 
% We can use |updateLeft| also to "close the zipper" from right, after permuting 
% the left and right legs of the ket and bra tensors.

Tovl = 1; % overlap
for itN = (N:-1:1)
    Tovl = updateLeft(Tovl,2,permute(M{itN},[2 1 3]), ...
        [],[],permute(M_L{itN},[2 1 3]));
end
Tovl/S_L
%% 
% We can further bring |M_L| into the right-canonical form.

[M_R,S_R] = canonForm_Ex(M_L,0,[],0);
Tovl = 1; % overlap
for itN = (1:N)
    Tovl = updateLeft(Tovl,2,M{itN},[],[],M_R{itN});
end
Tovl/S_R
%% 
% Now also the right-most tensors are compactified.

M(end-4:end)
M_R(end-4:end)
%% 
% As the last check, let's consider the bond-canonical form.

[M_B,S_B] = canonForm_Ex(M,25,[],0);
Tovl = 1; % overlap
for itN = (1:N)
    Tovl = updateLeft(Tovl,2,M{itN},[],[],M_B{itN});
    if itN == 25
        Tovl = Tovl*diag(S_B);
    end
end
%% 
% Here the overlap should equal to the squared norm, since the diagonal matrix 
% of singular values is contracted in the middle.

Tovl/(S_L^2)
%% Exercise (b): Truncate bond dimensions
% *Before solving this Exercise, it is necessary to solve Exercise 1 above, 
% since we will use the function |canonForm_Ex.m| that brings the MPS into canonical 
% forms.*
% 
% In this Exercise, we will compare how the truncation of the bond space affects 
% the norm of the MPS, for different canonical forms. Let's generate a random 
% MPS again.

clear

N = 50; % number of sites
d = 3; % local space dimension
D = 30; % bond dimension

M = cell(1,N); % MPS; M{n} is the tensor at site n

for itN = (1:N)
    % assign individual tensors
    % leg order: left, right, bottom
    if itN == 1
    % left end; left leg is of size 1
        M{itN} = rand(1,D,d);
    elseif itN == N
    % right end; right leg is of size 1
        M{itN} = rand(D,1,d);
    else
        M{itN} = rand(D,D,d);
    end
end
%% 
% Make a round trip, by first bringing |M| to the left- and then to the right-canonical 
% form.

M = canonForm_Ex(M,numel(M),[],0); % left-canonical
M = canonForm_Ex(M,0,[],0); % right-canonical
%% 
% *With this right-canonical MPS |M| as the input to |canonForm_Ex.m,| apply 
% the transformation into (i) left-canonical form, (ii) right-canonical form, 
% and (iii) bond-canonical form for the bond between |M{25}| and |M{26}|, with 
% different values of |Nkeep| < |D| = |30|, while keeping |Skeep| = |0*|. How 
% does the norm change? Can you explain the result?