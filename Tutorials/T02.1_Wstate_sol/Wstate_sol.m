%% [Solution] W state
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% Solution to Exercise (a): Tensor representation of the W state
% The tensor representation of the W state for a general qubit number $N$ can 
% be generated with just two lines of code!

A = zeros(2*ones(1,N)); % define a rank-N tensor
A(2.^(0:N-1)+1) = 1/sqrt(N); % assign coefficients
%% 
% Let's check whether it works.

N = 5;
A = zeros(2*ones(1,N));
A(2.^(0:N-1)+1) = 1/sqrt(N);
nnz(A) % number of nonzero elements
A(2,1,1,1,1)*sqrt(N)
A(1,2,1,1,1)*sqrt(N)
A(1,1,2,1,1)*sqrt(N)
A(1,1,1,2,1)*sqrt(N)
A(1,1,1,1,2)*sqrt(N)