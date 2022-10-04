function [Lambda,Gamma,Eiter] = iTEBD_GS_Vidal_Ex (Lambda,Gamma,H,Nkeep,taus)
% < Description >
%
% [Lambda,Gamma,Eiter] = iTEBD_GS_Vidal_Ex (Lambda,Gamma,H,Nkeep,taus)
%
% The iTEBD (infinite time-evolving block decimation) method to find the
% ground state of an infinite one-dimensional system, by applying imaginary
% time evolutions. Here we consider a unit cell of two sites.
%
% < Input >
% Lambda : [1 x 2 cell] Lambda{1} and Lambda{2} contain the singular values
%       at odd and even bond, respectively, as column vectors. An odd
%       (even) bond sits just on the right of an odd (even) site.
% Gamma : [1 x 2 cell] Gamma{1} and Gamma{2} are rank-3 "Gamma" tensors for
%       odd and even sites, respectively. Their legs are ordered as left-
%       right-physical(bottom).
%       In an infinite MPS, the tensors are repeated as follows (here the
%       numbers next to the legs indicate their orders):
%
% ->-Gamma{1}->-*->-diag(Lambda{1})->-*->-Gamma{2}->-*->-diag(Lambda{2})->- 
% 1    ^     2   1                 2   1     ^    2   1                 2 
%      |3                                    |3
%
% H : [tensor] Two-site interaction Hamiltonian. Its leg convention is as
%       below:
%
%    2         4        [ 2 (4) is to be contracted with the third leg
%    ^         ^          of Gamma{1} (Gamma{2}) ]
%    |   ...   |
%   [     H     ]
%    |   ...   |
%    ^         ^
%    1         3
%
% Nkeep : [integer] Maximum bond dimension.
% taus : [numeric] Vector of imaginary time step sizes. Each "outer"
%       iteration consists of two imaginary-time evolutions, the first for
%       odd bonds and the second for even bonds. Both time evolutions
%       within the m-th outer iteration take the times size taus(m).
%
% < Output >
% Lambda, Gamma : [1 x 2 cells each] Cell arrays of Lambda and Gamma
%       tensors, repectively, after the imaginary time evolution.
% Eiter : [(numel(taus) x 2 x 2 matrix] Eiter(m,n,k) is the measured energy
%       for an odd (k = 1) or even (k = 2) bond after odd (n = 1) or
%       even (n = 2) bonds are updated, at the m-th "outer" iteration.
%
% Written by S.Lee (Jun.18,2017)
% Updated by S.Lee (Jun.19,2017)
% Updated by S.Lee (Jun.04,2019): Revised for Sose 2019.
% Updated by S.Lee (Sep.30,2022): Revised for the course at SNU.


tobj = tic2;

Lambda = Lambda(:);
Gamma = Gamma(:);
Nstep = numel(taus);
ldim = size(H,1); % local space dimension
Skeep = 1e-8;

% % % check the integrity of input
if any([numel(Lambda) numel(Gamma)] ~= 2)
    error('ERR: # of sites per unit cell should be 2.');
end

if ndims(H) > 4
    error('ERR: H should be rank-4.');
elseif any(ldim ~= [size(H,2) size(H,3) size(H,4)])
    error('ERR: All the legs of H should have the same size.');
end

for it = (1:2)
    if ~isvector(Lambda{it})
        error(['ERR: Lambda{',sprintf('%i',it),'} should be vector.']);
    elseif numel(Lambda{it}) ~= size(Gamma{it},2)
        error(['ERR: Dimensions for Lambda{',sprintf('%i',it),'} and Gamma{', ...
            sprintf('%i',it),'} do not match.']);
    elseif numel(Lambda{mod(it,2)+1}) ~= size(Gamma{it},1)
        error(['ERR: Dimensions for Lambda{',sprintf('%i',mod(it,2)+1), ...
            '} and Gamma{',sprintf('%i',it),'} do not match.']);
    elseif size(Gamma{it},3) ~= ldim
        error(['ERR: The third leg of Gamma{',sprintf('%i',mod(it)), ...
            '} should be of size equal to the leg of H.']);
    end
end
% % % 

% show information
disptime(['iTEBD ground state search: Nkeep = ',sprintf('%i',Nkeep), ...
    ', # of imag. time steps = ',sprintf('%.4g',Nstep)]);

% energy expectation value at each step
Eiter = zeros(Nstep,2,2);

% diagonalize the Hamiltonian to exponentiate
Hmat = reshape(permute(H,[1 3 2 4]),(ldim^2)*[1 1]); % matrix representation
[VH,DH] = eig((Hmat+Hmat')/2);
DH = diag(DH);

for it1 = (1:Nstep)
    % exponentiate the matrix representation of Hamiltonian
    expH = VH*diag(exp(-taus(it1)*DH))*VH';

    % reshape matrix -> rank-4 tensor
    expH = reshape(expH,ldim*ones(1,4));
    
    for it2 = (1:2) % 1 (2): update odd (even) bonds
        % % % % TODO (start) % % % %

        % Contract Lambda's and Gamma's to construct a rank-4 ket tensor to
        % be updated
    
        % Contract a two-site gate exp(-taus(it1)*H) with the ket tensor
    
        % SVD; truncate singular values smaller than Skeep (= 1e-8 by
        % default)
    
        % Normalize the singular value vector (so that the norm becomes 1)
        
        % Update Gamma{1}, Gamma{2}
        
        % Measure energy per bond; consider the following ket:
        % Lambda{2}*Gamma{1}*Lambda{1}*Gamma{2}*Lambda{2}*Gamma{1}*Lambda{1}
        %  ----------------   --------------------------  ----------------
        %      = TA                 = TB                      = TC
        % TA, TB, and TC are useful choice of tensor contractions.
        % The physical legs of TA and TB contract to the legs of the
        % two-site gate that updates an odd bond; those of TB and TC for
        % updating an even bond.
        % And assume that the tensors on the left and the right of this
        % rank-5 ket tensor are precisely left- and right-normalized,
        % respectively. Of course they are not, but Vidal has found that
        % the tensors converge to be left- and right-normalized forms, as
        % the iteration goes.
        if it2 == 1
            
        else

        end

        % contract the nearest-neighbor interaction term for an odd bond,
        % i.e., its 1st and 2nd legs act on an odd site
        
        % contract the nearest-neighbor interaction term for an even bond,
        % i.e., its 1st and 2nd legs act on an even site
        
        % compute the squared norm of the rank-5 ket tensor as a denominator for normalization
        
        % assign normalized energy values
        Eiter(it1,it2,:) = 

        % % % % TODO (end) % % % %
    end

    if (mod(it1,500) == 0) && (it1 < Nstep)
        disptime(['#',sprintf('%i/%i',[it1 Nstep]),', E = ',sprintf('%.8g',mean(Eiter(it1,end,:),3))]);
    end
end

disptime(['#',sprintf('%i/%i',[it1 Nstep]),', E = ',sprintf('%.8g',mean(Eiter(end,end,:),3))]);

% check performance
toc2(tobj,'-v');
chkmem;


end
