function [Lambda,GA,GB,Es] = SimpleUp_Honeycomb_Ex (Lambda,GA,GB,H,Nkeep,betas)
% < Description >
%
% [Lambda,GA,GB,Es] = SimpleUp_Honeycomb_Ex (Lambda,GA,GB,H,Nkeep,betas)
%
% Find the ground state, as a projected entangled-pair state (PEPS), of a
% system on an infinite honeycomb lattice whose Hamiltonian consists of 
% nearest-neighbor interaction terms H. The state is described by a unit
% cell of two sites, and represented by the Gamma-Lambda notation
% generalized to two dimensions. Similarly to the infinite time-evolving
% block decimation (iTEBD) algorithm of finding the ground-state of
% infinite one-dimensional system, we use the imaginary time evolution to
% find the ground state.
%
% < Input >
% Lambda : [cell array] Each cell Lambda{n} is the vector of singular
%       values that describe the bonds connecting Gamma tensors.
% GA, GB : [rank-4 tensor] Gamma tensors that have physical legs onto A and
%       B sublattice sites, respectively. The detailed description of
%       Lambda{n}, GA, and GB is as below.
%
%  \1                                 2/
%   \                                 /
%   Lambda{2}                      Lambda{3}
%     \2                           1/
%      \2  1                   1   /3       (Physical legs of GA and GB are
%       GA ----- Lambda{1} ----- GB           the last, i.e., the 4th; they
%      /3     2             1      \2         are hidden for brevity)
%     /2                            \1
%   Lambda{3}                      Lambda{2}
%   /                                 \
%  /1                                  \2
%
%       The numbers next to the legs are the leg order. The first legs of
%       Lambda{..} head towards GB, and the second legs towards GA. The
%       legs of GA and GB are ordered in counter-clockwise way, and the
%       n-th legs of GA and GB contract together.
% H : [rank-4 tensor] Two-site interaction term. The leg order convention
%       is:
%      2      4
%      ^      ^
%      |      |
%   [     H      ]
%      |      |
%      ^      ^
%      1      3
% Nkeep : [numeric] Maximum bond dimension.
% betas : [numeric vector] The vector of imaginary time step sizes.
%
% < Output >
% Lambda, GA, GB : The result tensors after imaginary time evolution.
% Es : [matrix] Measured value of two-site term H across each bond *before*
%       applying the imaginary time evolution to the bond. Each column is
%       the value at each iteration, each row is the value for each
%       Lambda{n}.
%
% Written by S.Lee (Jul.08,2019)
% Updated by S.Lee (Nov.18,2022): Revised for the course at SNU.

tobj = tic2;

% sanity check
if ndims(GA) > 4
    error('ERR: Tensor GA should be rank-4.');
elseif ndims(GB) > 4
    error('ERR: Tensor GB should be rank-4.');
elseif ~iscell(Lambda) || (numel(Lambda) ~= 3)
    error('ERR: Lambda should be a cell aray with 3 elements.');
end
szA = [size(GA),ones(1,4-ndims(GA))];
szB = [size(GB),ones(1,4-ndims(GB))];
if ~all(szA(1:3) == szB(1:3))
    error('ERR: The bond dimensions of GA and GB do not match.');
end
for itl = (1:numel(Lambda))
    if ~isvector(Lambda{itl})
        error(['ERR: Lambda{',sprintf('%i',itl),'} is not a vector.']);
    elseif numel(Lambda{itl}) ~= szA(itl)
        error(['ERR: Length of Lambda{',sprintf('%i',itl), ...
            '} does not match with the corresponding bond dimensions of GA and GB']);
    end
end
% % % 

sz = [size(H),ones(1,4-ndims(H))]; % local space dimensions
z = numel(Lambda); % coordination number = number of nearest neighbors

% diagonalize the Hamiltonian to exponentiate
Hmat = reshape(permute(H,[1 3 2 4]),[sz(1)*sz(3) sz(2)*sz(4)]); % matrix representation
[VH,DH] = eig((Hmat+Hmat')/2);
DH = diag(DH);

% record of energy measurement
Es = zeros(numel(betas),z);

disptime('Start');

for it1 = (1:numel(betas))
    % imaginary time Trotter step
    expH = reshape(VH*diag(exp(-betas(it1)*DH))*VH',sz([1 3 2 4]));

    for it2 = (1:z) % index for the bond (or Lambda) to which the imaginary
                    % time evolution is applied
        % % % % TODO (start) % % % %

        % contract Lambda's and Gamma's
        

        % Bond projection: should be no truncation here!
        

        % imaginary time Trotter step & SVD

        % divide the singluar values (that are associated with the other
        % bonds than Lambda{it2})
        
        
        % % % % TODO (end) % % % %
    end
    
    if (mod(it1,3e3) == 0) || (it1 == numel(betas))
        disptime(['#',sprintf('%i/%i',[it1 numel(betas)]), ...
            ', Measured energy = ',sprintf('%.6g',Es(it1,it2))]);
    end
end

toc2(tobj,'-v');

end