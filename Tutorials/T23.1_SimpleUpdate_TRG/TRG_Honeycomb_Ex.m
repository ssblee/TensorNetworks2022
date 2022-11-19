function [Nsq,varargout] = TRG_Honeycomb_Ex (A,B,rgstep,Nkeep,varargin)
% < Description >
%
% Nsq         = TRG_Honeycomb_Ex (A,B,rgstep,Nkeep)
% [Nsq,OOval] = TRG_Honeycomb_Ex (A,B,rgstep,Nkeep,OA,OB)
%
% Compute the squared norm (e.g., < \Psi | \Psi >) per site, of a projected
% entangled-pair state (PEPS) on an exponentially large honeycomb lattice
% with periodic boundary conditions, by using the tensor renormalization
% group (TRG) method. Here the PEPS is defined by a unit cell of two sites,
% represented by local ket tensors A and B, which correspond to two
% bipartition sublattices, respectively.
% If local operators OA and OB are also given, this function computes the
% nearest-neighbour correlator of OA and OB, with respect to the PEPS.
%
% < Input >
% A, B : [rank-4 tensors] State tensors for each sublattice. The leg order
%       convention is counter-clockwise: 
%
%     \2            3/
%      \   1    1   /    (Physical legs are the last, i.e., the 4th;
%       A -------- B      they are hidden for brevity)
%      /            \ 
%     /3            2\
%
%       The n-th leg of A contracts with the n-th leg of B for n = 1,2,3.
% rgstep : [integer] # of RG steps. Here the original lattice has
%       6*(3^rgstep) number of sites.
% Nkeep : [numeric] Maximum bond dimension.
% OA, OB : [rank-2 or rank-3 tensors] (Optional) The operators acting onto
%       nearest-neighbor sites. OA (OB) acts on a site whose ket state is
%       represented by A (B). The first (second) legs of OA and OB act on
%       the bra (ket) tenor. If they are rank-3, their third legs are to be
%       contracted each other.
%
% < Output >
% Nsq : [numeric] The squared norm of the PEPS, per lattice site. The
%       total norm of the PEPS is given by Nsq^(6*(3^rgstep)).
% OOval : [numeric] (Optional) Only for when OA and OB are given. The
%       correlation function of local operators OA and OB acting onto
%       nearest-neighbor sites.
%
% Written by S.Lee (Jul.18,2017)
% Updated by S.Lee (Jul.08,2019): Revised for SoSe 2019.
% Updated by S.Lee (Nov.18,2022): Revised for the course at SNU.

tobj = tic2;
isOp = false; % default: no operator -> not compute correlator

% parsing optional input
if nargin > 4
    isOp = true;
    OA = varargin{1};
    OB = varargin{2};
end

% sanity check
if ndims(A) > 4
    error('ERR: Tensor A should be rank-4.');
elseif ndims(B) > 4
    error('ERR: Tensor B should be rank-4.');
end
szA = [size(A),ones(1,4-ndims(A))];
szB = [size(B),ones(1,4-ndims(B))];
if ~all(szA(1:3) == szB(1:3))
    error('ERR: The bond dimensions of A and B do not match.');
end

disp(['TERG on a honeycomb lattice: rgstep = ',sprintf('%g',rgstep), ...
    ', Nkeep = ',sprintf('%g',Nkeep)]);

% make reduced tensor
z = 3; % coordination number = number of nearest neighbours
ids = [(1:z);z+(1:z)]; % leg permutation
gA = contract(A,z+1,z+1,conj(A),z+1,z+1,ids(:));
gB = contract(B,z+1,z+1,conj(B),z+1,z+1,ids(:));

dimA = ones(2,z); dimB = ones(2,z); % bond dimensions
dimA(1:ndims(gA)) = size(gA);
dimB(1:ndims(gB)) = size(gB);

gA = reshape(gA,prod(dimA,1)); % make fat legs
gB = reshape(gB,prod(dimB,1));

if isOp
    % % % % TODO (start) % % % %

    % obtain "impurity" tensors that are made of rank-4 ket/bra 
    % tensors and local operators
    

    % % % % TODO (end) % % % %
    
    Timp = {gOA,gOB}; % impurity tensors
end

Nsq = 1; % initialize with 1

disptime('Start');

for it1 = (1:rgstep)
    Tnorm = sqrt(norm(gA(:))*norm(gB(:)));
    Nsq = Nsq*exp(log(Tnorm)/(3^(it1-1))); % each tensor (gA, gB, ...) represents 3^(it1-1) number of lattice sites
    
    % normalize tensors
    gA = gA/Tnorm;
    gB = gB/Tnorm;
    if isOp
        for it2 = (1:numel(Timp))
            Timp{it2} = Timp{it2}/Tnorm; % use the same normalization as the bulk tensors
        end
    end
    
    % % % % TODO (start) % % % %
    % contraction of gA and gB, along three different directions
    

    % SVD the contractions of gA and gB. And take the square root of
    % singular values, and contract the square root with the isometries.
    

    % contract the tensors around a plaquette to make the coarse-grained uniform tensor
    

    % after each RG step, the whole lattice rotated clockwise by 90 degree
    
    if isOp
        % when local operators to be measured are given, coarse-grain
        % "impurity" tensors also
        

    end
    % % % % TODO (end) % % % %
    
    disptime(['#',sprintf('%02i/%02i',[it1,rgstep])]);
end

% % % % TODO (start) % % % %
% now we have six tensors; contract exactly



if isOp
    % perform the contraction of six tensors, two of which are replaced
    % with Timp{1} and Timp{2}
    
    varargout{1} = 
end
% % % % TODO (end) % % % %

strs = ['Squared norm per site = ',sprintf('%.6g',Nsq)];
if isOp
    strs = [strs,', Correlation func. = ',sprintf('%.6g',varargout{1})];
end
disptime(strs);
toc2(tobj,'-v');

end