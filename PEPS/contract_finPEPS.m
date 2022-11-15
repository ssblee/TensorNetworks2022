function res = contract_finPEPS (T,Nkeep,Nsweep)
% < Description >
%
% res = contract_finPEPS (T,Nkeep,Nsweep)
%
% Contract the reduced tensors on a finite square lattice with open
% boundary conditions, by using the MPO-MPS method. That is, the first row
% of the tensor array is regarded as an MPS (with fat bonds and physical
% legs), and the second, third, etc. rows are absorbed into the firt row
% sequentially, by using the variational multiplication implemented in
% DMRG/mtimes_MPO.m. Note that in "mtimes_MPO" an MPO is multiplied with
% the other MPO, not MPS. But here an MPO whose up legs have all singleton
% dimentions can be regarded as an MPS, so we can use the function.
% 
% < Input >
% T : [cell array] T{m,n} is a rank-4 reduced tensor, obtained by 
%       contracting a rank-5 tensor at site (m,n) from the projected
%       entangled-pair state (PEPS) as a ket state and the tensor's complex
%       conjugate. If an expectation value of a local operator or a product
%       of local operators is to be measured, then those local operators
%       are sandwiched between the physical legs of the rank-5 tensor and
%       its conjugate.
%       The legs of each T{m,n} are ordered as left-up-down-right.
% Nkeep : [numeric] The maximum bond dimension along the horizontal
%       direction.
% Nsweep : [numeric] Number of sweeps to be performed in each of MPO-MPS
%       multiplication. This parameter is handed over to "mtimes_MPO" used
%       as a sub-function here.
%
% < Output >
% res : [numeric] Contraction result.
%
% Written by S.Lee (Nov.13,2022): for the lecture course at SNU.

tobj = tic2;

% sanity check
for it1 = (1:size(T,1))
    for it2 = (1:size(T,2))
        if (it1 < size(T,1)) && (size(T{it1,it2},3) ~= size(T{it1+1,it2},2))
            error(['ERR: The down leg of T{',sprintf('%i,%i',it1,it2), ...
                '} and the up leg of T{',sprintf('%i,%i',it1+1,it2),'} do not match.'])
        elseif (it2 < size(T,2)) && (size(T{it1,it2},4) ~= size(T{it1,it2+1},1))
            error(['ERR: The right leg of T{',sprintf('%i,%i',it1,it2), ...
                '} and the left leg of T{',sprintf('%i,%i',it1,it2+1),'} do not match.'])
        elseif (it1 == 1) && (size(T{it1,it2},2) ~= 1)
            error(['ERR: The up leg of T{',sprintf('%i,%i',it1,it2),'} has non-singleton dimension.']);
        elseif (it1 == size(T,1)) && (size(T{it1,it2},3) ~= 1)
            error(['ERR: The down leg of T{',sprintf('%i,%i',it1,it2),'} has non-singleton dimension.']);
        elseif (it2 == 1) && (size(T{it1,it2},1) ~= 1)
            error(['ERR: The left leg of T{',sprintf('%i,%i',it1,it2),'} has non-singleton dimension.']);
        elseif (it2 == size(T,2)) && (size(T{it1,it2},4) ~= 1)
            error(['ERR: The right leg of T{',sprintf('%i,%i',it1,it2),'} has non-singleton dimension.']);
        end
    end
end
% % %

disptime(['Contract tensors on a ',sprintf('%i x %i',[size(T,1) size(T,2)]), ...
    ' lattice, with Nkeep = ',sprintf('%i',Nkeep)]);

% % permute legs of T{..}, to use DMRG/mtimes_MPO
for itN = (1:numel(T))
    T{itN} = permute(T{itN},[3 2 1 4]); % down-up-left-right
end

% reduced tensors are not strictly normalized, so their contraction can
% lead to very small or large numbers. In such cases, numerical procedures
% can be unstable; sometimes one gets just 0 or Inf. To avoid this,
% separate the norm of the contracted MPO after each step, and collect the
% norms as the sum of their logarithms.
logNorm = 0;

% first row of reduced tensors as an MPO
T2 = T(1,:);

% % % % TODO (start) % % % %
for it1 = (2:size(T,1)-1)
    T2 = mtimes_MPO (T(it1,:),T2,Nkeep,Nsweep);

    % factor out the norm by bringing into a right-canonical form
    % First, convert rank-4 tensors into rank-3, by merging physical legs, to
    % use the canonForm function that canonicalize MPSs
    Aloc = cell(1,size(T,2)); % isometries for merging the bottom and top legs of MPO tensors
    for it2 = (1:size(T,2))
        Aloc{it2} = getIdentity(T2{it2},1,T2{it2},2);
        T2{it2} = contract(T2{it2},4,[1 2],Aloc{it2},3,[1 2]);
    end
    % Use canonForm for MPS
    [T2,S] = canonForm(T2,0,Nkeep,[]);
    % Bring back to rank-4
    for it2 = (1:size(T,2))
        T2{it2} = contract(T2{it2},3,3,conj(Aloc{it2}),3,3,[3 4 1 2]);
    end

    logNorm = logNorm + log(S);
end

% contract the contraction result of the rows so far (all but except the
% last row) and the last row
res = 1;
for it2 = (1:size(T,2))
    Ttmp = contract(res,2,2,T2{it2},4,3);
    res = contract(T{end,it2},4,[3 2 1],Ttmp,4,(1:3));
end
% % % % TODO (end) % % % %

% restore the separted norm
res = res*exp(logNorm);

toc2(tobj,'-v');

end