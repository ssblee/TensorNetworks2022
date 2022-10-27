function Inrg = getRhoFDM (Inrg,T)
% < Description >
%
% Inrg = getRhoFDM_Ex (Inrg,T)
%
% Construct the full density matrix (FDM) in the basis of both discarded
% and kept states, for given temperature T.
%
% < Input >
% Inrg : [struct] NRG information obtained after running NRG_IterDiag.
% T : [number] Temperature. Here we set \hbar = k_B = 1.
%
% < Ouput >
% Inrg : [struct] NRG result. It keeps the result of NRG_IterDiag. In
%       addition to the result, this function adds two more fields to Inrg:
%   .RD, .RK : [cell] Full density matrix in the discarded and kept state
%       basis, respectively. Each cell element Inrg.RD{n} is a column
%       vector whose elements are the density matrix elements associated
%       with the discarded energy eigenstates at iteration n-1. (Refer to
%       the documentation of NRG_IterDiag for the interation indexing
%       convention.)
%       Inrg.RK{n} is a matrix in the basis of the kept energy eigenstates
%       at the iteration n-1.
%
% Written by S.Lee (May 22,2017)
% Updated by S.Lee (May 12,2019): Revised for SoSe 2019.
% Updated by S.Lee (Jun.20,2020): Revised for SoSe 2020.
% Updated by S.Lee (Oct.22,2020): Revised for the 2022 Fall semester at SNU.

tobj = tic2;
disptime(['Construct full density matrix @ T = ',sprintf('%.4g',T),'.']);

N = numel(Inrg.E0);

% extract the local space dimension from ket tensors
locdim = zeros(N,1);
for itN = (1:N)
    if ~isempty(Inrg.AK{itN})
        locdim(itN) = size(Inrg.AK{itN},3);
    else
        locdim(itN) = size(Inrg.AD{itN},3);
    end
end

RD = cell(1,N); % FDM in the discarded state basis; row vector
RK = cell(1,N); % FDM in the kept state basis; matrix
Ztot = zeros(1,N); % sum of Boltzmann weights

% % % % TODO (start) % % % %
% the shift of energy in each shell measured from the lowest-energy of the
% last iteration
E0r = [Inrg.EScale(2:end).*Inrg.E0(2:end),0];
E0r = fliplr(cumsum(fliplr(E0r)));

% obtain the Boltzamann weights
for itN = (1:N)
    % Obtain the column vector RD{itN} whose elements are the Boltzmann
    % weights
    RD{itN} = exp(-(Inrg.ED{itN}*Inrg.EScale(itN)-E0r(itN))/T)*prod(locdim(itN+1:end));
    Ztot(itN) = sum(RD{itN});
end

Ztot = sum(Ztot);

% normalize the Boltzmann weights to get the elements of the density matrix
% in the discarded basis
for itN = (1:N)
    RD{itN} = RD{itN}/Ztot;
end

% update the FDM in the kept basis
for itN = (N:-1:2)
    % Construct RK{itN-1} as the sum of RD{itN} and RK{itN}, with the local
    % Hilbert space for the site s(itN-1).
    AD2 = permute(conj(Inrg.AD{itN}),[2 1 3]);
    AK2 = permute(conj(Inrg.AK{itN}),[2 1 3]);
    RK{itN-1} = updateLeft(diag(RD{itN}),2,AD2,[],[],AD2) ...
        + updateLeft(RK{itN},2,AK2,[],[],AK2);
    % NOTE: AK and AD are in left-canonical form, not right-canonical.

end
% % % % TODO (end) % % % %

if sum(RD{end}) > 1e-2
    disptime('WRN: sum(Inrg.RD{end}) > 1e-2 ; chain length may not be long enough.');
end

Inrg.T = T; % record T in Inrg
Inrg.Ztot = Ztot; % record the sum of Boltzmann weights
Inrg.RK = RK;
Inrg.RD = RD;

toc2(tobj,'-v');

end