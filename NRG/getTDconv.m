function [T,Tchi,Sent] = getTDconv (Inrg,Sz,beta0,varargin)
% < Description >
%
% [T,Tchi,Sent] = getTDconv (Inrg,Sz,beta0 [, S_A0L])
%
% Compute thermodynamic properties T*chi (= temperature * static spin
% susceptibility) and entropy by using the conventional NRG method (see
% Sec. III A 1. in Bulla2008 [R. Bulla et al., Rev. Mod. Phys. 80, 395
% (2008)]).
% Here each set of energy levels in the iteration 'n' (so-called "shell")
% is considered as the restricted but effective set of energy eigenstates
% responsible for temperature T = Inrg.EScale/beta0. Then we can compute
% the expectation value < .. > of operators (e.g., spin-z operator S_z)
% and the partition function Z as the sum of the Boltzmann factors, for the
% spin susceptibility chi multiplied by T is given by
%
%   T*chi = < S_z^2 > - < S_z >^2 
%     (see the equation just before Eq. (47) in Bulla2008)
%
% and the entropy is given by
%
%   S = < H > / T + ln (Z)
%     (see Eq. (50) in Bulla2008)
%
% These properties are for the whole system. So the impurity contribution
% to the spin susceptiblity and the entropy can be obtained by substracting
% the properties of only the bath (by running the NRG for only the bath in
% the absence the impurity) from those of the full system (by runing the
% NRG for the target impurity model).
%
% < Input >
% Inrg : [struct] Struct containing NRG result, as the output from function
%        NRG_IterDiag_Ex.
% Sz : [numeric] Spin-z operator acting onto the physical space associated
%       with the third leg of A0 (the input isometry to NRG_IterDiag_Ex)
%       and the spaces associated with the bath sites that are added along
%       the iterative diagonalization.
% beta0 : [number] Prefactor to define the effective tempearture T for all
%       the iterations which is an output as well; see below for detail.
%       The values of temperature associated with the iterations are given
%       by:
%       T = [Inrg.EScale(2)*sqrt(Inrg.Lambda), Inrg.EScale(2:end)]/beta0;
%       Here Inrg.EScale is the energy scale of the iterations which is
%       used to rescale the energy eigenvalues along the iterative
%       diagonalization. And Inrg.Lambda is the logarithmic discretization
%       parameter Lambda. By the convention of NRG_IterDiag_Ex,
%       Inrg.EScale(1) is always 1 so that one can easily read off the
%       energy eigenvalues of the impurity without considering prefactor.
%       Therefore, to have a consistent logarithmic scaling along T, the
%       first element T(1) is determined by multiplying sqrt(Inrg.Lambda)
%       to Inrg.EScale(2).
% S_A0L : [numeric] (Optional) Spin-z operator acting onto the space
%       associated with the left leg (i.e., first leg) of A0. If the
%       left-leg space of A0 is not vacuum (i.e., the associated dimension
%       is larger than 1), the input Sz is regarded as the spin operator
%       acting (only) onto the left-leg space. If the left-leg space is
%       vacuum (i.e., the associated dimension is 1), the spin operator is
%       not considered.
%
% < Output >
% T : [vector] Temperature values corresponding to shell indices. It is in
%       descending order.
% Tchi : [vector] Temperature * static spin susceptibility. In the Kondo
%       regime, it should decay linearly with T.
% Sent : [vector] Entropy.
%
% Written by S.Lee (May 05,2017); edited by S.Lee (Aug.09,2017)
% Updated by S.Lee (May 12,2019): Revised for SoSe 2019.
% Updated by S.Lee (Jun.20,2020): Revised for SoSe 2020.
% Updated by S.Lee (Oct.15,2022): Revised for the semester at SNU.

tobj = tic2;
disptime('Compute thermodynamic properties with the conventional NRG method.');

% parsing optional input
S_A0L = [];
if ~isempty(varargin)
    S_A0L = varargin{1};
end

% % % sanity check of input
if (size(Sz,1) ~= size(Sz,2)) || (ndims(Sz) > 3)
    error('ERR: ''Sz'' should be a square matrix.');
elseif (size(Inrg.AK{1},1) > 1) && isempty(S_A0L)
    error('ERR: If the first leg (= left leg) of A0 is non-singleton, ''S_A0L'' should be set.');
elseif ~isempty(S_A0L) && ((ndims(S_A0L) > 3) || ~all(size(Inrg.AK{1},1) == size(S_A0L,1)))
    error('ERR: ''S_A0L'' should be a square matrix whose dimensions match with the first leg (= left leg) of A0.');
end
for itN = (1:numel(Inrg.AK))
    if (~isempty(Inrg.AK{itN}) && (size(Sz,2) ~= size(Inrg.AK{itN},3))) || ...
        (~isempty(Inrg.AD{itN}) && (size(Sz,2) ~= size(Inrg.AD{itN},3)))
        error('ERR: Dimensions of ''Sz'' are not consistent with the third-leg dimension of isometries.');
    end
end
% % % % % 

N = numel(Inrg.EScale); % total number of iterations
T = [Inrg.EScale(2)*sqrt(Inrg.Lambda), Inrg.EScale(2:end)]/beta0; % temperature

% properties before even-odd averaging
Tchi0 = zeros(1,N); % T*(spin susceptibility)
Sent0 = zeros(1,N); % entropy

% % % % TODO (start) % % % %
for itN = (1:N)
    if itN == 1
        % % For the bottom-leg space of A0
        % Sz operator in kept state basis
        SzK = updateLeft([],[],Inrg.AK{itN},Sz,2,Inrg.AK{itN});
        % Sz operator in discarded state basis
        SzD = updateLeft([],[],Inrg.AD{itN},Sz,2,Inrg.AD{itN});
        % NOTE: here we consider all the states in a single iteration; both
        % kept and discarded states. By doing this, we can cover a wider
        % range of energy levels to reduce the numerical artifact
        % originating from finite energy threshold.
        
        % % For the left-leg space of A0:
        % If not vacuum, add the spin operator acting on the space
        if size(Inrg.AK{itN},1) > 1
            SzK = SzK + updateLeft(S_A0L,3,Inrg.AK{itN},[],[],Inrg.AK{itN});
        end
        if size(Inrg.AD{itN},1) > 1
            SzD = SzD + updateLeft(S_A0L,3,Inrg.AD{itN},[],[],Inrg.AD{itN});
        end
    else
        % Obtain the Sz operators in kept/discarded state basis. Use
        % SzKprev which is the Sz operator in the kept state basis at the
        % last iteration (see below for its definition), and the Sz
        % operator at the current iteration.
        
        % SzK: Sz operator in kept state basis
        SzK = updateLeft(SzKprev,2,Inrg.AK{itN},[],[],Inrg.AK{itN}) + ...
            updateLeft([],[],Inrg.AK{itN},Sz,2,Inrg.AK{itN});
        % SzD: Sz operator in discarded state basis
        SzD = updateLeft(SzKprev,2,Inrg.AD{itN},[],[],Inrg.AD{itN}) + ...
            updateLeft([],[],Inrg.AD{itN},Sz,2,Inrg.AD{itN});
    end
    
    % matrix of the Sz operator for the kept+discarded space
    Sztot = blkdiag(SzK,SzD);
    % third dimensions of SzK and SzD are singleton
    E = [Inrg.EK{itN};Inrg.ED{itN}]; % rescaled energy eigenvalues within a shell
    
    % Obtain the temperature * spin susceptibility Tchi0(itN) and the
    % entropy Sent0(itN) for the current iteration
    % NOTE: be aware of that the E above is *rescaled* energy values.
    
    Zfac = exp(-E*(Inrg.EScale(itN)/T(itN))); % Boltzmann weights
    Zsum = sum(Zfac); % partiton function
    rho = Zfac(:)/Zsum; % diagonal of density matrix

    Tchi0(itN) = (rho.'*diag(Sztot*Sztot)) - (rho.'*diag(Sztot))^2;
    Sent0(itN) = (rho.'*E)*(Inrg.EScale(itN)/T(itN)) + log(Zsum);
    
    SzKprev = SzK;
end
% % % % TODO (end) % % % %

% % even-odd averaging
Tchi = (interp1(T(1:2:end),Tchi0(1:2:end),T,'linear','extrap') + ...
    interp1(T(2:2:end),Tchi0(2:2:end),T,'linear','extrap'))/2;
Sent = (interp1(T(1:2:end),Sent0(1:2:end),T,'linear','extrap') + ...
    interp1(T(2:2:end),Sent0(2:2:end),T,'linear','extrap'))/2;

toc2(tobj,'-v');

end