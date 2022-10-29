function [ocont,RhoV2s,iscvg] = DMFT_Hubbard_Ex (U,mu,T,D,ozin,RhoV2in,Lambda,N,Nkeep)
% < Description >
%
% [ocont,RhoV2s] = DMFT_Hubbard_Ex (U,mu,T,D,ozin,RhoV2in,Lambda,N,Nkeep)
%
% Perform the DMFT calculation of the one-band Hubbard model with the NRG
% as an impurity solver. Here we consider a lattice (on which the Hubbard
% model is defined) whose non-interacting density of states is semi-ellipse
% of half-bandwidth D. Such density of states is known to come from the
% Bethe lattice with infinite coordination number.
%
% < Input >
% U, mu, T, D : [numeric] Parameters that define the Hubbard model.
%       Interaction strength, chemical potential, temperature, and the
%       half-bandwidth of the non-interacting density of states,
%       respectively.
% ozin, RhoV2in : [vector] Initial hybridization function to start the
%       DMFT self-consistency loop. The function is assumed to be the
%       linear interpolation among the data points (ozin vs. RhoV2in), and
%       it is re-evaluted on the frequency grid "ocont" (see below); then
%       the re-evaluated data is fed into the logarithmic discretization
%       routine ("doCLD").
% Lambda, N, Nkeep : [numeric] NRG parameters. Logarithmic discretization
%       parameter, the maximum of the Wilson chain length, and the maximum
%       number of kept states, respectively. "Lambda" and "N" are directly
%       fed into "doCLD" and "Nkeep" into "NRG_IterDiag".
%
% < Output >
% ocont : [vector] Frequency grid on which the hybridization functions
%       ("RhoV2s"; see below) are defined. Here we use the default output
%       of "getAcont".
% RhoV2s : [matrix] Hybridization functions that are updated during the
%       DMFT self-consistency loop. At the n-th iteration of the loop (n =
%       1, 2, 3, ...), RhoV2s(:,n) is used as the input (fed into "doCLD").
%       The new hybridization function, as the result of the iteration by
%       using "getAdisc", "getAcont", "SEtrick", etc., is stored as
%       RhoV2s(:,n+1), which is directly used as the input for the next
%       iteration. The size of "RhoV2s" is numel(ocont)-by-(m+1), where m
%       is the number of DMFT iterations taken.
% iscvg : [logical] If true, the DMFT self-consistency is achieved up to
%       a tolerance (set by "cvgth"; see inside the code). If false, the
%       self-consistency loop is terminated since it reached the maximnum
%       iteration index (set by "Ndmft"; see inside the code).
%
% Rewritten by S.Lee (Oct.28,2022): Revised for the TN lecture course at
%       SNU.

% % default values of parameters
Ndmft = 30; % maximum number of DMFT iterations

% DMFT convergence criterion. If the difference between the old
% hybridization function (the input to "doZLD") and the new hybridization 
% function (the output from a DMFT iteration) is smaller than "cvgth" for
% all frequencies, the DMFT self-consistency is achieved.
cvgth = 1e-3; 
% % % % %

% % % % TODO (start) % % % %


ocont = getAcont(0,0,0,0);
RhoV2s = zeros(numel(ocont),Ndmft+1);
iscvg = false;

for itd = (1:Ndmft)
    [ff,gg] = doCLD(ocont,RhoV2s(:,itd),Lambda,N);

    % iterative diagonalization
    Inrg = NRG_IterDiag(H0,A0,Lambda,ff,F,gg,sum(NF,3),Z,Nkeep);

    Inrg = getRhoFDM(Inrg,T);

    [odisc,Adisc1] = getAdisc(Inrg,F(:,:,1),F(:,:,1),Z);
    [~    ,Adisc2] = getAdisc(Inrg,FHU(:,:,1),F(:,:,1),Z);
    
    [ocont,Acont1] = getAcont(odisc,Adisc1,log(Lambda),T/5);
    [~    ,Acont2] = getAcont(odisc,Adisc2,log(Lambda),T/5);

    SE = SEtrick(ocont,Acont1,Acont2);

    xi = ocont + mu - SE;
    Glat = (2/(D^2))*(xi - 1i*sqrt(D^2 - xi.^2));
    Alat = (-1/pi)*imag(Glat);
    RhoV2s(:,itd+1) = ((D/2)^2)*Alat;

    if max(abs(RhoV2s(:,itd+1)-RhoV2s(:,itd)),[],1) < cvgth
        disptime('DMFT converged; exit the loop.');
        iscvg = true;
        break;
    end
end

RhoV2s(:,itd+2:end) = [];

% % % % TODO (end) % % % %


end