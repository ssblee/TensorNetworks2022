function [odisc,Adisc] = getAdisc_Ex (Inrg,O1,O2,Z,varargin)
% < Description >
%
% [odisc,Adisc] = getAdisc_Ex (Inrg,O1,O2,Z  [, option]); % for fermionic operators O1 & O2
% [odisc,Adisc] = getAdisc_Ex (Inrg,O1,O2,[],[, option]); % for bosonic operators O1 & O2
%
% Obtain the discrete spectral function of the correlation function, by
% evaluating its Lehmann representation via the full density matrix NRG
% (fdmNRG). The correlation function is defined by a pair of operators O1
% and O2,
% 
%          G = -1i*\theta(t) Tr [ \rho * [O1(t), O2']_{+-} ],
%
% where \theta(t) means the theta function, Tr means the trace, and O1(t)
% is a time evolution of O1 in the Heisenberg picture. The commutator takes
% + (-) if O1 and O2 are fermionic (bosonic). Note that O2 is Hermitian-
% conjugated in the definition of the correlator we use.
%
% In this fdmNRG routine, we assume that both O1 and O2 act on the third (=
% bottom) leg of A0 (the isometry for the impurity site, also an input to
% NRG_IterDiag). The resulting delta functions from the Lehmann
% representation are binned into a histogram vector Adisc.
%
% < Input >
% Inrg : [struct] NRG information obtained after running NRG_IterDiag.
% O1, O2 : [tensor] Operator acting at the impurity, corresponding to the
%       third (= bottom) leg of A0.
% Z : [tensor] Fermionic sign operator at each chain site. If O1 and O2 are
%       bosonic, then put Z as empty []. 
%
% < Option >
% 'emin', .. : [number] Minimum (in absolute value) frequency limit.
%               (Default: 1e-12)
% 'emax', .. : [number] Minimum (in absolute value) frequency limit.
%               (Default: 1e4)
% 'estep', .. : [number] Number of bins per decade, i.e., number of bins
%               from frequency 1 to 10.
%               (Default: 500)
%
% < Output >
% odisc : [vector] Logarithmic grid of frequency.
% Adisc : [vector] Discrete spectral weights corresponding to odisc.
%       Adisc(n) is the sum of the spectral weights whose frequency
%       positions are close to the bin center position odisc(n).
%
% Written by S.Lee (May 22,2017)
% Updated by S.Lee (May 11,2019): Revised for SoSe 2019.
% Updated by S.Lee (Jun.20,2020): Revised for SoSe 2020.
% Updated by S.Lee (Oct.25,2020): Revised for the 2022 Fall semester at SNU.

tobj = tic2;

% default option values
emin = 1e-12;
emax = 1e4;
estep = 500;

% % parsing option
while ~isempty(varargin)
    switch varargin{1}
        case 'emin'
            emin = varargin{2};
            varargin(1:2) = [];
        case 'emax'
            emax = varargin{2};
            varargin(1:2) = [];
        case 'estep'
            estep = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: Unknown option.');
    end
end
% % % 

% sanity check
if ~all([size(Inrg.AK{1},1),size(Inrg.AD{1},1)] == 1)
    error('ERR: The first (= left) leg of the first isometries (for the impurity site) should have singleton dimensions.')
end
if ~isempty(Z)
    disptime('Correlation function for anti-commuting operators.');
    opsign = 1; % anti-commutation
else
    disptime('Correlation function for commuting operators');
    opsign = -1; % commutation
end

loge = (floor(log10(emin)*estep):ceil(log10(emax)*estep)).'/estep;
odisc = 10.^(loge);
odisc = [-flipud(odisc);0;odisc]; % bin centers; zero frequency in the middle
Adisc = zeros(size(odisc)); % to be binned


% % % % TODO (start) % % % %


for itN = (1:N)
    
    
    disptime(['#',sprintf('%02i/%02i',[itN-1, N-1]),' : sum(Adisc) = ', ...
        sprintf('%.4g',sum(Adisc(:)))]); % show accumulative sum of the discrete weights
end
% % % % TODO (end) % % % %

toc2(tobj,'-v');
chkmem;

end