function [odisc,Adisc] = getAdisc (Inrg,O1,O2,Z,varargin)
% < Description >
%
% [odisc,Adisc] = getAdisc (Inrg,O1,O2,Z  [, option]); % for fermionic operators O1 & O2
% [odisc,Adisc] = getAdisc (Inrg,O1,O2,[],[, option]); % for bosonic operators O1 & O2
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
% Inrg : [struct] NRG information obtained after running NRG_IterDiag and
%       then running getRhoFDM. 
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
logemin = loge(1);
nloge = numel(loge);
N = numel(Inrg.E0);

% cell array of operators, for convenience
Os = {O1,O2};

% Rerouting Z string: by doing this, we don't need to contract Z tensors
% when we update operators iteratively. Refer to the appendix of [A.
% Weichselbaum, Phys. Rev. B 86, 245124 (2012)]. (Especially Figs. 10 and
% 11)
if ~isempty(Z)
    for ito = (1:numel(Os))
        Os{ito} = contract(Z,2,2,Os{ito},3,1);
    end
end

for itN = (1:N)
    Atmp = {Inrg.AD{itN},Inrg.AK{itN}};
    Rtmp = {diag(Inrg.RD{itN}),Inrg.RK{itN}};
    Etmp = {Inrg.ED{itN},Inrg.EK{itN}};
    
    for it1 = (1:2) % D(iscarded), K(ept)
        for it2 = (1:2) % D, K
            if ~isempty(Atmp{it1}) && ~isempty(Atmp{it2})
                Onow = cell(size(Os));
                if itN == 1
                    for ito = (1:numel(Os))
                        Onow{ito} = updateLeft([],[],Atmp{it1},Os{ito},3,Atmp{it2});
                    end
                else
                    for ito = (1:numel(Os))
                        Onow{ito} = updateLeft(Os{ito},3,Atmp{it1},[],[],Atmp{it2});
                    end
                end

                if (it1 == 2) && (it2 == 2)
                    Os = Onow;
                else
                    % compute discrete spectral weights and add them to the
                    % result histogram 'Adisc'
                    Adisc = Adisc + ...
                        getAdisc_1shell(contract(Rtmp{it1},2,2,Onow{1},3,1),Onow{2}, ...
                        Etmp{it1},Etmp{it2},Inrg.EScale(itN), ...
                        logemin,nloge,estep);

                    Adisc = Adisc + opsign*...
                        getAdisc_1shell(Onow{1},contract(Onow{2},3,2,Rtmp{it2},2,1,[1 3 2]), ...
                        Etmp{it1},Etmp{it2},Inrg.EScale(itN), ...
                        logemin,nloge,estep);
                end
            end
        end
    end
    
    disptime(['#',sprintf('%02i/%02i',[itN-1, N-1]),' : sum(Adisc) = ', ...
        sprintf('%.4g',sum(Adisc(:)))]); % show accumulative sum of the discrete weights
end
% % % % TODO (end) % % % %

toc2(tobj,'-v');
chkmem;

end


function Adisc = getAdisc_1shell (O1,O2,E1,E2,EScale,logemin,nloge,estep)
% Compute spectral weights for given operators O1 and O2 in a given
% subspace configuration.

E21 = E2(:).'-E1(:);
M = sum(O1.*conj(O2),3);
E21 = E21(:);
M = M(:);

% indexing for frequency
ids = round((log10(abs(E21))+log10(EScale)-logemin)*estep)+1;
ids(ids<0) = 0;
ids(ids>nloge) = nloge;

% indexing for sign
oks = (E21 >= 0);
ids(oks) = ids(oks) + (nloge+1);
ids(~oks) = (nloge+1) - ids(~oks);

Adisc = accumarray(ids,M,[2*nloge+1,1]);

end