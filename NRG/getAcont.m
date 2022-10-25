function [ocont,Acont] = getAcont (odisc,Adisc,sigmab,g,varargin)
% < Description >
%
% [ocont,Acont] = getAcont (odisc,Adisc,sigmab,g [,option])
%
% Broadening discrete spectral data of correlation functions given by
% getAdisc. Here we apply two broadenings (i.e., convolutions) over all
% frequencies sequentially.
% The primary broadening is logarithmic broadening. It uses logarithmic
% Gaussian or Gaussian kernel, where the broadening width (in linear
% frequency scale) is proportional to the frequency bin position. For a
% discrete weight at frequency bin of w', the value of kernel at w is given
% by the symmetric logarithmic Gaussian,
%
%   L(w,w') = exp(-(log(w'/w)/sigmab - sigmab/4)^2)/(sqrt(pi)*sigmab*abs(w')),
%                 for sign(w) == sign(w')),
%   L(w,w') = 0,  otherwise.
%
% Here sigmab is an input parameter.
% Then, appply the secondary linear broadening of width g to remove the
% artifacts at low frequency. For the secondary broadening kernel, we use
% the derivative of Fermi-Dirac function,
%
%   f(w,w') = 1/(1+cosh((w - w')/g))*(1/2/g)
%           = - d(Fermi-Dirac function of temperature g)/d(energy)
%
% For detail, refer to [S.-S. B. Lee and A. Weichselbaum, Phys. Rev. B 94,
% 235127 (2016)].
%
% < Input >
% odisc : [numeric vector] Logarithimic frequency bins. Here the original
%         frequency values from the differences b/w energy eigenvalues are
%         shifted to the closest bins.
% Adisc : [numeric vector] Spectral function.
% sigmab : [numeric vector] The width parameter for the primary broadening
%       kernel.
% g : [numeric] The width parameter for secondary linear broadening kernel.
%
% < Option >
% 'emin', .. : [numeric] Minimum absolute value of frequency grid. Set this
%       be equal to or smaller than the minimum of finite elements of
%       odisc, to properly broaden the spectral weights at frequencies
%       lower than 'emin'. The spectral weights binned at frequencies
%       smaller that emin are *not* broadened by the primary logarithmic
%       broadening; they are broadened only by the secondary linear
%       broadening.
%       (Default : 1e-12)
% 'emax', .. : [numeric] Maximum absoulte value of frequency grid.
%       (Default : 1e4)
% 'estep', .. : [integer] Number of frequency grid points per decade, i.e.,
%       between a frequency and the frequency times 10 (e.g., between 1 and
%       10).
%       (Default: 200)
% 'tol', .. : [numeric] Minimum value of (spectral weight)*(value of
%       broadening kernel) to consider. The spectral weights whose
%       contribution to the curve Acont are smaller than tol are not
%       considered. Also the far-away tails of the broadening kernel, whose
%       resulting contribution to the curve are smaller than tol, are not
%       considered also.
%       (Default : 1e-14)
% '-v' : Show details.
%
% < Output >
% ocont : [numeric vector] Logarithimic frequency grid.
% Acont : [numeric vector] Smoothened spectral function, on the frequency
%       grid 'ocont'.
%
% Written by S.Lee (May 22,2017)
% Rewritten by S.Lee (Jun.20,2020)

% default option values
emin = 1e-12;
emax = 1e4;
estep = 200;
tol = 1e-14;

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
        case 'tol'
            tol = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: Unknown option.');
    end
end
% % % 

odisc = real(odisc(:)); % column vector

if isempty(Adisc)
    error('ERR: ''Adisc'' is empty');
elseif isempty(odisc)
    error('ERR: ''odisc'' is empty');
elseif isempty(sigmab)
    error('ERR: ''sigmab'' is empty');
elseif estep ~= round(estep)
    error('ERR: ''nstep'' needs to be integer');
elseif numel(odisc) ~= size(Adisc,1)
    error('ERR: Input dimensions do not match; numel(odisc) ~= size(Adisc,1)');
elseif numel(g) > 1
    error('ERR: numel(g) > 1');
end

if ~(emin > 0)
    error('ERR: ''emin'' should be positive and finite.');
elseif ~(emax > emin)
    error('ERR: ''emax'' should be larger than''emin''.');
elseif tol <= 0
    error('ERR: ''tol'' should be positive and finite.');
elseif estep < 1
    error('ERR: ''estep'' should be equal to or larger than 1.');
end

% temporary frequency grid; includes 0 in the middle
xs = (floor(log10(emin)*estep):ceil(log10(emax)*estep)).'/estep; % log-frequencies; increasing, column vector
ocont = 10.^xs;
ocont = [flipud(-ocont);0;ocont]; % center of frequency bins
docs = [ocont(2)-ocont(1);(ocont(3:end)-ocont(1:end-2))/2;ocont(end)-ocont(end-1)]; % width of frequency bins
Acont = zeros(numel(ocont),1); % temporary result

oks1 = (odisc >= emin);
if any(oks1)
    odtmp = odisc(oks1);
    Adtmp = Adisc(oks1,:);
    [ots,dots,yts] = getAcont_logBroaden(odtmp,Adtmp,sigmab,tol,emin,emax,estep);
    Acont = getAcont_linBroaden(ots,dots,yts,g,ocont,docs,Acont,tol);
end

oks2 = (odisc <= -emin);
if any(oks2)
    odtmp = -odisc(oks2); % negative -> positive
    Adtmp = Adisc(oks2,:);
    [ots,dots,yts] = getAcont_logBroaden(odtmp,Adtmp,sigmab,tol,emin,emax,estep);
    Acont = getAcont_linBroaden(-ots,dots,yts,g,ocont,docs,Acont,tol); % -ots: return to negative frequency
end

oks3 = ~any([oks1,oks2],2);
Acont = getAcont_linBroaden(0,docs((end+1)/2),sum(sum(Adisc(oks3,:))),g,ocont,docs,Acont,tol);

end


function [ots,dots,yts] = getAcont_logBroaden (odisc,Adisc,sigmab,tol,emin,emax,estep)
% primary broadening: log-Gaussian broadening

% Determine the additional range of frequencies smaller than emin. It is
% necessary since the log-Gaussian broadening of the spectral weights at
% low-frequency bins result in the tail towards small frequency, and the
% tail should be fully included in the frequency range. Otherwise, the
% spectral weight of the truncated tail is missing in the secondary linear
% broadening, which can lead to wrong result.
Atmp = (tol*sqrt(pi))*(odisc.*sigmab./abs(Adisc));
okA = (Atmp > 1);
Atmp(okA) = 1;

% widthx : the maximum half-range of the kernel in the log-frequency
widthx = abs(sqrt(-log(Atmp))).*sigmab;

xtmp = log(odisc)-widthx-((sigmab.^2)/4);

% xmin : the log of the minimum frequency of the additional range, which
%       is given by the lower edge of the broadening kernel (times spectral
%       weight), which is cut off by tol. 
xmin = min(xtmp(:))/log(10);
% modify xmin so that the grid of xts matches with the grid of the result
% ocont. The mismatch between xts and ocont can lead to an error later in
% the linear broadening step!
xdiff = (xmin - log10(emin))*estep;
xmin = xmin - (xdiff - floor(xdiff))/estep;

% % temporary frequency grid
xts = ((xmin*estep):(log10(emax)*estep)).'/estep; % exponents of frequencies (base 10); increasing, column vector

if numel(xts) > 1
    xtse = xts*log(10); % exponents of frequencies (base e)
    ots = 10.^xts; % center of frequency bin
    dots = [ots(2)-ots(1);(ots(3:end)-ots(1:end-2))/2;ots(end)-ots(end-1)]; % width of frequency bin
    yts = zeros(numel(ots),1); % temporatry result; weights are to be binned
    
    omat = odisc + zeros(1,numel(sigmab));
    smat = zeros(numel(odisc),1) + sigmab;
    
    cenx = log(odisc) - ((sigmab.^2)/4);

    omat(okA) = [];
    smat(okA) = [];
    cenx(okA) = [];
    widthx(okA) = [];
    Adisc(okA) = [];
    
    for ito = (1:numel(omat))
        okt = (abs(xtse - cenx(ito)) <= widthx(ito));
        if any(okt)
            ytmp = -( (cenx(ito)-xtse(okt))/smat(ito) ).^2;
            ytmp = exp(ytmp)*(Adisc(ito)/sqrt(pi)/smat(ito)/omat(ito));

            yts(okt) = yts(okt) + ytmp.*dots(okt); % height -> weight
        end
    end
else
    ots = []; yts = [];
end
        
end

function Acont = getAcont_linBroaden (ots,dots,yts,g,ocont,docs,Acont,tol)
% secondary linear broadening

if ~isempty(ots) && ~isempty(yts)
    Atmp = (tol*2*g)./abs(yts);
    okA = (Atmp > 0.5);
    Atmp(okA) = 0.5;
    widtho = abs(acosh((1./Atmp(~okA))-1))*g;
    
    yts(okA) = [];
    ots(okA) = [];
    dots(okA) = [];
    
    if ~isempty(yts)
        widtho = max(widtho,dots); % in case of widtho == 0

        % in case when min(abs(octmp(okt) - ots(ito))) is much larger than g
        g2 = max(g,dots/25);

        for ito = (1:numel(ots))
            okt = (abs(ocont - ots(ito)) <= widtho(ito));
            if any(okt)
                ytmp = 1./(1+cosh((ocont(okt) - ots(ito))/g2(ito))); %/(2*g);
                ytmp = ytmp.*docs(okt); % height -> weight
                ytmp = ytmp/sum(ytmp); % normalize
                Acont(okt) = Acont(okt) + (ytmp./docs(okt))*yts(ito); % weight -> height
            end
        end
    end
end

end