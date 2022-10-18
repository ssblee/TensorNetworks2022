function [ff,gg] = doCLD (ozin,RhoV2in,Lambda,N,varargin)
% < Description >
%
% [ff,gg] = doCLD (ozin,RhoV2in,Lambda,N [, option] )
%
% This function performs the Campo--Oliveira logarithmic discretization [V.
% L. Campo and L. N. Oliveira, Phys. Rev. B 72, 104432 (2005)] of the input
% hybridization function (paramterized by the inputs "ozin" and "RhoV2in")
% and maps the resulting star-geometry Hamiltonian onto the chain-geometry
% Hamiltonian, so-called the Wilson chain. The output "ff" and "gg"
% describe the hopping amplitudes and the on-site energies along the chain.
%
% < Input >
% ozin, RhoV2in : [numeric vector] RhoV2in(n) is the value of the spectral
%       function part of the hybiridzation function (-1/pi)*imag[ \Delta(
%       \omega)] evaluated at \omega = ozin(n), where \Delta(\omega) is the
%       hybridization function.
% N : [integer] Wilson chain length.
% Lambda : [numeric] Discretization parameter.
%
% < Option >
% 'estep', ... : [numeric] "estep", "emax", and "emin" are the paramters
%       that define the logrithmic frequency grid, say "oz", to be used in
%       the discretization (which can be different from "ozin"). "estep" is
%       the number of steps between a frequency grid point, say oz(n), and
%       oz(n)*Lambda.
%       (Default: 40)
% 'emax', ... : [numeric] Maximum frequency value of the grid "oz".
%       (Default: max(abs(ozin)))
% 'emin', ... : [numeric] Maximum frequency value of the grid "oz".
%       (Default: 1e-20)
% 'minrho', ... : [numeric] Parameter to set effective band edges when the
%       input "RhoV2in" has long tails. The discretization grid points are
%       defined, starting those band edges, dividing \Lambda's.
%       In this "doCLD" function, the hybridization for positive
%       frequencies and that for negative frequencies are discretized 
%       separately. For example, on the positive side, the largest (in
%       terms of absolute value) frequency grid point, say ozin(m), at
%       which RhoV2in(m) >= rhomin*max(RhoV2in(ozin>0)) is chosen as the
%       effective band edge. Similar for the negative side.
%       (Default: 0.01)
% 'fftol', ... : [numeric] Tolerance for the hopping amplitudes "ff"
%       within the Lanczos tridiagonalization. If ff(n) < Lambda^(-n/2)*
%       fftol, then the tridiagonalization stops and outputs only ff(1:n-1)
%       and gg(1:n-1).
%       (Default: 1e-2)
% 
% < Output >
% ff, gg : [numeric vectors] Hopping amplitudes and on-site energies of
%       the Wilson chain, respectively. The hopping amplitudes correspond
%       to the superdiagonals [diag(..,+1) and diag(..,-1)] of the
%       tridiagonal matrix representation of a single-particle Hamiltonian;
%       the on-site energies correspond to the diagonals of the tridiagonal
%       matrix.
%
% Written by S.Lee (May 5,2017); edited by S.Lee (May 9,2017)
% Updated by S.Lee (Jun.15,2020): Revised for SoSe 2020.
% Updated by S.Lee (Oct.15,2022): Revised for the semester at SNU.

% default parameter
estep = 10;
emax = max(abs(ozin));
emin = 1e-20;
minrho = 0.01;
fftol = 0.01;

while ~isempty(varargin)
    switch varargin{1}
        case 'estep'
            estep = varargin{2};
            varargin(1:2) = [];
        case 'emax'
            emax = varargin{2};
            varargin(1:2) = [];
        case 'emin'
            emin = varargin{2};
            varargin(1:2) = [];
        case 'minrho'
            minrho = varargin{2};
            varargin(1:2) = [];
        case 'fftol'
            fftol = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: check input!');
    end
end

% parsing input
ozin = ozin(:);
RhoV2in = RhoV2in(:);

if isempty(ozin)
    error('ERR: Empty frequency input (1st input)');
elseif isempty(RhoV2in)
    error('ERR: Empty hybridization input (2nd input)');
elseif numel(ozin) ~= numel(RhoV2in)
    error('ERR: Different # of elements between frequency and hybridization inputs (1st & 2nd)');
end

% logarithmic frequency grid on which the hybridization function to be discretized is defined
xs = (ceil(log(emax)/log(Lambda)*estep):-1:floor(log(emin)/log(Lambda)*estep))/estep;
xs = flipud(xs(:)); % increasing, column
oz = Lambda.^xs; % increasing

rho1 = interp1(ozin,RhoV2in,+oz,'linear','extrap');
rho1(rho1<0) = 0; % to ensure non-negativity
[repE1,repT1] = doCLD_1side (oz,rho1,estep,minrho);

rho2 = interp1(ozin,RhoV2in,-oz,'linear','extrap');
rho2(rho2<0) = 0; % to ensure non-negativity
[repE2,repT2] = doCLD_1side (oz,rho2,estep,minrho);

if (numel(repE1)+numel(repE2)) < N
    fprintf(['WRN: Number of discretization intervals (= ', ...
        sprintf('%i',numel(repE1)+numel(repE2)),') is smaller than the chain length N (= ', ...
        sprintf('%i',N),'\n']);
    N2 = numel(repE1) + numel(repE2);
else
    N2 = N;
end

ff = zeros(N2,1); % hopping amplutudes; corresponds to the super diagonal 
                  % (diag(..,+1) or diag(..,-1)) in the tridiagonal matrix
gg = zeros(N2,1); % hopping amplutudes; corresponds to the elements (2:end) 
                  % of the main diagonal (diag(..)) in the tridiagonal matrix

% % % % TODO (start) % % % %
% % Lanczos tridiagonalization
% star-geometry Hamiltonian
Xis = [repE1; -repE2];
Gammas = [sqrt(repT1); sqrt(repT2)];
H = [0 Gammas'; Gammas, diag(Xis)];

U = zeros(size(H,1)*[1 1]);
U(1,1) = 1;

for itN = (1:N2)
    v = H*U(:,itN);

    % Orthgonalize
    v = v-U(:,1:itN)*(U(:,1:itN)'*v);
    v = v-U(:,1:itN)*(U(:,1:itN)'*v); % twice for numerical reason

    ff(itN) = norm(v);

    if ff(itN) < (Lambda^(-itN/2)*fftol)
        itN = itN-1;
        break;
    end

    U(:,itN+1) = v/ff(itN);
    gg(itN) = U(:,itN+1)'*H*U(:,itN+1);
end

ff(itN+1:end) = [];
gg(itN+1:end) = [];

% % % % TODO (end) % % % %

end


function [repE,repT] = doCLD_1side (oz,rho,estep,minrho)
% Obtain the representative energies (repE, \mathcal{E} in Campo2005) and
% the integral of the hybridization function (repT) for each discretization
% interval, for either positive or negative energy side.

ids0 = find(rho >= minrho*max(rho),1,'last')-estep;

% oz(ids) are the discretization grid points at which the input
% hybridization function is splot
ids = [numel(oz),(ids0:-estep:1)];

% % % % TODO (start) % % % %
repT = zeros(numel(ids)-1,1);
repE = zeros(size(repT));

for itx = (1:numel(repT))
    ozp = oz(ids(itx+1):ids(itx));
    rhop = rho(ids(itx+1):ids(itx));

    % perform numerical integration
    repT(itx) = sum((rhop(2:end)+rhop(1:end-1)).*(ozp(2:end)-ozp(1:end-1)))/2;
    repE(itx) = (rhop(end)-rhop(1)) + ...
            sum( (ozp(2:end).*rhop(1:end-1) - ozp(1:end-1).*rhop(2:end)) ./ ...
                (ozp(2:end) - ozp(1:end-1)) .* log(abs(ozp(2:end)./ozp(1:end-1))) );
end

repE = repT./repE;
% % % % TODO (end) % % % %

end