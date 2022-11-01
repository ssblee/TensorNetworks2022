function [SE,varargout] = SEtrick (ocont,Acont1,Acont2,varargin)
% < Description >
%
% SE = SEtrick (ocont,Acont1,Acont2) % self-energy only
% [SE,Aimp] = SEtrick (ocont,Acont1,Acont2,epsd,RhoV2in) % also obtain
%                % the improved estimate of the impurity spectral function
%
% Compute the impurity self-energy due to interactions and (if "epsd" and
% "RhoV2in" are given) the improved estimate of the impurity spectral
% function by using Bulla's "self-energy trick", developed in [R. Bulla, A.
% C. Hewson, and Th. Pruschke, J. Phys. Condens. Matter 10, 8365 (1998)].
% Based on the equation of motion (EoM) for the impurity Green's function
% G[F, F'] that is derived by differentiating the function with respect to
% the first time argument for the first operator F, the self-energy SE can
% be given by the ratio of two correlation functions,
%       SE = G[FH, F'] / G[F, F'],                         ----- (1)
% where F is the particle annihilation operator at the impurity, F' is the 
% Hermitian conjugate of F, FH = [F, HU] is the commutator, and HU means
% the interaction (i.e., non-quadratic) term of the impurity Hamiltonian.
%
% < Input >
% ocont : [numeric vector] Frequency grid on which "Acont1" and "Acont2"
%       (also "RhoV2in", if given) are defined.
% Acont1 : [numeric] The spectral function part of a correlator G[F, F'].
%       This 'Acont1' is the broadned curve (by using 'getAcont') from the
%       discrete spectral data of the spectral function of G[F,F'] (by
%       running 'getAdisc').
% Acont2 : [numeric] The spectral function part of a correlator G[FH, F'].
%       It is obtained in a similar way as "Acont1", while using FH instead
%       of F.
% epsd : (Optional) [numeric] Impurity level.
% RhoV2in : (Optional) [numeric] The spectral function part of the
%       hybridization function defined on the frequency grid "ocont".
%
% < Output >
% The data structure (in terms of dimensions) of the following results are
% the same as the input Acont1.
% SE : [numeric] Impurity self-energy.
% Aimp : (Optional) [numeric] Improved impurity spectral function, obtained
%       by susbstituting the self-energy "SE" and the optional input
%       "RhoV2in" into the Dyson equation.
%
% Written by S.Lee (Oct.28,2022): for the TN lecture course at SNU.

% default values of numerical parameter and optional inputs
iSEmax = -1e-14; % maximum of the imaginary part of the self-energy, needed 
                 % to post-process the self-energy; see below
epsd = [];
RhoV2in = [];

% parse optional input
if ~isempty(varargin)
    if numel(varargin) < 2
        error('ERR: If optional inputs are given, there should be two optional inputs.');
    end
    epsd = varargin{1};
    RhoV2in = varargin{2};
end

% sanity check
if ~isvector(ocont)
    error('ERR: ''ocont'' should be a vector.');
elseif ~ismatrix(Acont1) || ~ismatrix(Acont2)
    error('ERR: ''Acont1'' and ''Acont2'' should be matrices.');
elseif ~all(numel(ocont) == [size(Acont1,1) size(Acont2,1)])
    error('ERR: size(Acont1,1) and size(Acont2,1) should be equal to numel(ocont).');
elseif size(Acont1,2) ~= size(Acont2,2)
    error('ERR: size(Acont1,2) and size(Acont2,2) should be the same.');
elseif isempty(RhoV2in) ~= isempty(epsd)
    error('ERR: Only one optional input was set?');
elseif ~isempty(RhoV2in)
    if ~ismatrix(RhoV2in)
        error('ERR: RhoV2in should be a matrix.');
    elseif numel(ocont) ~= size(RhoV2in,1)
        error('ERR: size(RhoV2in,1) should be equal to numel(ocont).');
    elseif size(epsd,1) ~= 1
        error('ERR: ''epsd'' should be a row vector.');
    elseif ~all(size(Acont1,2) == [size(epsd,2) size(RhoV2in,2)])
        error('ERR: size(epsd,2) and size(RhoV2in,2) should be equal to size(Acont1,2).');
    end
end
% % % % 

% % % % TODO (start) % % % %

% complete correlators
G1 = (-pi)*(KKi2r(ocont,Acont1)+1i*Acont1);
G2 = (-pi)*(KKi2r(ocont,Acont2)+1i*Acont2);
SE = G2./G1; % self-energy

% The imaginary part of the self-energy should be negative, to ensure the
% causality. However, in this numerical calculation, the imaginary part can
% be slightly positive over a narrow interval. To make the numerical result
% of SE physically sound, we cut off the imaginary part larger than
% "iSEmax" (defined above) to be bounded by "iSEmax".
oks = (imag(SE) > iSEmax);
SE(oks) = real(SE(oks)) + 1i*iSEmax;

if ~isempty(RhoV2in) % compute the improved estimate of the impurity spectral function
    % complete hybridization function, having both real and imaginary parts
    Deltaw = (-pi)*(KKi2r(ocont,RhoV2in)+1i*RhoV2in);
    Gimp = 1./(ocont - epsd - Deltaw - SE);
    varargout = {(-1/pi)*imag(Gimp)};
end

% % % % TODO (end) % % % %

end