function [Nsq,Sv] = TRG_GILT_Square_Ex (T,rgstep,Nkeep,varargin)
% < Description >
%
% [Nsq,Sv] = TRG_GILT_Square_Ex (T,rgstep,Nkeep [, options])
%
% Perform the tensor renormalization group (TRG) calculation of contracting
% a tensor network on an exponentially large square lattice. If the option
% 'GILT' is set 'true', the graph-independent local truncations (GILT) are
% made to effectively remove local correlations, largely enhancing the
% accuracy of the TRG-type coarse graining. When GILT is used, the
% contraction scheme is called GILT-TNR. Refer to Hauru2018 [M. Hauru, C.
% Delcamp, and S. Mizera, Phys. Rev. B 97, 045111 (2018)] for the details
% of the GILT algorithm.
%
% < Input >
% T : [rank-4 tensor] As a tensor network to be contracted, this function
%       considers only the case that it has a one-site unit cell. The input
%       "T" is a rank-4 tensor that represents the unit cell. If one
%       computes the norm of a quantum state in its PEPS representation,
%       then "T" is the reduced tensor made of a rank-5 ket tensor and its
%       conjugate (i.e., rank-5 bra tensor). The legs of T are ordered as
%       left-up-down-right.
% rgstep : [numeric] Number of TRG coarse graining steps. This function
%       considers a lattice of (4*2^rgstep) sites with periodic boundary
%       conditions.
% Nkeep : [numeric] Maximum bond dimension (\chi) for tensors.
%
% < Output >
% 'GILT', .. : [numeric] Maximum number of GILT iterations per bond. If set
%       zero, then only pure TRG calculations will be done.
%       (Default: 0, no GILT)
% 'epsilon', .. : [numeric] A parameter for determining the elements
%       t^\prime in Eq. (31) of Hauru2018.
%       (Default: 8e-7, as noted in the caption of Fig. 5 of Hauru2018)
% 'Skeep_GILT', .. : [numeric] The lower bound to the singular values to be
%       kept, in the SVD of the R^\prime tensor (sitting on the bond under
%       a GILT) into two tensor that are to be absorbed into the
%       neighboring tensors.
%       (Default: 8e-10, based on 'epsilon' as noted in the caption of Fig.
%       5 of Hauru2018 and the authors' code implementation given in
%       https://github.com/Gilt-TNR/Gilt-TNR/blob/master/GiltTNR2D.py )
% 'tol', .. : [numeric] At each iteration within GILT, a new part of the
%       R^\prime tensor is split into the left and right parts [given by
%       U*sqrt(S) and sqrt(S)*Vd, respectively, where [U,S,Vd] =
%       svdTr(...)]. Then those left and right parts are absorbed into the
%       existing left and right parts. If the singular values S(n) are
%       close to 1, then it means that the GILT iteration converges, since
%       the R^\prime tensor does not change much. The optional parameter
%       "tol" is used to check whether all the singular values satisfy
%       abs(S - 1) < tol. Once converged, the left and right parts are
%       absorbed into the neighboring tensors.
%       (Default: 1e-2)
%
% < Output >
% Nsq : [numeric] (contraction of the tensor network)^(1/(number of lattice
%       sites)). If the input "T" is a reduced tensor for a PEPS, then
%       "Nsq" indicates the squared norm of the PEPS per lattice site.
% Sv : [numeric] Sv(:,n) indicates the singular values in the SVD of a
%       tensor on the upper-left corner of a plaquette, at the n-th TRG
%       coarse-graining step.
%
% Written by S.Lee (Nov.21,2022): for the lecture course at SNU.


tobj = tic2;

% lattice parameters
z = 4; % coordination number = # of nearest neighbors
nuc = 4; % # of lattice sites that define a unit cell; here 2 x 2 = 4

% default values of optional parameters
N_GILT = 0;
epsilon = 8e-7;
Skeep_GILT = 8e-10;
tol = 1e-2;


% parse options
while ~isempty(varargin)
    if numel(varargin) < 2
        error('ERR: Option should be set by a pair of option name and value.');
    end
    switch varargin{1}
        case 'GILT'
            N_GILT = varargin{2};
            varargin(1:2) = [];
        case 'epsilon'
            epsilon = varargin{2};
            varargin(1:2) = [];
        case 'Skeep_GILT'
            Skeep_GILT = varargin{2};
            varargin(1:2) = [];
        case 'tol'
            tol = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: Unknown input.');
    end
end


% use a 2 x 2 unit cell to use GILT
Ts = cell(nuc,1);
Ts(:) = {permute(T,[4 2 1 3])}; % re-order the legs to a counter-clockwise order, i.e., right-up-left-down

% We use the following convention of the 2 x 2 unit cell and the leg order:
%
%      2 |            | 4
%  3     |     1      |     3
% ---- Ts{1} ------ Ts{2} ----
%        |            |
%      4 |            | 2
%  1     |     3      |     1
% ---- Ts{4} ------ Ts{3} ----
%        |            |
%      2 |          4 |

% result data
Nsq = 1; % initialize with 1
Sv = cell(1,rgstep);

if N_GILT > 0
    disp(['GILT-TNR on a square lattice: rgstep = ',sprintf('%g',rgstep), ...
        ', Nkeep = ',sprintf('%g',Nkeep),', epsilon = ',sprintf('%.4g',epsilon)]);
else
    disp(['TRG on a square lattice: rgstep = ',sprintf('%g',rgstep), ...
        ', Nkeep = ',sprintf('%g',Nkeep)]);
end

for itr = (1:rgstep)
    % pull out normalization prefactors so that the contraction result does
    % not diverge
    normfac = max(abs(Ts{1}(:)));
    Nsq = Nsq*exp(log(normfac)/(2^(itr-1))); % contribution per site; currently each tensor represents 2^(itr-1) sites

    % normalize with the same factor
    for it2 = (1:numel(Ts))
        Ts{it2} = Ts{it2}/normfac;
    end

    if N_GILT > 0 % GILT
        itgs = zeros(1,numel(Ts)); % # of GILT iterations
        for it1 = (1:numel(Ts)) % update each of four independent bonds along the plaquette
            [Ts,itgs(it1)] = GILT_Square (Ts,epsilon,Skeep_GILT,tol,N_GILT);

            % rotate the plaquette counter-clockwise by 90 degree
            Ts = Ts([(2:end) 1]);
            for it2 = (1:numel(Ts))
                Ts{it2} = permute(Ts{it2},[(2:z) 1]);
            end
        end
    end

% "bipartitions" of Ts; we consider the following bipartitions:
%
%                1 |               | 3
%                  |               | 
%                Vs{1} ---   --- Vs{2}
%               2 /     3    1      \ 2
%                /                   \
%        3    2 /                     \ 2    1
%        --- Us{1}                   Us{2} ---
%              |                       |
%              | 1                   3 |
%
%              | 3                   1 |
%        1     |                       |    3
%        --- Us{4}                   Us{3} ---
%             2 \                   2 /
%                \                   /
%               2 \     1   3     2 /
%                Vs{4} ---   --- Vs{3}
%                  |               | 
%                3 |               | 1 

    Us = cell(numel(Ts),1);
    Vs = cell(numel(Ts),1);

    % % % % TODO (start) % % % %
    

    % % % % TODO (end) % % % %

    strtmp = ['#',sprintf('%02i/%02i',[itr,rgstep])];
    
    if N_GILT > 0
        strtmp = [strtmp,', GILT iterations = [',sprintf('%i, ',itgs)];
        strtmp = [strtmp(1:end-2),']'];
    end

    disptime(strtmp);
end

% % % % TODO (start) % % % %
% last step: contract four tensors exactly

% % % % TODO (end) % % % %

Sv = cell2mat(cellfun(@(x) [x;nan(Nkeep-numel(x),1)], Sv, 'UniformOutput', false));

strs = ['Contraction result per site = ',sprintf('%.6g',Nsq)];
disptime(strs);

toc2(tobj,'-v');

end

function [U,V,S] = TRG_Square_split (T,Nkeep,idU,pU,pV)
% < Description >
%
% [U,V,S] = TRG_Square_split (T,Nkeep,idU,pU,pV)
%
% Wrapper for svdTr. After the SVD, the singular value tensor S is split
% into sqrt(S)*sqrt(S), and each par is absorbed into U and V. And then U
% and V are permuted according to pU and pV, respectively.

% % % % TODO (start) % % % %

% % % % TODO (end) % % % %

end

function [Ts,itg] = GILT_Square (Ts,epsilon,Skeep_GILT,tol,N_GILT)
% < Decription >
%
% [Ts,itg] = GILT_Square (Ts,epsilon,Skeep_GILT,tol,N_GILT)
%
% Optimize a bond between Ts{1} and Ts{2} by GILT. (Ts{n}'s are placed on a
% square plaquette in the clockwise order; see the diagram inside the main
% function above.)
%
% < Input >
% Ts : [cell array] The tensors sitting along the plaquette.
% epsilon, Skeep_GILT, tol, N_GILT : GILT parameters. Refer to the
%       documentation of the main function for details.
%
% < Output >
% Ts : [cell array] Ts(1:2) are updated, while the rest remains the same.
% itg : [numeric] Number of GILT iterations.

% Partial results of the environment tensor, made of all of Ts(:)
[Uenv,Senv] = GILT_Square_env (Ts);

% Cumulative products of left- and right-part tensors. They will be
% contracted to the neighboring tensors attach to the bond that is
% GILT-updated.
Ucum = eye(size(Uenv,1));
Vcum = eye(size(Uenv,1));

iscvg = false; % true if the iterative GILT converges

for itg = (1:N_GILT)
    Senv = Senv/sum(Senv); % normalize

    % % % % TODO (start) % % % %

    % Trace of each slice of "Uenv"
    

    % Determine the elements of t^\prime_i, according to Eq. (31) of
    % Hauru2018
    

    % A new part of R^\prime, to be split and absorbed

    
    % Decompose "Rp" into the left and right parts via SVD


    % Absorb the left and right parts into the existing tensors


    % Check convergence
    if max(abs(S2-1)) < tol
        iscvg = true;
        break;
    end

    % Update the environment contributions (Uenv, Senv): here the full
    % computation of the environment tensor is not necessary!
    

    % % % % TODO (end) % % % %
end

if ~iscvg
    disp(['WRN: GILT iterations did not converge after ',sprintf('%i',N_GILT),' iteration(s).']);
end

% Absorb the left and right parts of R^\prime into Ts{1} and Ts{2}
Ts{1} = contract(Ts{1},4,1,Ucum,2,1,[4 1 2 3]);
Ts{2} = contract(Ts{2},4,1,Vcum,2,2,[4 1 2 3]);

end

function [Uenv,Senv] = GILT_Square_env (Ts)
% < Description >
%
% [Uenv,Senv] = GILT_Square_env (Ts)
%
% Compute the partial result of SVD of the environment tensor E, described
% in Eqs. (23) and (24) of Hauru2018. As specified in the footnote 4 in the
% paper, since we only need U and S, we can simplify the calculation by
% eigen-decomposing E*E', where in the contraction E*E' all the legs
% (except for the legs for the update to be updated) are contracted.

% % % % TODO (start) % % % %

% % % % TODO (end) % % % %

end
