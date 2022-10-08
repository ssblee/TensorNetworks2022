function [ts,M,Ovals,EE,dw] = tDMRG_Ex (M,Hs,O,Nkeep,dt,tmax)
% < Description >
%
% [ts,M,Ovals,EE,dw] = tDMRG_Ex (M,Hs,O,Nkeep,dt,tmax)
%
% The time-dependent DMRG (tDMRG) method that simulates real-time evolution
% of matrix product state (MPS), for an one-dimensional system whose
% Hamiltonian contains only nearest-neighbor terms, described by Hs{..}.
% The function uses the second-order Trotter decomposition: the time
% evolution operator for time step dt is decomposed into exp(-dt/2*Hodd) * 
% exp(-dt*Heven) * exp(-dt/2*Hodd). Those exponential terms in the
% decomposition are applied bond by bond; after acting an exponential term,
% the corresponding bond is truncated via SVD.
% This function also computes the expectation value of the local operator O
% is evaluated for every site and every time instances. 
%
% < Input >
% M : [cell] The initial state as the MPS. The length of M, i.e., numel(M),
%       defines the chain length. The leg convention of M{n} is as follows:
%
%    1      2   1      2         1        2
%   ---M{1}---*---M{2}---* ... *---M{end}---
%       |          |                 |
%       ^3         ^3                ^3
%
% Hs : [cell] Hamiltonian. Each cell element Hs{n} describes the two-site
%       interaction between site n and n+1. Thus, Hs(1:2:end) acts on odd
%       bonds; Hs(2:2:end) on even bonds. It should satisfy numel(M) ==
%       numel(Hs) + 1.
%       The leg convention of Hs{n} are as follows:
%
%       2      4      [legs 1 and 2 are for site n;
%       |      |       legs 3 and 4 are for site n+1]
%      [ Hs{n}  ]
%       |      |
%       1      3
%
% O : [matrix] Rank-2 tensor as a local operator acting on a site.
% Nkeep : [integer] Maximum bond dimension.
% dt : [numeric] Real time step size. Each real-time evolution by step dt
%       consists of three Trotter steps, exp(-dt/2*Hodd) * exp(-dt*Heven) *
%       exp(-dt/2*Hodd).
% tmax : [numeric] Maximum time range.
%
% < Output >
% ts : [numeric] Row vector of discrete time values.
% M : [cell] The final MPS after real-time evolution.
% Ovals : [matrix] Ovals(m,n) indicates the expectation value of local
%       operator O (input) at the site n and time ts(m).
% EE : [matrix] EE(m,n) indicates the entanglement entropy (with base 2) of
%       the MPS with respect to the bipartition cutting the bond between
%       the sites n and n+1, after applying the m-th row of time evolution
%       gates. For m's that are multiples of 3, EE(m/3,:) indicates the
%       entanglement entropy at time ts(m/3). Since the base 2 is chosen,
%       the value 1 of the entanglement entropy means one "ebit".
% dw : [matrix] Discarded weights (i.e., the sum of the squares of the
%       discarded singular values) after appying a row of time evolution
%       gates. dw(m,n) corresponds to the same bond and Trotter step
%       associated with EE(m,n).
%
% Written by S.Lee (Jun.19,2017); updated by S.Lee (Jun.22,2017)
% Updated by S.Lee (Jun.07,2019): Revised for SoSe 2019.
% Updated by S.Lee (Oct.08,2022): Revised for the course at SNU.


tobj = tic2;

% % % check the integrity of input
if numel(M) ~= (numel(Hs)+1)
    error('ERR: it should be: numel(M) == (numel(H)+1)');
elseif ~ismatrix(O)
    error('ERR: local operator O should be rank 2.');
end

for itN = (1:numel(Hs))
    if ~all(size(Hs{itN},1) == [size(Hs{itN},2) size(M{itN},3)])
        error(['ERR: The first and second legs of Hs{', ...
            sprintf('%i',itN),'} and the third leg of M{',sprintf('%i',itN), ...
            ' should have the same dimensions.']);
    elseif ~all(size(Hs{itN},3) == [size(Hs{itN},4) size(M{itN+1},3)])
        error(['ERR: The third and fourth legs of Hs{', ...
            sprintf('%i',itN),'} and the third leg of M{',sprintf('%i',itN+1), ...
            ' should have the same dimensions.']);
    end
end
% % % 

Nstep = ceil(tmax/dt);

% results
ts = dt*(1:Nstep);
Ovals = zeros(Nstep,numel(M));
EE = zeros(3*Nstep,numel(M)-1);
dw = zeros(size(EE));

% show information
fprintf('tDMRG : Real-time evolution with local measurements\n');
fprintf(['N = ',sprintf('%i',numel(M)),', Nkeep = ',sprintf('%i',Nkeep), ...
    ', dt = ',sprintf('%.4g',dt),', tmax = ',sprintf('%g',ts(end)), ...
    ' (',sprintf('%.4g',Nstep),' steps)\n']);

% generate the unitray operator exp(-it*H) for each two-site pairs
expH = cell(1,numel(Hs));
for it1 = (1:numel(Hs))
    if ~isempty(Hs{it1})
        sdim = [size(Hs{it1},1) size(Hs{it1},3)];
        Htmp = permute(Hs{it1},[1 3 2 4]);
        Htmp = reshape(Htmp,(sdim(1)*sdim(2))*[1 1]);
        if mod(it1,2) == 1
            ttmp = dt/2; % half time step for odd bonds, as the time evolution steps for odd bonds will happen twice
        else
            ttmp = dt;
        end
        [VH,DH] = eig(Htmp);
        eH = VH*diag(exp((-1i*ttmp)*diag(DH)))*VH';
        expH{it1} = reshape(eH,[sdim sdim]);
    end
end

disptime('Transform the MPS into right-canonical form.');
% since the first sweep is left-to-right, bring the input into
% right-canonical form, *without* truncation.
M = canonForm(M,0,[],[]);

disptime('Trotter steps: start');

for it1 = (1:3*Nstep)
% Here we use the 2nd order Trotter step exp(-dt/2*Hodd) * exp(-dt*Heven) *
% exp(-dt/2*Hodd). That is, for the case mod(it1,3) == 2, we act the
% unitary on even bonds. Otherwise, on odd bonds.
    expHtmp = cell(1,numel(Hs));
    if mod(it1,3) == 2 % even bonds
        expHtmp(2:2:end) = expH(2:2:end);
    else % odd bonds
        expHtmp(1:2:end) = expH(1:2:end);
    end
    
    % % % % TODO (start) % % % %
    % call local function tDMRG_1sweep which is written below in this file.
    
    
    if mod(it1,3) == 0
        % evaluate local expectation values
        
    end
    % % % % TODO (end) % % % %
        
    if (mod(it1,round(3*Nstep/10)) == 0) || (it1 == (3*Nstep))
        disptime(['#',sprintf('%i/%i',[it1/3,Nstep]), ...
            ' : t = ',sprintf('%g/%g',[ts(it1/3),ts(end)])]);
    end
end

toc2(tobj,'-v');
chkmem;

end

function [M,EE,dw] = tDMRG_1sweep (M,expH,Nkeep,isright)
% Apply exp(-it*Heven/odd), which is a row of two-site gates acting on
% either even or odd bonds, and then truncate bonds by using SVD. After
% applying this function, left-canonical state becomes right-canonical, and
% vice versa.
%
% < Input >
% M : [cell] Input MPS.
% expH : [cell] exp(-i*H*T) unitary operators for each bond. The length
%       should satisfy numel(expH) == numel(M)-1. And the every first (or
%       second) elements should be empty, since we act either even or odd
%       bonds at once.
% Nkeep : [numeric] Maximum bond dimension.
% isright : [logical] If true, we take left-to-right sweep. Otherwise, take
%       right-to-left sweep.
% 
% < Output >
% M : [cell] MPS after applying exp(-it*H) and truncating bonds.
% EE : [numeric vector] Entanglement entropy at each bond.
% dw : [numeric vector] Discarded weights when truncating the bond
%       dimensions.
    
N = numel(M);
EE = zeros(1,N-1);
dw = zeros(1,N-1);
Skeep = 1e-8;

% % % % TODO (start) % % % %
if isright % left -> right
    for it = (1:N-1)
        
    end
    
else % right -> left
    for it = (N-1:-1:1)
        
    end
    
end
% % % % TODO (end) % % % %

end

function Ovals = tDMRG_expVal (M,O,isleft)
% Expectation values of local operator O (acting on only one site) for
% given MPS.
%
% < Input >
% M : [cell] Input MPS.
% O : [matrix] Rank-2 operator acting on one site.
% isleft : [logical] If true, it means that the MPS M is in left-canonical
%       form. Otherwise, right-canonical form.
%
% < Output >
% Ovals : [vector] Expectation value of operator O. It will substitute the
%       rows of Ovals in the main function tDMRG.

N = numel(M);
Ovals = zeros(1,N);

MM = 1; % contraction of bra/ket tensors
if isleft % left-normalized
    for itN = (N:-1:1)
        T = permute(M{itN},[2 1 3]); % permute left<->right to make use of updateLeft
        T2 = contract(MM,2,2,T,3,1);
        T2 = contract(T2,3,3,O,2,2);
        Ovals(itN) = contract(conj(T),3,(1:3),T2,3,(1:3));
        MM = updateLeft(MM,2,T,[],[],T);
    end
else % right-normalized
    for itN = (1:N)
        T2 = contract(MM,2,2,M{itN},3,1);
        T2 = contract(T2,3,3,O,2,2);
        Ovals(itN) = contract(conj(M{itN}),3,(1:3),T2,3,(1:3));
        MM = updateLeft(MM,2,M{itN},[],[],M{itN});
    end 
end

end

