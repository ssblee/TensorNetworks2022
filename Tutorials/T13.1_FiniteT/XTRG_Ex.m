function [taus,lnZs,rho] = XTRG_Ex (Hs,dt,tmax,Nkeep,Nsweep)
% < Description >
%
% [taus,lnZs,rho] = XTRG_Ex (Hs,dt,tmax,Nkeep,Nsweep)
%
% Exponential tensor renormalization group (XTRG) method for simulating
% thermal density matrix as a matrix product operator (MPO). This function
% takes logarithmic imaginary-time steps and evaluates the logarithm of the
% partition function at those time instances.
% 
% < Input >
% Hs : [1 x N cell array] MPO representation of the Hamiltonian. Each Hs{n}
%       is a rank-4 tensor acting on site n. The order of legs of Hs{n} is
%       bottom-top-left-right, where the bottom (top) leg contracts to the
%       physical leg of bra (ket) tensor.
% dt : [numeric] Initial value of the imaginary time. In this function, the
%       first thermal density matrix in its MPO form is constructied via
%       linearization, i.e., I - dt*H, where I is the identity for the
%       whole Hilbert space and H is the Hamiltonian.
% tmax : [numeric] Maximum time to be reached at the end of iterations.
%       "tmax" equals to the target inverse temperature.
% Nkeep : [integer] Maximum bond dimension of the MPO form of the purified
%       thermal state.
% Nsweep : [integer] Number of round trips in the variational MPO
%       multiplication. "Nkeep" and "Nsweep" are directly forwarded to
%       "DMRG/mtimes_MPO.m"; see the description of those inputs therein.
%
% < Output >
% taus : [numeric] Time instances taken by logarithmically within the XTRG.
%       The actual imaginary time instances are -1i*taus. The elements of
%       "taus" equal to the inverse temperature values.
% lnZs : [numeric] Logarithms of the partition function measured at time
%       instances "taus".
% rho : [cell] The MPO representation of the thermal density matrix at the
%       last time instance.
% 
% Written by S.Lee (Oct.10,2022)

tobj = tic2;

N = numel(Hs);

% % % % TODO (start) % % % %
% initialize the first thermal density matrix, given by linearizing the
% exponential, I - dt*H, where H is the MPO Hamiltonian
rho = cell(1,N);
rho{1} = cat(4,getIdentity(Hs{1},2),-dt*Hs{1});
for itN = (2:N-1)
    
end
% % % % TODO (end) % % % %

Nstep = round(log2(tmax)-log2(dt));
taus = dt*(2.^(1:Nstep));
lnZs = zeros(size(taus)); % result

% show information
fprintf('Finite T: XTRG\n');
fprintf(['N = ',sprintf('%i',N),', Nkeep = ',sprintf('%i',Nkeep), ...
    ', Nsweep = ',sprintf('%i',Nsweep), ...
    ', dt = ',sprintf('%.4g',dt),', tmax = ',sprintf('%.4g',taus(end)), ...
    ' (',sprintf('%.4g',Nstep),' steps)\n']);

% main part of XTRG; iterative "squaring" of density matrix
for it1 = (1:Nstep)
    % Hermitian conjugate
    rho2 = rho;
    for itN = (1:N)
        rho2{itN} = permute(conj(rho2{itN}),[2 1 3 4]);
    end

    % MPO multiplication
    rho = mtimes_MPO_Ex(rho2,rho,Nkeep,Nsweep);

    % % % % TODO (start) % % % %
    % compute the trace, which leads to the partition function
    

    % add the log of the trace to the result array lnZs


    % add the result from the last iteration, since the MPO gets normalized
    % at every iteration (see below) to avoid divergence    
    

    % normalize the MPO


    % % % % TODO (end) % % % %

    if (mod(it1,round(Nstep/10)) == 0) || (it1 == Nstep)
        disptime(['#',sprintf('%i/%i',[it1,Nstep]), ...
            ' : t = ',sprintf('%.4g/%.4g',[taus(it1),taus(end)])]);
    end
end

toc2(tobj,'-v');
chkmem;

end
