%% [Solutions] Controlled bond expansion (CBE) DMRG for ground state search
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% Solution to Exercise (a): Complete the function for CBE-DMRG
% The complete version is uploaded as |DMRG/DMRG_GS_CBE.m|. Compare with your 
% implementation!
%% Solution to Exercise (b): Error analysis
% We repeat the CBE-DMRG calculations for the free spinful fermions on a tight-binding 
% chain, for different values of |Nkeep|. We use the same initial MPS obtained 
% from the iterative diagonalization, with the lowest value of |Nkeep|.

clear

% system parameters
t = 1; % hopping amplitude
L = 30; % number of sites in a chain

% CBE-DMRG parameters
Nkeep = (50:50:300);
Nsweep = 8;
delta = 0.1;

E0_exact = nonIntTB(ones(L-1,1)*t)*2; % *2 due to two spins

% Local operators
[F,Z,S,I] = getLocalSpace('FermionS');

% % MPO formulation of Hamiltonian
% Hamiltonian tensor for each chain site
Hloc = cell(6,6);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = Z*F(:,:,1);
Hloc{3,1} = Z*F(:,:,2);
Hloc{4,1} = Hloc{2,1}';
Hloc{5,1} = Hloc{3,1}';
Hloc{6,2} = -t*F(:,:,1)';
Hloc{6,3} = -t*F(:,:,2)';
Hloc{6,4} = -t*F(:,:,1);
Hloc{6,5} = -t*F(:,:,2);
Hloc{6,6} = I;

Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)]));

% full chain
Hs = cell(1,L);
Hs(:) = {Hloc};
Hs{1} = Hs{1}(:,:,end,:); % choose the last components of the left leg
Hs{end} = Hs{end}(:,:,:,1); % choose the first components of the right leg

% initial MPS
Minit = cell(1,L);

Hprev = 1; % initialize Hamiltonian with 1, as we will use MPO
Aprev = 1; % identity tensor for the dummy leg

for itN = (1:L)
    % add new site
    Anow = getIdentity(Aprev,2,I,2,[1 3 2]);
    Hnow = updateLeft(Hprev,3,Anow,Hs{itN},4,Anow);

    Hmat = Hnow(:,:,1);
    [V,D] = eig((Hmat+Hmat')/2);
    [D,ids] = sort(diag(D),'ascend');
    if itN < L
        Ntr = min([numel(D);Nkeep(1)]);
    else
        Ntr = 1;
    end
    V = V(:,ids(1:Ntr));
    
    Anow = contract(Anow,3,2,V,2,1,[1 3 2]);

    Minit{itN} = Anow;
    
    Hprev = contract(Hnow,3,2,V,2,1);
    Hprev = contract(V',2,2,Hprev,3,1,[1 3 2]);
    Aprev = Anow;
end
%% 
% Run calculations.

% result arrays
E0s = zeros(1,numel(Nkeep));
Eiters = cell(1,numel(Nkeep));
dws = zeros(1,numel(Nkeep));
varEs = zeros(1,numel(Nkeep));

for itk = (1:numel(Nkeep))
    [M,E0s(itk),Eiters{itk},~,dw] = DMRG_GS_CBE (Minit,Hs,Nkeep(itk),Nsweep,delta);
    dws(itk) = sum(dw(:,end));
    varEs(itk) = varE_2site(M,Hs);
end
%% 
% Let's see how the variational energy decreases with iterations, for different 
% values of |Nkeep|.

legs = cell(1,numel(Nkeep));
clrs = parula(round(numel(Nkeep)*1.2));

figure;
hold on;
for itk = (1:numel(Nkeep))
    plot((1:numel(Eiters{itk}))/(L-1),Eiters{itk}(:)-E0_exact, ...
        'Color',clrs(itk,:),'LineWidth',1);
    legs{itk} = sprintf('Nkeep = %i',Nkeep(itk));
end
hold off;
set(gca,'XScale','Linear','YScale','log','FontSize',13,'LineWidth',1);
xlim([0 2*Nsweep]);
grid on;
xlabel('# of sweeps');
ylabel('Ground-state energy error');
legend(legs);
%%
figure;
hold on;
plot(dws,E0s-E0_exact,'LineWidth',1,'Marker','x','MarkerSize',12);
plot(varEs,E0s-E0_exact,'LineWidth',1,'Marker','+','MarkerSize',12);
set(gca,'XScale','log','YScale','log','FontSize',13,'LineWidth',1);
grid on;
ylabel('Ground-state energy error');
legend({'vs. discarded weight','vs. 2-site variance'},'Location','northwest');
%% 
% We find that the ground-state energy error scales linearly against both the 
% discarded weight and the two-site variance.