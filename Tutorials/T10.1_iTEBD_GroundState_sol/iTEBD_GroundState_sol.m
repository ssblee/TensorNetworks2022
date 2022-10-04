%% [Solution] iTEBD: Ground state search
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% Solution to Exercise (a): Complete the function for Vidal's original iTEBD
% Check out the funciton |iTEBD_GS_Vidal.m| under the |DMRG| sub-directory. 
% Compare with your implementation of |iTEBD_GS_Vidal_Ex.m|!
%% Solution to Exercise (b): Check the energy convergence
% In the demonstration, we use the following choice of parameters:

clear

% iTEBD parameters
Nkeep = 30;
tau_ini = 1; % initial imaginary time step size
tau_fin = 1e-6; % final imaginary time step size
Nstep = 2e3;
%% 
% In this solution, I will see how iTEBD's performance changes with chaging 
% (i) |tau_fin|, (ii) |Nstep|, or (iii) |Nkeep|. We use the same random MPS as 
% the initial state, for fair comparison.

% Local operators
[S,I] = getLocalSpace('Spin',1);
% Heisenberg interaction as two-site gate S*S'
H = contract(S,3,3,permute(conj(S),[2 1 3]),3,3);

% Initialize with random Lambda and Gamma
Lambda_init = cell(1,2);
Gamma_init = cell(1,2);
for itn = (1:numel(Lambda_init))
    Lambda_init{itn} = rand(Nkeep,1);
    Gamma_init{itn} = rand(Nkeep,Nkeep,size(I,2));
end

% "numerically exact" value by White & Huse
Eexact = -1.401484039;
%% 
% First, let's see the dependence of |tau_fin|.

tau_fins = 10.^(-4:-1:-8);
Eiters = cell(size(tau_fins));

for itt = (1:numel(Eiters))
    taus = tau_ini*((tau_fins(itt)/tau_ini).^linspace(0,1,Nstep));
    % iTEBD ground state search
    [~,~,Eiters{itt}] = iTEBD_GS_Vidal(Lambda_init,Gamma_init,H,Nkeep,taus);
end
figure;
legs = cell(size(Eiters));
hold on;
for itt = (1:numel(Eiters))
    Eiter2 = reshape(permute(Eiters{itt},[2 1 3]), ...
        [size(Eiters{itt},2)*size(Eiters{itt},1) size(Eiters{itt},3)]);
    plot((1:size(Eiter2,1)).'/2,mean(Eiter2,2)-Eexact,'LineWidth',1);
    legs{itt} = ['$\tau_\mathrm{fin} = 10^{', ...
        sprintf('%i',log10(tau_fins(itt))),'}$'];
end
hold off;
set(gca,'YScale','log','LineWidth',1,'FontSize',13);
xlabel('Step');
ylabel('Energy per bond - exact energy');
legend(legs{:},'Interpreter','latex');
grid on;
%% 
% We find that for the largest |tau_fin|'s ($10^{-4}$ and $10^{-5}$) the energy 
% has not fully converged. It is because larger step sizes in the last iterations 
% lead to larger Trotterization error. For smaller |tau_fin|'s ($10^{-7}$ and 
% $10^{-8}$) , the energy approaches faster, but the final value is larger than 
% that for $10^{-6}$. It is because the total imaginary time interval is smaller 
% than that for $10^{-6}$:

disp(sum(tau_ini*((tau_fins(3)/tau_ini).^linspace(0,1,Nstep)))); % 1e-6
disp(sum(tau_ini*((tau_fins(5)/tau_ini).^linspace(0,1,Nstep)))); % 1e-8
%% 
% How about the dependence of |Nstep|?

Nsteps = (1e3:500:3e3);
Eiters = cell(size(Nsteps));

for itt = (1:numel(Eiters))
    taus = tau_ini*((tau_fin/tau_ini).^linspace(0,1,Nsteps(itt)));
    % iTEBD ground state search
    [~,~,Eiters{itt}] = iTEBD_GS_Vidal(Lambda_init,Gamma_init,H,Nkeep,taus);
end
figure;
legs = cell(size(Eiters));
hold on;
for itt = (1:numel(Eiters))
    Eiter2 = reshape(permute(Eiters{itt},[2 1 3]), ...
        [size(Eiters{itt},2)*size(Eiters{itt},1) size(Eiters{itt},3)]);
    plot((1:size(Eiter2,1)).'/2,mean(Eiter2,2)-Eexact,'LineWidth',1);
    legs{itt} = ['$N_\mathrm{step} = ',sprintf('%i',Nsteps(itt)),'$'];
end
hold off;
set(gca,'YScale','log','LineWidth',1,'FontSize',13);
xlabel('Step');
ylabel('Energy per bond - exact energy');
legend(legs{:},'Interpreter','latex');
grid on;
%% 
% For smaller |Nstep|'s, the energy approaches faster, but the final value is 
% larger than that for the largest |Nstep| (= 3000). It is because the total imaginary 
% time interval is smaller.
%% 
% As the last part of this analysis, we change |Nkeep|.

Nkeeps = (10:10:50);
Eiters = cell(size(Nkeeps));
taus = tau_ini*((tau_fin/tau_ini).^linspace(0,1,Nstep));

for itt = (1:numel(Eiters))
    % iTEBD ground state search
    [Lambda,Gamma,Eiters{itt}] = iTEBD_GS_Vidal(Lambda_init,Gamma_init,H,Nkeeps(itt),taus);
end
figure;
legs = cell(size(Eiters));
hold on;
for itt = (1:numel(Eiters))
    Eiter2 = reshape(permute(Eiters{itt},[2 1 3]), ...
        [size(Eiters{itt},2)*size(Eiters{itt},1) size(Eiters{itt},3)]);
    plot((1:size(Eiter2,1)).'/2,mean(Eiter2,2)-Eexact,'LineWidth',1);
    legs{itt} = ['$N_\mathrm{keep} = ',sprintf('%i',Nkeeps(itt)),'$'];
end
hold off;
set(gca,'YScale','log','LineWidth',1,'FontSize',13);
xlabel('Step');
ylabel('Energy per bond - exact energy');
legend(legs{:},'Interpreter','latex');
grid on;
figure;
hold on;
plot((1:numel(Lambda{1})).',Lambda{1},'x','LineWidth',1,'MarkerSize',10);
plot((1:numel(Lambda{2})).',Lambda{2},'+','LineWidth',1,'MarkerSize',10);
hold off;
set(gca,'YScale','log','LineWidth',1,'FontSize',13);
ylabel('Singluar values');
grid on;
%% 
% We find that the energy gets lower as |Nkeep| increases. The differences in 
% the final energy values are due to the contibution from small singular values 
% and the corresponding singular vectors (that would be discarded if |Nkeep| is 
% smaller).