%% MPO representation of time evolution operators
% Author: <https://cqm.snu.ac.kr Seung-Sup Lee>
%% 
% In the tutorial for tDMRG, we have time-evolved MPSs by applying rows of two-site 
% gates. Often, it is useful to represent a time evolution operator for a time 
% step as an MPO.
% 
% In this tutorial, we use two different approaches for constructing MPOs that 
% represent time evolution operators. As a concrete example, we consider the XY 
% spin-1/2 chain of even length $L$,
% 
% $$\hat{H}_{\mathrm{XY}}= -\sum_{\ell=1}^{L-1} (\hat{S}_{\ell,x} \hat{S}_{\ell+1,x} 
% + \hat{S}_{\ell,y} \hat{S}_{\ell+1,y})= -\frac{1}{2} \sum_{\ell=1}^{L-1} (\hat{S}_{\ell,+} 
% \hat{S}_{\ell+1,-} + \hat{S}_{\ell,-} \hat{S}_{\ell+1,+}) ,$$
% 
% and a time step $\Delta t = 0.01$.
%% Exercise (a): MPO for the first-order Trotterization
% In the first-order Trotter decomposition, the time evolution operator for 
% time step $\Delta t$ is split into $\exp (-\mathrm{i} \hat{H}_\mathrm{odd} \Delta 
% t) \exp (-\mathrm{i} \hat{H}_\mathrm{even} \Delta t)$, which has an error of 
% the order of $O(\Delta t^2)$. The decomposition is represented by two rows of 
% time evolution gates, as depicted in the upper part of the figure below:
% 
% 
% 
% We can decompose each two-site gate into two single-site tensors at different 
% sites, and contract the single-site tensors from different rows, to obtain an 
% MPO shown in the lower part of the figure above. Write a script that constructs 
% such an MPO for general even $L$. Verify your result by explicitly computing 
% $\exp (-\mathrm{i} \hat{H}_\mathrm{odd} \Delta t) \exp (-\mathrm{i} \hat{H}_\mathrm{even} 
% \Delta t)$ for a small system, say $L = 6$.
%% Exercise (b): MPO for the first-order Taylor expansion
% On the other hand, we can make another first-order approximation of $\exp(-\mathrm{i} 
% \hat{H}_\mathrm{XY} \Delta t)$, which as an error of the order of $O(\Delta 
% t^2)$, by using the Taylor expansion: namely, $\exp(-\mathrm{i} \hat{H}_\mathrm{XY} 
% \Delta t) \approx \hat{I} - \mathrm{i} \hat{H}_\mathrm{XY} \Delta t$. Write 
% a script that constructs such an MPO for general even $L$. Verify your result 
% by explicitly computing $\hat{I} - \mathrm{i} \hat{H}_\mathrm{XY} \Delta t$ 
% for a small system, say $L = 6$.
% 
% (_Hint_: Refer to Sec. 5.2 of Schollwoeck2011 [<https://www.sciencedirect.com/science/article/abs/pii/S0003491610001752?via%3Dihub 
% U. Schollwöck, Ann. Phys. *326*, 96 (2011)>].)