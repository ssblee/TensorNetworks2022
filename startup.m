% < Description >
% startup.m for the lecture course "Tensor Networks". It
% (i) resets the path; 
% (ii) add the subdirectories of the directory (in which this 'startup.m'
%      file lies) to the path;
% (iii) clear variables in workspace;
% (iv) and show a startup message.
%
% Written by S.Lee (Apr.21,2017)
% Rewritten by S.Lee (Aug.31,2022)
% Updated by S.Lee (Nov.23,2022): Add one more sub-directory "ML".

% % Reset path
% temporarily turn off warning; if warning is not turned off above, MATLAB
% might give warning
warning('off','MATLAB:dispatcher:pathWarning');
restoredefaultpath; % reset path to default
warning('on','MATLAB:dispatcher:pathWarning'); % turn on warning back
% NOTE: if you want to keep custom path besides TN folder, comment out the
% above three commands.

% % Add to path
fpath = fileparts(mfilename('fullpath')); % the TN directory in which this "startup.m" lies
addpath(fpath);
dirnames = {'MyWork','Util','Tensor','NRG','DMRG','PEPS','ML','Tutorials'}; % sub-directories
for it1 = (1:numel(dirnames))
    fpath2 = [fpath,filesep,dirnames{it1}];
    if exist(fpath2,'dir')
        addpath(genpath(fpath2)); % add sub-directories to path
    end
    addpath(genpath(fpath2));
end

% % Clear variables in memory
clear

% % startup messeage
fprintf('startup.m | "Tensor Networks" tutorial material\n');
chkmem;
