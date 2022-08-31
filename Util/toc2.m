function [realt,cput,avgcore] = toc2 (tobj,varargin)
% < Description >
%
% [realt,cput,avgcore] = toc2 (tobj [,'-v'])
%
% Count real time (in sec), CPU time (in sec), and the average number of
% cores used after the timer object 'tobj' is created by function 'tic2'.
%
% < Input >
% tobj : [struct] Timer object. Refer to 'tic2' for detail.
%
% < Option >
% '-v' : Option to display results.
%
% < Output >
% realt : [double] Real time in seconds.
% cput : [double] CPU time in seconds.
% avgcore : [double] Average number of cores used. avgcore = cput/realt.
%
% Written by S.Lee (May 05,2017)
% Updated by S.Lee (Apr.17,2019): Cleaned code up.

oshow = false;

% % parsing option
while ~isempty(varargin)
    switch varargin{1}
        case '-v'
            oshow = true;
            varargin(1) = [];
        otherwise
            error('ERR: Unknown option.');
    end
end
% % % 

realt = toc(tobj.real);
cput = cputime - tobj.cpu;
avgcore = cput/realt;

if oshow
    fprintf(['Elapsed time: ',sprintf('%.4g',realt),'s, CPU time: ', ...
        sprintf('%.4g',cput),'s, Avg # of cores: ', ...
        sprintf('%.4g',avgcore),'\n']);
end

end