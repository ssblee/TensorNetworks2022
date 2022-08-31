function tobj = tic2
% < Description >
%
% tobj = tic2
%
% Create a timer object (struct) to measure both real and CPU times.
%
% < Output >
% tobj : [struct] Timer object. This will be the input for 'toc2' function.
%   .real : [uint8] Internal timer for built-in functions 'tic' and 'toc'.
%   .cpu : [double] Time in seconds after the current session of MATLAB
%           started.
%
% Written by S.Lee (May 05,2017)

tobj = struct;
tobj.real = tic;
tobj.cpu = cputime;

end