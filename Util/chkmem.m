function varargout = chkmem (varargin)
% < Description >
%
% chkmem (['thval', ..] [, 'fmt', ..]) % display message on screen
% membyte = chkmem (['thval', ..] [, 'fmt', ..]) % suppress message
%
% Check the memory usage of a current MATLAB program (via system command
% 'ps' for Mac and Linux, or matlab command 'memory' in Windows). If the
% usage is larger than a threshold, show the message.
%
% < Option >
% 'thval', .. : [numeric] Threshold of memory usage (in byte) to show the
%           memory usage. (Default: 0, i.e., show always.)
% 'fmt', .. : [string] Format spec to convert number to string. (Default : '%.2f')
%
% < Output >
% membyte : [numeric] Memory usage in bytes. If output is specified,
%       the function does not display the message on screen.
%
% Written by S.Lee (2015)
% Updated by S.Lee (May 05,2017)
% Updated by S.Lee (Apr.17,2019): Cleaned code up.

% default option
thval = 0;
kunit = 1024; % currently, all the plaforms are using 1024 for RAM, to my knowledge
fmt = '%.2f';

% parse option
while ~isempty(varargin)
    switch varargin{1} 
        case 'fmt'
            fmt = varargin{2};
            varargin(1:2) = [];
        case 'thval'
            thval = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: cannot interpret option.');
    end
end

if strcmp(computer,'MACI64') || strcmp(computer,'GLNXA64') % for Mac or Linux
    mypid = feature('getpid');
    [isok,val] = system(['ps -p ',sprintf('%i',mypid),' -o rss | tail -n 1']);

    if isok == 0
        membyte = str2double(val)*kunit;
    else
        warning('WRN: failed to retreive memory usage from ''ps''.');
        membyte = 0;
    end
else % for Windows
    meminfo = memory;
    membyte = meminfo.MemUsedMATLAB;
end

if nargout
    varargout = {membyte};
else
    str = disptime(['Memory usage : ',byte2read(membyte,'k',kunit,'fmt',fmt),'\n']);
    if membyte > thval
        fprintf(str);
    end
end

end