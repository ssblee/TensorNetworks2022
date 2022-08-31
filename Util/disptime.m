function varargout = disptime (varargin)
% < Description >
% 
% disptime                  % display timestamp only
% disptime (str)            % display message (timestamp + str) on screen
% timeval = disptime (str)  % suppress message
%
% Display "current time (i.e. timestamp) + input string". Write the string
% on the file if  file name is given. Timestamp format is 'yy-mm-dd
% HH:MM:SS'.
%
% < Input >
% str : [string] String to display together with timestamp.
%
% < Output >
% timeval : [string] Input string + current time. If output is specified
%       (e.g. command such as "res = disptime;" is used), the function does
%       not display the message on screen.
%
% Written by S.Lee (2015)
% Updated by S.Lee (May 05,2017)
% Updated by S.Lee (Apr.17,2019): Cleaned code up.

if nargin
    if ischar(varargin{1})
        str = varargin{1};
    elseif isnumeric(varargin{1})
        str = num2str(varargin{1});
    else
        error('ERR: Unknown input.');
    end
        
    res = [datestr(now,'yy-mm-dd HH:MM:SS'),' | ',str];
else
    res = datestr(now,'yy-mm-dd HH:MM:SS');
end

if nargout
    varargout = {res};
else
    disp(res);
end

end