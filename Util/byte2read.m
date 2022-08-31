function str = byte2read (inbyte, varargin)
% < Description >
%
% str = byte2read (inbyte [, 'fmt', ..] [, 'k', ..])
%
% Convert input data size in bytes to human readable size using TB, GB, MB,
% KB, B (decimal, 1000) or 'TiB','GiB','MiB','KiB','B' (binary, 1024).
%
% < Input >
% inbyte : [numeric] Byte.
%
% < Option >
% 'fmt', .. : [string] Format spec to convert number to string. (Default : '%.2f')
% 'k', .. : [integer] Set 1000 to use kB, MB, GB, ...; or set 1024 to use
%       KiB, MiB, GiB, ...) 
%       (Default: 1000 for linux & mac, 1024 for Windows)
%
% < Output >
% str : [string] String describing data size.
%
% Written by S.Lee (2015)
% Updated by S.Lee (May 05,2017)
% Updated by S.Lee (Apr.17,2019): Cleaned code up.

% default value of options
fmt = '%.2f';
if ismember(computer,{'GLNXA64','MACI64'})
    kunit = 1000;
else % for windows
    kunit = 1024;
end

% parse options
while ~isempty(varargin)
    switch varargin{1}
        case 'fmt'
            fmt = varargin{2};
            varargin(1:2) = [];
        case 'k'
            kunit = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: cannot interpret option.');
    end
end

% sanity check for kunit
if ~any(kunit == [1000 1024])
    error('ERR: Value for ''k'' should be either 1000 or 1024.');
end

if kunit == 1000 % Linux and Mac
    strunits = {'B','kB','MB','GB','TB'};
else %if kunit == 1024 % Windows
    strunits = {'B','KiB','MiB','GiB','TiB'};
end

kpow = floor(log(inbyte)/log(kunit));
if kpow > (numel(strunits)-1)
    kpow = numel(strunits)-1;
elseif kpow < 0
    kpow = 0;
end

str = [sprintf(fmt,inbyte/(kunit^kpow)),strunits{kpow+1}];

end