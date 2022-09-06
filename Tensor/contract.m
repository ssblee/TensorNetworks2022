function C = contract(A,rankA,idA,B,rankB,idB,varargin)
% < Description >
%
% C = contract(A,rankA,idA,B,rankB,idB [,idC] [,'-v']) 
%
% Contract tensors A and B. The legs to be contracted are given by idA
% and idB.
%
% < Input >
% A, B : [numeric array] Tensors.
% rankA, rankB : [integer] Rank of tensors. Since MATLAB removes the last
%       trailing singleton dimensions, it is necessary to set rankA and
%       rankB not to miss the legs of size 1 (or bond dimension 1,
%       equivalently).
% idA, idB : [integer vector] Indices for the legs of A and B to be
%        contracted. The idA(n)-th leg of A and the idB(n)-th leg of B will
%        be contracted, for all 1 <= n <= numel(idA). idA and idB should
%        have the same number of elements. If they are both empty, C will
%        be given by the direct product of A and B.
% 
% < Option >
% idC : [integer array] To permute the resulting tensor after contraction,
%       assign the permutation indices as idC. If the dummy legs are
%       attached (see the description of C below), this permutation is
%       applied *after* the attachment.
%       (Default: no permutation)
% '-v' : Show details. Set this option to see how the legs of C are related
%       to the legs of A and B.
%
% < Output >
% C : [numeric array] Contraction of A and B. If idC is given, the
%       contracted tensor is permuted accordingly. If the number of open
%       legs are smaller than 2, the dummy legs are assigned to make the
%       result array C be two-dimensional.
%
% Written by S.Lee (Apr.30,2017)
% Updated by S.Lee (Apr.25,2019): Revised code.


% % default values of options
idC = []; % permutation of the contracted tensor (default: no permutation)
oshow = false; % option to show details

% % % parsing options
while ~isempty(varargin)
    if isnumeric(varargin{1})
        idC = varargin{1};
        varargin(1) = [];
    elseif isequal(varargin{1},'-v')
        oshow = true;
        varargin(1) = [];
    else
        disp(varargin{1});
        error('ERR: Unkown option.')
    end
end
% % %

% % % option -v
if oshow
    tobj = tic2;
end
% % % 

% % check the integrity of input and option
Asz = size(A); Bsz = size(B); % size of tensors
if rankA < numel(Asz)
    error('ERR: Input ''rankA'' is smaller than the rank of other input ''A''.');
elseif rankB < numel(Bsz)
    error('ERR: Input ''rankB'' is smaller than the rank of other input ''B''.');    
end
% append 1's, in case that the tensor is indeed high-rank but the trailing
% singleton dimensions (i.e. dimensions of size 1)
Asz = [Asz,ones(1,rankA-numel(Asz))];
Bsz = [Bsz,ones(1,rankB-numel(Bsz))];

if numel(idA) ~= numel(idB)
    error('ERR: Different # of leg indices to contract for tensors A and B.');
end
if ~isempty(idA) % sanity check of idA and idB
    % reshape idA and idB as row vectors
    idA = idA(:).'; 
    idB = idB(:).';
    
    % logical array to check that idA has unique elements
    oks = bsxfun(@eq, (1:rankA).', idA);
    if ~all(any(oks,1)) ... % if elements are not in the range [1, rankA]
            || any(sum(oks,2) > 1) % if elements are not unique
        error(['ERR: idA = [',sprintf('%.4g ',idA), ...
            '] should consist of unique integers in the range (1:rank(A)).']);
    end
    % likewise for idB
    oks = bsxfun(@eq, (1:rankB).', idB);
    if ~all(any(oks,1)) ... % if elements are not in the range [1, rankA]
            || any(sum(oks,2) > 1) % if elements are not unique
        error(['ERR: idB = [',sprintf('%.4g ',idB), ...
            '] should consist of unique integers in the range (1:rank(B)).']);
    end
    
    if ~all(Asz(idA) == Bsz(idB))
        error(['ERR: Dimensions of A to be contracted = [',sprintf('% i',Asz(idA)), ...
            '] do not match with those of B = [',sprintf('% i',Bsz(idB)),'].']);
    end 
end

% check whether permutation option is correct
if ~isempty(idC) && (numel(idC) < (rankA + rankB - 2*numel(idA)))
    error(['ERR: # of indices for the permutation after contraction (= ', ...
        sprintf('%i',numel(idC)), ...
        ') is different from # of legs after contraction (= ', ...
        sprintf('%i',(rankA + rankB - 2*numel(idA))),').']);
end
% % % % % 

% % % % Main computational part (start) % % % %
% indices of legs *not* to be contracted
idA2 = (1:rankA); idA2(idA) = []; 
idB2 = (1:rankB); idB2(idB) = [];

% reshape tensors into matrices with "thick" legs
A2 = reshape(permute(A,[idA2 idA]),[prod(Asz(idA2)) prod(Asz(idA))]); % note: prod([]) == 1
B2 = reshape(permute(B,[idB idB2]),[prod(Bsz(idB)) prod(Bsz(idB2))]);
C2 = A2*B2; % matrix multiplication

% size of C
if (numel(idA2) + numel(idB2)) > 1
    Cdim = [Asz(idA2),Bsz(idB2)];
else
    % place dummy legs x of singleton dimension when all the legs of A (or
    % B) are contracted with the legs of B (or A)
    Cdim = [1,1];
    if ~isempty(idA2)
        Cdim(1) = Asz(idA2);
    end
    if ~isempty(idB2)
        Cdim(2) = Bsz(idB2);
    end
end

% reshape matrix to tensor
C = reshape(C2,Cdim);

if ~isempty(idC) % if permutation option is given
    C = permute(C,idC);
end
% % % % Main computational part (end) % % % %

% % % option -v
if oshow
    % display result
    % text to contain leg and size information
    Astr = cell(1,numel(idA2)); 
    Bstr = cell(1,numel(idB2)); 
    for it = (1:numel(idA2))
        Astr{it} = ['A',sprintf('%i (%i)',[idA2(it) Asz(idA2(it))])]; % leg index (size)
    end
    for it = (1:numel(idB2))
        Bstr{it} = ['B',sprintf('%i (%i)',[idB2(it) Bsz(idB2(it))])]; % leg index (size)
    end
    
    if (numel(idA2) + numel(idB2)) > 1
        Cstr = [Astr,Bstr];
    else
        Cstr = cell(1,2); Cstr(:) = {'x (1)'};
        if ~isempty(idA2)
            Cstr(1) = Astr;
        end
        if ~isempty(idB2)
            Cstr(2) = Bstr;
        end
    end
    
    if ~isempty(idC) % if permutation option is given
        Cstr = Cstr(idC);
    end
    
    for it = (1:numel(Cstr))
        Cstr{it} = ['C',sprintf('%i',it),'=',Cstr{it}];
    end
    fprintf(['Result leg (size): ',strjoin(Cstr,', '),'\n']);
    
    % count elapsed time
    toc2(tobj,'-v');
end
% % % 

end