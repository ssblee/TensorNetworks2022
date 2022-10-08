%% MATLAB 101
% Author: <https://cqm.snu.ac.kr/ Seung-Sup Lee>
%% 
% We explain the basics of MATLAB. The strengths of the MATLAB are
%% 
% * Intuitive and efficient linear algebra operations.
% * Supports tons of mathematical functions.
% * Platform independent.
% * Many toolboxes (e.g., parallel computing, signal processing, big data, ...)
% * Graphic user interface.
% * Powerful editor, debugger, and profiler.
% * Supports compiler.
% * Last but not least: QSpace, one of the most powerful tensor network libraries 
% and used by my group, runs on MATLAB.
%% 
% You can try out the following commands in three different ways:
% 
% (1) First click a section. The selected section will be highlighted with light 
% blue borders. Then click "Live Editor" (top left), then "Run Section" (top right) 
% (see a screenshot below). It can be done by shortcut command Ctrl+Enter for 
% Windows or Command+Enter for MacOS.
% 
% 
% 
% (2) Type the commands into the "Command Window".
% 
% (3) Write a MATLAB script or function (.m file) under the |MyWork| directory. 
% You can run it (i) by typing the name of the script or function (only the file 
% name, without the |.m| extension) in the "Command Window", or (ii) by clicking 
% "Editor" then "Run" (shortcuts: F5 in Windows or Option+Command+R in MacOS) 
% (see a screenshot below).
% 
% 
% 
% *For this tutorial, there is no exercise problem. Experienced users of MATLAB 
% may skip this tutorial, and move on to the next tutorials.*
%% Initialization
% First, clear memory to avoid any possible collision.

clear % clear memory
%% Basic algebra
% Algebra for numbers. MATLAB basically uses double type variables. In practice, 
% this policy helps us to avoid data conversion mistakes (which frequently happen 
% in, e.g., C programming).

A = 1; % assign A = 1
B = 2; % assign B = 2
A+B
A-B
A*B
A/B
%% 
% Also there are some predefined constants.

i % imaginary number \sqrt{-1}
1i % the same
pi % pi
%% 
% e is not a predefined constant.

e % doesn't work
%% 
% When a live in a section contains an error, as here, subsequent lives are 
% not executed. Therefore, fix the error (here: comment it out) to be able to 
% proceed.

exp(1) % Euler constant
%% 
% For unusual calculations, MATLAB also supports |Inf| (infinity) and |NaN| 
% (not-a-number).

1/0 % + infinity
-1/0 % - infinity
0/0 % not-a-number
%% 
% To suppress displaying result, put |;| at the end of command.

A = 1; B = 2;
C = A+B  % substitute to C
C = A+B; % suppress displaying result by putting ';' to the end
%% Numerical precision
% Most of numeric variables in MATLAB have so-called double precision, unless 
% specified. It means that the numbers beyond 16 digits to the right of the leading 
% significant number are rounded. A simple example to see this is:

clear
sqrt(2)^2 - 2
%% 
% Analytically, $(\sqrt{2})^2 - 2$ should be strictly zero, but in numerics, 
% the digits below $\sim 10^{-16}$ were rounded at each of |sqrt| and |^2|, hence 
% finite value $\sim 10^{-16}$. Thus in most of numerical methods, such small 
% numbers are _de facto_ zeros. Of course, the precision is relative to the magnitude 
% of numbers. To see the precision for a numeric variable, you can use the |eps| 
% command.

eps(2)
eps(1e10)
eps([1 1e3 1e5])
%% Vectors and matrices
% We can create vectors and matrices.

clear
A = [1 2 3] % row vector
A = [1,2,3] % row vector (space and , work in the same way)
A = [1;2;3] % column vector
A = [1 2; 3 4] % matrix
%% 
% The vector whose elements constitute the arithmetic series can be generated 
% easily.

A = (1:3) % row vector, arithmetic series
A = (1:3:10) % start from 1, step size = 3, up to <= 10
A = (10:12:100) % starting term can be different
A = (1:3:2) % 1
A = (1:3:0) % empty
A = (1:-3:-10) % also negative step size possible
%% 
% Functions are available for generating commom types of vectors, matices, and 
% multi-dimensional arrays with specified sizes.

A = rand(3,2) % 3*2 matrix with random elements in interval (0,1)
A = ones(3,2) % 3*2 matrix with all ones
A = zeros(3,2) % 3*2 matrix with all zeros
A = rand(3) % 3*3 matrix with random elements
A = ones(3) % 3*3 matrix with all ones
A = zeros(3) % 3*3 matrix with all zeros
A = rand(3,2,3) % multi-dimensional array
A = eye(3) % 3*3 identity matrix
A = eye(3,4) % 3*4 matrix whose diagonal elements are 1 and others are 0
A = eye(4,2) % 4*2 matrix whose diagonal elements are 1 and others are 0
A = eye(3,3,3) % doesn't work
%% 
% You can also define the matrix of |NaN|'s and |Inf|'s.

A = nan(2,3)
A = inf(2,3)
%% Matrix operations, element-wise operations
% In MATLAB, matrix operations are intuitive (and efficient).

clear
A = rand(3,2)
B = rand(2,4)
A*B % matrix multiplication
B*A % doesn't work
A+B % doesn't work
A-B % doesn't work
%% 
% The reason why the last three commands fail is that |A| and |B| have different 
% sizes. But they do work if |B| and |A| have the same dimensions.

B = rand(3,2)
A+B
A-B
A*B % doesn't work
A.*B % element-wise multiplication
A./B % element-wise division
A/B % equivalent to mrdivide function, giving a solution C such that C*B = A
C = A/B;
C*B-A % output will be 0, up to double precision
%% Size commands
% The size of vectors, matrices, and multi-dimensional arrays can be retrieved 
% by the following functions.

clear
A = rand(3,2)
size(A) % dimensions of A
size(A,1) % 1st dimension of A
size(A,2) % 2nd dimension of A
numel(A) % total number of elements
A = rand(3,1) % vector
size(A)
A = rand(3,2,3)
size(A)
numel(A)
%% Transpose and Hermitian conjugation

clear
A = rand(3,2)+1i*rand(3,2) 
A' % Hermitian conjugate
A.' % transpose
(1i)' % complex conjugate for a number
%% Logical variables and operations
% In addition to double type and character type (shortly mentioned in the cell 
% array section), there is also logical data type.

clear
true % logical variable
false
double(true) % 1
double(false) % 0
A = true(3,2) % logical array is possbile
A'
2 > 1 % true
2 == 2 % true (== : the same)
2 ~= 2 % false (~= : not the same)
2 >= 1 % true (>= : left is larger than or the same as right)
2 <= 1 % false (<= : left is smaller than or the same as right)
0 > 1 % false
~(2 > 1) % logical NOT operation
(2 > 1) && (3 > 1) % logical AND operation
(2 > 1) && (0 > 1)
(2 > 1) || (3 > 1) % logical OR operation
(2 > 1) || (0 > 1)
%% 
% Logical variables also can consitute an array.

true(3,3)
%% 
% The array of logical variables is useful in various operations. For example, 
% we can extract the subvector or submatrix whose elements satisfy certain condition, 
% e.g., being larger than 0.5.

A = rand(1,6)
B = (A > 0.5) % logical vector
A(B) % only the elements of A larger than 0.5
C = find(B) % indices of B which is true
A(C) % the same as A(B), but A(B) is faster
A = rand(4,3)
B = (A > 0.5) % logical 4*3 matrix
A(B) % the vector of the elements of A larger than 0.5
C = find(B) % indices of B which is true (linear indexing)
A(C) % the same as A(B), but A(B) is faster
%% Accessing the elements and submatrices of matrices
% We can access the elements and submatrices as follows. Note that the indexing 
% in MATLAB starts from 1, not 0.

clear
A = rand(5,5)
A(1,3) % Element at row 1, column 3
A(1,:) % Row vector at row 1 (':' means all possible indices)
A(:,3)
A(:,:)
A(:) % column vector with all the elements of A
A(0,1) % doesn't work (indexing starts from 1 in MATLAB)
A(2:3,3:5) % submatrix at the intersection of row 2-3 and column 3-5
A(2:3,3:end) % 'end' means the last index
A(1:2:5,2:4) % intersections of rows 1, 3, 5 with columns 2 to 4
A(1:2:end,2:end) % intersections of rows 1,3,5 with columns 2 to 5
A(10,7) % doesn't work since the index is out of range
%% 
% We can also access the elements and submatrices by using logical indexing.

ok = [true false true true false];
A(:,ok) % columns 1, 3, 4
%% 
% If a logical indexing vector is shorter than the corresponding size, the parts 
% (e.g. rows, columns) from the first are considered.

A(:,[true true]) % columns 1, 2
%% 
% The position indexing (using integer indices) and the logical indexing can 
% be combined.

A(1:2,ok) % combine position and logical indexing
%% 
% In the same way, we can substitute the values.

A = rand(5,5)
A(2,3) = 100 % substitute to a single element
A(3:4,:) = 200 % substitute to two columns by the uniform value 200
A(2:4,4:end) = 2 % substitute to a submatrix by the uniform value 2
A(:) = 0 % substitute all the elements by 0
%% 
% In addition to indexing as (row index, column index), linear indexing is also 
% available. The linear index is 1 for the upper left corner elements. Then the 
% index increases from top to bottom, then from left to right.

A = rand(3,3)
A(1) % = A(1,1)
A(3) % = A(3,1)
A(4) % = A(4,1)
B = A(:) % column vector with all the elements of A
B(1)
B(3)
B(4)
%% 
% It is consistent also for multi-dimensional arrays.

A = rand(3,2,3)
A(1)
A(10)
A(1,:,(1:2))
A(1,:,[true true false]) % same as A(1,:,(1:2))
%% Reshape and permute matrices
% We can reshape matrices, which keeps the total number of elements while changing 
% |size(A)|.

clear
A = (1:9) % row vector
B = reshape(A,[3 3])
B(:) % same as A, except that B is column vector
C = permute(B,[2 1]) % permute dimensions
B.' % for matrices, the permutation is the same as transpose
%% Sum and product
% The sum and the product of vectors, matrices, and multi-dimensional array 
% can be obtained by the following.

clear
sum(1:9) % sum integers from 1 to 9
A = reshape((1:9),[3 3])
sum(A) % row vector whose elements are the sum of individual columns
sum(A,2) % sum over the 2nd dimension (along rows) -> result: column vector
sum(A,1) % the same as sum(A) for matrix A
sum(sum(A)) % the sum of all the elements
sum(A(:)) % the same
prod(A) % row vector whose elements are the product of individual column
prod(A,2)
prod(A,1)
prod(A(:))
%% Eigendecomposition
% One of the most important functions for physics is the eigendecomposition 
% (also called spectral decomposition).

clear
A = rand(3,3);
A = (A+A') % symmetrize to make a Hermitian matrix
D = eig(A) % the vector of eigenvalues
[U,D] = eig(A) % spectral decomposition A = U*D*U'
%% 
% |U| is the unitary matrix whose columns are eigenvectors, and |D| is the diagonal 
% matrix whose diagonal elements are eigenvalues. Check the accuracy of the eigendecomposition.

U*D*U' - A % should be zero up to numerical double precision ~ 1e-16
abs(U*D*U' - A) < 10*max(eps(A(:))) % error comparable with the precision of A
U'*U % left-unitarity
U*U' % right-unitarity
%% 
% Often, the zeros above and below the diagonal of |D| can be redundant. Then 
% one may use:

[U,D] = eig(A,'vector')
%% Singular value decomposition (SVD)
% SVD is a key concept in the MPS methods. If the SVD is applied to the matrix 
% whose elements are coefficient of bipartite quantum state, it provides the Schmidt 
% decomposition.

clear
A = rand(3,3)
S = svd(A) % the vector of singular values in decreasing order
[U,S,V] = svd(A) % Singular value decomposition A  = U*S*V'
%% 
% |U| is the unitary matrix whose columns are left-singular vectors, |S| is 
% the diagonal matrix whose diagonal elements are singular values, and |V| is 
% the unitary matrix whose columns are right-singular vectors.

U*S*V' - A % should be zero up to numerical double precision ~ 1e-16
U*S*V - A % non-zero
U*U' % identity; U is unitary
U'*U % identity; U is unitary
V*V' % identity; V is unitary
V'*V % identity; V is unitary
%% 
% Note that, simiarly to |eig| function without |'vector'| option, the shape 
% of |S| differs depending on how we request the output of |svd|. Contrary to 
% the eigendecomposition, SVD is applicable to non-symmetric or non-Hermitian 
% matrix, even to non-square matrix.

A = rand(4,3)
S = svd(A) % three singular values
[U,S,V] = svd(A)
%% 
% Note that the sizes of |U|, |S|, and |V| are different.

U*S*V' - A
U*U' % identity; U is unitary
U'*U % identity; U is unitary
V*V' % identity; V is unitary
V'*V % identity; V is unitary
%% 
% But we see that the last column of |U| is redundant, since it is associated 
% with zero singular value. The |'econ'| option of |svd| function makes the result 
% compact, by removing redundant singluar vectors.

[U,S,V] = svd(A,'econ') 
U*S*V' - A % should be zero up to numerical double precision ~ 1e-16
U*U' % not identity; U is left-unitary
U'*U % identity; U is left-unitary
V*V' % identity; V is unitary
V'*V % identity; V is unitary
%% 
% Consider another matrix which has different size from the above case.

A = rand(3,4);
[U,S,V] = svd(A,'econ')
U*S*V' - A % should be zero up to numerical double precision ~ 1e-16
U*U' % identity; U is unitary
U'*U % identity; U is unitary
V*V' % not identity; V is right-unitary
V'*V % identity; V is right-unitary
%% QR decomposition
% QR decomposition is used in transforming tensors into canonical forms.

clear
A = rand(4,3)
[Q,R] = qr(A) % QR decomposition A = Q*R
%% 
% |Q| is 4-by-4 unitary matrix and |R| is 4-by-3 upper triangular matrix.

Q*R - A % should be zero up to numerical double precision ~ 1e-16
Q*Q'
Q'*Q
%% 
% Also there is a similar option as |'econ'| in |svd|.

[Q,R] = qr(A,0)
Q*R - A
Q*Q'
Q'*Q
%% Cell array
% Cell is a data type which can contain other general data types.

clear
A = cell(4,1) %  [] means empty
A{1} = rand(3,2) % substitute matrix to a cell array element
A{2} = rand(4,1) % cell array elements can have different size
A{4} = 'Hello' % even different data type
%% 
% Accessing the elements or subarrays of cell array:

celldisp(A) % shows the content of cell array.
A{1} % curly bracket: read or access the "content" of cell
A(1) % round bracket: 1*1 sub-array of cell array
A(1) = rand(3,2) % doesn't work
A{1:3}
A(1:3)
A{1:3} = 1 % doesn't work
%% 
% Multiple cells can be replaced with a 1-by-1 cell.

A(1:3) = {rand(2,2)}
celldisp(A)
%% 
% Cell array can be constructed by using the following syntax.

A = {rand(3,2),[1 2 3];'Hi',[]}
celldisp(A)
A{2} = {rand(3,2),rand(2,1)}; % cell can contain another cell array
%% 
% Reshaping and permuting work the same as numerical array.

B = reshape(A,[2 2])
B = B'
B = B.'
B = permute(B,[1 2])
%% Time counters
% To measure computational cost, we need to check the elapsed time.

clear
% real time counter
tic
A = rand(100,100)*rand(100,100);
toc
% CPU time counter
cput = cputime % elapsed CPU time after the current MATLAB session started
B = rand(100,100)*rand(100,100);
cputime - cput % difference in time (in seconds)
%% 
% We see that usually the elapsed CPU time is larger than the elapsed real time, 
% since MATLAB automatically parallelizes operations! For this course, we will 
% use the customized time counters, |tic2| and |toc2|, written by Seung-Sup Lee. 
% For example:

tobj = tic2;
C = rand(100,100)*rand(100,100);
toc2(tobj,'-v');
%% Conditional operations: if and switch
% |if| statement checks whether the following expression is |true| or |false|, 
% and execute the following commands until |else| or |end| appear.

clear
A = 1;
if A > 0
    A = A+1; % will happen
end
%% 
% |elseif| is also available.

A = 1;
if A < 2
    A = A+1; % will happen
else
    A = A-1; % not happen
end

A = 3;
if A > 5
    A = A+1; % not happen
elseif A > 2
    A = 2*A; % will happen
else
    A = A-1; % not happen
end
%% 
% Note that for the above if-elseif-else, only one of commands will be executed. 
% |switch| checks the value of a variable (e.g., |A| here) and execute the commands 
% for matching case (only one case at most).

A = 3;
switch A
    case 1
        B = 1; % will happen
    case 2
        B = 0; % not happen
    otherwise % if A does not match with anyone above with 'case'
        B = 100; % not happen
end
%% For-loops
% Below we substitute to the elements to |A| with the value as the triple of 
% index.

clear
A = zeros(3,1);
for it = (1:3)
    A(it) = it*3;
end
A
%% 
% For-loops can be nested.

A = zeros(3,2);
for it1 = (1:size(A,1))
    for it2 = (1:size(A,2))
        A(it1,it2) = it1*2+it2;
    end
end
A
%% 
% *Try to avoid for-loops as possible, and use linear algebra operations instead!* 
% MATLAB consists of multiple parts: it uses LAPACK or its relatives for low-level 
% operations (e.g., linear algebra operations), and uses Java for high-level operations 
% (e.g., for-loops, drawing figures). Of course, the former is much faster, since 
% it has less computational overhead and can be better parallelized.

A = rand(10,8);
B = rand(8,10);

tobj = tic2;
C = A*B % matrix multiplication (MATLAB built-in)
toc2(tobj,'-v');

tobj = tic2;
% 'manual' implementation of matrix multiplication using for-loops
D = zeros(size(A,1),size(B,2));
for it1 = (1:size(A,1))
    for it2 = (1:size(B,2))
        for it3 = (1:size(A,2))
            D(it1,it2) = D(it1,it2) + A(it1,it3)*B(it3,it2);
        end
    end
end
toc2(tobj,'-v') % much longer time!
D - A*B % results are the same
%% Save and load
% One can save the variables into .mat file, which is the MATLAB format of saving 
% numerical data, with |save| function.

clear % clear variables
A = rand(3,4);
save('test.mat','A') % creates test.mat in the current working directory
%% 
% Beware of wrapping the names of file and variables with |' '|, to treat them 
% as char array. When variable names are not specified, |save| function saves 
% all the variables in the workspace. Also, you can specify the absolute path 
% of the .mat file.

whos('-file','test.mat') % data in 'test.mat'
%% 
% To add variables, set |'-append'| option.

B = rand(1,2);
save('test.mat','B','-append')
whos('-file','test.mat')
%% 
% If variables of size larger than 2GB need to be saved, then use |'-v7.3'| 
% option.

save('test.mat','-v7.3')
%% 
% On the other hand, one can load the variables from .mat files.

clear
load('test.mat') % load all variable
whos
clear
load('test.mat','A') % load A only
whos
%% Other functionalities
% MATLAB provides *much more* functionalities beyond what we have explained 
% above. To explore them, use the MATLAB documentation. You can see the documentation 
% page of e.g., eig, by (i) typing in the command window:

help svd
%% 
% or by (ii) select a text |eig| in the command window or the editor window 
% and press F1. Of course, as MATLAB is very popular tool, there are many useful 
% websites, blogs, books, forums, etc. If you have any question, simply search 
% it from internet!