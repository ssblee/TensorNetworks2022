function [M,cfun,err,lid,cfun_test,err_test,lid_test] = ML_MPS_Ex (M,F,y,F_test,y_test,Nkeep,estep)
% < Description >
%
% [M,cfun,err,cfun_test,err_test] = ML_MPS_Ex (M,F,y,F_test,y_test,Nkeep,alpha)
%
% Machine learning method based on matrix product states (MPS), proposed by
% Stoudenmire2016 [E. M. Stoudenmire and D. J. Schwab, Adv. Neural Inf.
% Process. Syst. 29, 4799 (2016), or arXiv:1605.05775]. Two versions of the
% paper have minor differences; here we follow the notation of the
% published NeurIPS version.
%
% Note that there are some important technical details that are not
% discussed in Stoudenmire2016 and devised by Seung-Sup Lee to achieve
% stability and performance. Such parts in the code are denoted by the
% comments "[Unpublished; devised by S.Lee]".
%
% < Input >
% M : [cell array] MPS. Each cell contains rank-3 tensor for each chain
%       site. If given as empty (i.e. []), this function initializes the
%       MPS based on the feature vectors from the training dataset.
% F : [rank-3 tensor] Collection of feature vectors for the training
%       dataset.  F(m,:,n) is the feature vector for the n-th site (i.e.
%       pixel) and the m-th data (i.e. image).
% y : [matrix] Collection of correct decision functions for the training
%       dataset. y(m,n) is 1 if the m-th data (i.e. image) is labeled by
%       the n-th label, 0 otherwise.
% F_test : [rank-3 tensor] Collection of feature vectors for the test
%       dataset, which has a similar structure as "F" described above.
% y_test : [matrix] Collection of correct decision functions for the test
%       dataset, which has a similar structure as "y" described above.
% Nkeep : [numeric] Maximum bond dimension for the MPS.
% estep : [vector] Prefactors to the gradient of a B tensor, i.e. \eta
%       mentiond in between Eqs. (7) and (8) in Stoudenmire2016. The length
%       of estep defines the number of sweeps; there will be numel(estep)
%       pairs of sweeps and each pair consists of one left and one right
%       sweeps. Each value, estep(m), applies to the m-th pair of sweeps.
%
% < Output >
% M : [cell array] Result MPS optimized after sweeps.
% cfun : [matrix] Record of measured cost functions (see the second
%       paragraph of Sec. 4 of Stoudenmire2016) per training data at each
%       iteration. Its n-th column corresponds to the n-th sweep. The 
%       values are recorded following the linear indexing order of MATLAB.
%       That is, cfun(:) contains the values in the order of occurences.
%       Note that in Stoudenmire2016, the cost function is defined by
%       summing over the training dataset. But here, we use the *average*
%       cost function for each training data; that is, "mean" instead of
%       "sum". By doing this, the value of the cost function becomes
%       comparable when we use different sizes of datasets.
% err : [matrix] Record of classification error rates for predicting labels
%       at each iteration, for the training dataset. "err" follows the same
%       indexing as "cfun".
% lid : [vector] lid(m) is the label index [ranging from 1 to size(y,2)]
%       for the m-th training data, determined by the final MPS. Of course,
%       it can be wrong, and if wrong, it counts to the error rate "err".
% cfun_test : [matrix] Record of measured cost functions for the test
%       dataset, which has a similar structure as "cfun" described above.
% err_test : [matrix] Record of error rate of predicting labels for the
%       test dataset, which has a similar structure as "err" described
%       above.
% lid_test : [vector] Label indices for the test dataset, which has a
%       similar structure as "lid_test" described above.
%
% Written by S.Lee (Jul.17,2019)
% Updated by S.Lee (Jul.08,2020): Implemented a stable way of initializing
%       MPS, and the normalization of feature vectors.
% Updated by S.Lee (Jul.28,2020): Typo fixed.
% Updated by S.Lee (Nov.23,2022): Minor clean-up. Revised for the lecture
%       course at SNU.

tobj = tic2;

% sanity check of input
if size(F,1) ~= size(y,1)
    error('ERR: # of training data in F and y are inconsistent.');
elseif size(F_test,1) ~= size(y_test,1)
    error('ERR: # of test data in F and y are inconsistent.');
elseif size(y_test,2) ~= size(y,2)
    error('ERR: # of labels for training and test datasets are inconsistent.');
end

% result matrices
N = size(F,3); % number of sites
cfun = nan(N,numel(estep)*2); % cost function for the training dataset
err = nan(N,numel(estep)*2); % error rate for the training dataset
cfun_test = nan(N,numel(estep)*2); % cost function for the test dataset
err_test = nan(N,numel(estep)*2); % error rate for the test dataset

fprintf('Machine learning using MPS\n');
fprintf(['  Length = ',sprintf('%i',N),', # of training data = ', ...
    sprintf('%.4g',size(F,1)),', # of test data = ', ... 
    sprintf('%.4g',size(F_test,1)),'\n']);
fprintf(['  # of labels = ',sprintf('%.4g',size(y,2)),', Nkeep = ', ...
    sprintf('%.4g',Nkeep),', ',sprintf('%i',numel(estep)),' x 2 sweeps\n']);


if isempty(M)
% initialize MPS based on the training dataset, and update feature vectors
% in the effective basis
    disptime('Initialize MPS');
    M = cell(1,N); % MPS tensors
end

% Effective feature vectors in the bond space spanned by MPS. They are
% similar to "Hlr" used in the DMRG and TDVP codes, but here Flr{..}
% represent "vectors" while "Hlr" represent Hamiltonians that are matrices.
Flr = cell(1,N+2); % for training data
% Each cell element Flr{m} is a matrix. Each Flr{m}(n,:) corresponds to the
% feature vectors for the n-th data (i.e. n-th image).
Flr{1} = 1; Flr{end} = 1; % initialize the effective feature vectors at the left and right ends

for itN = (1:N)
    % contract Flr{itn} from the previous iteration and the local feature
    % vectors at the itn-th site
    T = reshape(F(:,:,itN),[size(F,1) 1 size(F,2)]).*Flr{itN};
    % T(n,:,:) is for the n-th data (i.e. n-th image); the 2nd and 3rd
    % dimensions of T are for the bond and physical spaces.
    % Here we use elementwise multiplication (.*) to distinguish the
    % contributions from different data (i.e. images)
    
    if isempty(M{itN})
        % [Unpublished; devised by S.Lee]
        % Initialize each tensor based on the training data set. We choose
        % the MPS tensors to better target the space, which is spanned by
        % the feature vectors.
        
        if itN < N
            % all the tensors except for the last one has three legs, as
            % usual: left-right-physical
            V = ML_MPS_rightSV(T(:,:),Nkeep); % find the dominant singular vectors
            M{itN} = permute(reshape(V,[size(T,2) size(T,3) size(V,2)]),[1 3 2]);

        else % itn == N
            % the last tensor has four legs: left-right-physical-label
            M{itN} = zeros(size(M{itN-1},2),1,size(F,2),size(y,2));
            for itl = (1:size(y,2))
                V = ML_MPS_rightSV(T(:,:).*y(:,itl),1); % find the most dominant singular vector, depending on the labels
                M{itN}(:,1,:,itl) = reshape(V,[size(T,2),1,size(T,3)]);
            end
        end
    end
    
    if itN < N
        % Flr{end+1} is not necessary, as it will not be used until it is
        % overwritten
        Flr{itN+1} = contract(T,3,[2 3],M{itN},3,[1 3]);
        
        % [Unpublished; devised by S.Lee] normalize feature vectors
        Flr{itN+1} = ML_MPS_normalize(Flr{itN+1});
    end
end

% construct Flr for the test data
Flr_test = cell(1,N+2);
Flr_test{1} = 1; Flr_test{end} = 1; % initialize the effective feature vectors at the left and right ends

for itN = (1:(N-2)) % skip for the last two sites, since Flr_test for them will not be used
    T = reshape(F_test(:,:,itN),[size(F_test,1) 1 size(F_test,2)]).*Flr_test{itN};
    Flr_test{itN+1} = contract(T,3,[2 3],M{itN},3,[1 3]);
    
	% [Unpublished; devised by S.Lee] normalize feature vectors
    Flr_test{itN+1} = ML_MPS_normalize(Flr_test{itN+1});
end

disptime('Start sweeping');

for its = (1:numel(estep))
    % left <- right
    for itN = (N:-1:2)
        B = contract(M{itN-1},3,2,M{itN},4,1,[1 2 4 3 5]);
        % leg order of M{itn-1}: left-right-physical
        % leg order of M{itn}: left-right-physical-label
        % leg order of B: left-physical(itn-1)-physical(itn)-right-label
        
        % for test data, compute only cost function and error rate
        [cfun_test(end+1-itN,2*its-1),err_test(end+1-itN,2*its-1)] = ...
            ML_MPS_1step (B,Flr_test{itN-1},F_test(:,:,itN-1),F_test(:,:,itN),Flr_test{itN+2},y_test);

        % compute cost function, error rate, gradient of B
        [cfun(end+1-itN,2*its-1),err(end+1-itN,2*its-1),dB] = ...
            ML_MPS_1step (B,Flr{itN-1},F(:,:,itN-1),F(:,:,itN),Flr{itN+2},y);
        
        % % % % TODO (start) % % % %
        % update B tensor
        
        % SVD and update M{itn-1}, M{itn}
        

        % update Flr{itn+1} in accordance with the updated M{itn}


        % [Unpublished; devised by S.Lee] normalize feature vectors
        Flr{itN+1} = ML_MPS_normalize(Flr{itN+1});
        
        % update Flr_test{itn+1} in accordance with the updated M{itn}


        % [Unpublished; devised by S.Lee] normalize feature vectors
        Flr_test{itN+1} = ML_MPS_normalize(Flr_test{itN+1});
        % % % % TODO (end) % % % %
    end
    
    disptime(['Sweep #',sprintf('%02i/%02i',2*its-1,2*numel(estep)), ...
        ' | left <- right,  eta = ',sprintf('%.3g',estep(its))]);
    fprintf(['  Training: cost fun = ',sprintf('%.3e',cfun(N-1,2*its-1)), ...
        ', error rate = ',sprintf('%.2f',err(N-1,2*its-1)*100),'%%\n']);
    fprintf(['     Test : cost fun = ',sprintf('%.3e',cfun_test(N-1,2*its-1)), ...
        ', error rate = ',sprintf('%.2f',err_test(N-1,2*its-1)*100),'%%\n']);

    % left -> right
    for itN = (2:N)
        B = contract(M{itN-1},4,2,M{itN},3,1,[1 2 5 4 3]);
        % leg order of M{itn-1}: left-right-physical-label
        % leg order of M{itn}: left-right-physical
        % leg order of B: left-physical(itn-1)-physical(itn)-right-label
        
        % for test data, compute only cost function and error rate
        [cfun_test(itN-1,2*its),err_test(itN-1,2*its)] = ...
            ML_MPS_1step (B,Flr_test{itN-1},F_test(:,:,itN-1),F_test(:,:,itN),Flr_test{itN+2},y_test);

        % compute cost function, error rate, gradient
        [cfun(itN-1,2*its),err(itN-1,2*its),dB] = ...
            ML_MPS_1step (B,Flr{itN-1},F(:,:,itN-1),F(:,:,itN),Flr{itN+2},y);

        % % % % TODO (start) % % % %
        % update B tensor
        
        % SVD and update M{itn-1}, M{itn}
        
        % update Flr{itn} in accordance with the updated M{itn-1}
        
        % [Unpublished; devised by S.Lee] normalize feature vectors
        Flr{itN} = ML_MPS_normalize(Flr{itN});

        % update Flr_test{itn} in accordance with the updated M{itn-1}
        
        % [Unpublished; devised by S.Lee] normalize feature vectors
        Flr_test{itN} = ML_MPS_normalize(Flr_test{itN});
        % % % % TODO (end) % % % %
    end
    
    disptime(['Sweep #',sprintf('%02i/%02i',2*its,2*numel(estep)), ...
        ' | left -> right,  eta = ',sprintf('%.3g',estep(its))]);
    fprintf(['  Training: cost fun = ',sprintf('%.3e',cfun(N-1,2*its)), ...
        ', error rate = ',sprintf('%.2f',err(N-1,2*its)*100),'%%\n']);
    fprintf(['     Test : cost fun = ',sprintf('%.3e',cfun_test(N-1,2*its)), ...
        ', error rate = ',sprintf('%.2f',err_test(N-1,2*its)*100),'%%\n']);
end

% label the images from the training and test datasets according to the
% final MPS
T = reshape(F(:,:,N),[size(F,1) 1 size(F,2)]).*Flr{N};
T = contract(T,3,[2 3],M{end},4,[1 3],[1 3 2]); % put the right leg with singleton dimension to the end
[~,lid] = max(T,[],2);

T = reshape(F_test(:,:,N),[size(F_test,1) 1 size(F_test,2)]).*Flr_test{N};
T = contract(T,3,[2 3],M{end},4,[1 3],[1 3 2]); % put the right leg with singleton dimension to the end
[~,lid_test] = max(T,[],2);

fprintf(['Final error rates = ', ...
    sprintf('%.2f',100*(1 - sum( y((1:size(y,1)).'+(lid-1)*size(y,1)) )/size(y,1))), ...
    '%% (training), ', ...
    sprintf('%.2f',100*(1 - sum( y_test((1:size(y_test,1)).'+(lid_test-1)*size(y_test,1)) )/size(y_test,1))), ...
    '%% (test)\n']);


toc2(tobj,'-v');
chkmem;

end


function [cfun,err,varargout] = ML_MPS_1step (B,F1,F2,F3,F4,y)
% < Description >
%
% [cfun,err [, dB] ] = ML_MPS_1step (B,F1,F2,F3,F4,y)
%
% < Input >
% B : [rank-5 tensor] The contraction of two local tensors of the MPS,
%       associated with the current orthogonality center.
% F1, F2, F3, F4 : [tensors] Feature vector data. F2 and F3 are local
%       feature vectors at two sites associated with the current
%       orthogonality center. F1 (F4) is the contraction of feature vectors
%       at the left (right) parts of chain with the MPS tensors; F1 and F4
%       are the feature vectors in effective basis.
% y : [matrix] Collection of correct decision functions.
%
% < Output >
% cfun, err : [numeric] Cost function per data and error rate in predicting
%       labels, respectively.               
% dB : (Optional) [rank-5 tensor] Gradient for the B tensor. It is \Delta B
%       in Eq. (7) of Stoudenmire2016.
%

% % % % TODO (start) % % % %
% % evaluate decision function f^l (x) in Eq. (6) of Stoudenmire2016

% Goal: generate matrix fx, such that fx(n,l) means f^l (x_n)

% insert one dimension (leg) to the front, which corresponds to data set indices

% % deviation of the decision function from the correct value;
% % y_n^l - f^l (x_n) in Fig. 6(d) of Stoudenmire2016

% % cost function per data set; see the second paragraph of Sec. 4 of
% % Stoudenmire2016
cfun = sum(abs(ydiff(:)).^2)/2/size(y,1);
% [Unpublished; devised by S.Lee] Normalize "cfun" with the number of data,
% so that one can plot the cost function values with the error rates on the
% same scale

if nargout > 2 %~isempty(estep)
    % % compute the gradient \Delta B in Eq. (7) of Stoudenmire2016
    
    % Goal: generate rank-5 tensor dB that corresponds to \Delta B
    
    % [Unpublished; devised by S.Lee]
    % In Eq. (7) of Stoudenmire2016, \Delta B is obtained by the sum over
    % the training dataset, which can be translated to sum(...) in MATLAB.
    % But here we use mean(...), i.e. dividing by the size of the training
    % dataset. By doing so, we can use similar choices of \eta's (estep)
    % for different training dataset sizes.
    dB = mean(dB,1);
    dB = reshape(dB,[size(dB,2) size(dB,3) size(dB,4) size(dB,5) size(dB,6)]);
    
    varargout{1} = dB;
end
% % % % TODO (end) % % % %

end

function V = ML_MPS_rightSV (M,Nkeep)
% [Unpublished; devised by S.Lee]
% Obtain the right singular vectors of a matrix M associated with largest
% singular values. When there are more than Nkeep singular vectors, the
% function chooses only the singular vectors associated with the Nkeep
% largest singular values

% as left singular vectors are not needed and the row dimension is huge, we
% first contract M to use standard eigendecomposition.
M2 = M'*M;
[V,D] = eig(M2+M2');
[~,ids] = sort(diag(D),'descend');
ids((Nkeep+1):end) = [];
V = V(:,ids);

end

function F = ML_MPS_normalize(F)
% [Unpublished; devised by S.Lee]
% Normalize each row of matrix F so that the row vector has norm 1. Within
% the function ML_MPS, each row means the feature vectors of each data set.

F = F./sqrt(sum(abs(F).^2,2));
% in case of divide-by-zero
F(isinf(F)) = 0;
F(isnan(F)) = 0;

end
