function [w,s,C,nor] = LDAbeamformer(P,dat,gamma,verbose)
%Linear-discriminant analysis (LDA) beamformer. Given a spatial pattern P 
%and sensor data DAT, a spatial filter w is constructed. 
%If more than one pattern is provided, the additional patterns are used as
%nulling constraints and an LDA nulling beamformer is constructed.
%
%Usage:
% lda_beamformer(P,DAT)
% lda_beamformer(P,DAT,gamma,verbose)
%
%Inputs:
% P: [vector or matrix] column vector specifying the spatial pattern  
%          estimated from the data. The spatial pattern will typically correspond
%          to a mean ERP peak in one experimental condition; alternatively, it can 
%          correspond to a peak in a difference-ERP that contrasts two different experimental
%          conditions.
%          P can also be a matrix of spatial patterns as columns. 
%          In this case, the additional patterns are 
%          used as nulling constraints, i.e. the LDA beamformer is set to 
%          give a zero response to the additional patterns.
%
% DAT: [matrix or FieldTrip struct] matrix with neural data (2D or 3D). 
%          Used to calculate the empirical covariance matrix 
%          and estimate the source time series.
%          Alternatively, the covariance matrix can be provided directly.
%          DAT can take three different forms:
%          (1) 2D matrix (channels x time samples) representing continuous
%          (unepoched) data. Empirical covariance is calculated on the
%          whole data.
%          (2) 3D matrix (channels x time samples x epochs) representing
%          epoched data. The empirical covariance is calculated for each
%          epoch separately and then averaged across epochs.
%          (3) Covariance matrix. If the dimension of DAT is channels x
%          channels it is assumed that a covariance matrix has been
%          provided. 
%          Instead of a matrix, a FieldTrip struct with continuous or
%          epoched data can be provided. In this case, the .trial field 
%          is used.
%
%Optional:
% gamma [double] regularization parameter. Blends between the
%                empirical covariance matrix (gamma=0, no regularization) and a
%                diagonal matrix of equal variance (gamma=1, maximum regularization). 
%                Default 0.0001.
% verbose [bool] give output on the console (default 1)
%
%Outputs:
% w     [vector] LDA beamformer / spatial filter
% C     [matrix] regularized covariance matrix
% nor   [double] normalization factor that scales the estimated time series
% s     [vector or matrix or FieldTrip struct] 
%                estimated source time series, calculated as w'*DAT. 
%                If DAT is 3D, s is a time points x number of epochs matrix
%                with the beamformed source time series for each epoch. If
%                a FieldTrip struct was given as input, s is a FieldTrip
%                struct and the source time series is given in the .trial
%                field.
%
%Remarks: 
% LDA is usually used to classify between two or more classes/experimental
%conditions. However, the LDA beamformer is an approach for estimating time
%series that does not require multiple classes. It only
%needs a spatial pattern which can either be a difference pattern between 2
%classes (p = pattern1 - pattern2) as in ordinary LDA, or just a single pattern
%(p = pattern1) if there is only 1 class. In the latter case, the task is
%equivalent to maximising the difference between this class and a
%"zero class" which has the same covariance but zero mean (p = pattern1 - 0).
%
%Just like the LCMV beamformer, the LDA beamformer is given by the solution 
%of the quadratic programming problem
%
%         arg min       w' C w 
%         subject to    w'p = 1
%
%where w is the beamformer (spatial filter), C is the covariance matrix, p 
%is the target spatial pattern, and w' is the transpose of w. In other
%words, the spatial filter is designed such that it permits the signal
%characterised by the topography p with unit gain (w'p = 1) while at the
%same time minimising the variance of the estimated source time series.
%The latter can be seen as follows. Assume X is the [time samples x channels]
%data matrix, then C = X'X is proportional to the covariance matrix.
%Furthermore, the estimated source time series vector s
%is obtained by projecting the data on the spatial filter: s = (Xw). Then
%we have
%
% w' C w = w' X'X w = (Xw)' (Xw) = s^2
%
%So the beamformer minimises the variance of the source time series. Since
%the linear constraint assures that the target signal is preserved, we
%essentially minimise the noise.
%
%Previously, beamforming has been extended (see e.g. Hui et al, 2010) to 
%include nulling constraints. This extends the optimisation problem to
%
%         arg min       w' C w 
%         subject to    w'p = 1
%                       w'n = 0
%
%The additional constraint (w'n = 0) assures that a signal with topography n
%is fully suppressed. The LDA beamformer can be extended with nulling
%constraints by adding columns with spatial topographies to the argument P
%[not mentioned in the paper]. In this case, the first column of P (P(:,1)) 
%is assumed to be the target topography p while all the other columns P(:,2:end) 
%are assumed to be noise topographies n. Nulling can be useful if eg time
%series of several ERP components are to be compared to each other and one
%wants to make sure that activity from one source does not leak
%into the estimate of the other source. 
%
%Reference:
% MS Treder, AK Porbadnigk, F Shahbazi Avarvand, K-R Müller, B Blankertz. The LDA beamformer:
% Optimal estimation of ERP source time series using linear discriminant analysis. NeuroImage.
% Volume 129, 1 April 2016, Pages 279–291.
%
%Other relevant references:
% M van Vliet, N Chumerin, S De Deyne, JR Wiersema, W Fias, G Storms, MM Van Hulle.
% Single-Trial ERP Component Analysis Using a Spatiotemporal LCMV Beamformer.
% IEEE Transactions on Biomedical Engineering, Volume 63, Issue 1, Pages 55-66.
%
% HB Hui, D Pantazis, SL Bressler, RM Leahya. Identifying True Cortical 
% Interactions in MEG using the Nulling Beamformer. Neuroimage, 
% 2010 Feb 15; 49(4): 3161. 
%
% (c) Matthias Treder 2013-2016

%% Process input arguments
if isvector(P)
  % Assure that P is a column vector
  P = P(:);
end

if ~exist('gamma','var') || isempty(gamma)
  gamma = 0.0001;
end
if ~exist('verbose','var') || isempty(verbose)
  verbose = 1;
end

% Check if we have a FieldTrip struct
isFT = isstruct(dat);
if isFT,  X= dat.trial;
else,     X= dat; clear dat;
end

if iscell(X)
    nEpochs= numel(X);
    nChan = size(X{1},1);
else
    [nChan, ~, nEpochs] = size(X);
end

if verbose, fprintf('Processing data with %d channels and %d epochs\n',nChan,nEpochs); end

%% Calculate covariance matrix
if size(P,1) ~= nChan
    error('First dimension of data [%d] does not match the number of channels inferred from P [%d]',nChan, size(P,1))
end
if iscell(X)
    if verbose, fprintf('Data provided as cell array. Calculating the covariance across the cells\n'), end
    C = zeros(nChan);
    for ii=1:nEpochs
        C = C + cov(squeeze(X{ii}'));
    end
    C = C / nEpochs;
elseif ndims(X) == 2
    if all(size(X) == size(X'))
        if verbose, fprintf('Covariance matrix has been provided directly\n'), end
        C= X;
    else
        if verbose, fprintf('Continuous data has been provided and is used for calculating the empirical covariance\n'), end
        C = cov(X');
    end
elseif ndims(X) == 3
    if verbose, fprintf('Epoched data provided. Calculating the covariance across the epochs\n'), end
    C = zeros(nChan);
    for ii=1:nEpochs
        C = C + cov(squeeze(X(:,:,ii)'));
    end
    C = C / nEpochs;
else error('Number of dimensions of X must be 2 or 3')
end

%% Regularise covariance matrix
if gamma>0
  % Define shrinkage target as diagonal matrix with with entries equal to
  % trace/(number sensors) (=mean variance) of the covariance matrix. The
  % scaling factor nu puts it on equal footing (same total variance) 
  % with the empirical covariance matrix
  nu = trace(C)/nChan;
  I  = nu*eye(nChan);

  % The regularized covariance matrix is a convex combination of the empirical
  % covariance matrix and a diagonal covariance matrix with equal total variance.
  C = C*(1-gamma) + I*gamma;
end

%% Calculate spatial filter w
% f: vector representing the linear constraint(s) 
f = zeros(size(P,2),1); 
f(1)= 1;

% nor: Normalization factor
nor= (P'* (C\P) );

% w: spatial filter
w = (C\P) * (nor\f);

%% Calculate source time series s
if nargout>1
    if isFT
        s= dat;
        clear dat
        s.label= {'LDAsource'};
        s.w= w;
        if iscell(X)
            for ii = 1:nEpochs
                s.trial{ii}= w'*X{ii};
            end
        else
            s.trial= zeros(1,size(s.trial,2),nEpochs);
            for ii = 1:nEpochs
                s.trial(1,:,ii)= w'* squeeze(X(:,:,ii));
            end
        end
    elseif iscell(X)
        s= cell(nEpochs,1);
        for ii = 1:nEpochs
            s{ii}= w'*X{ii};
        end
    elseif ndims(X)==2
        s = w'*X;
    elseif ndims(X)==3
        s = zeros(size(X,1),size(X,3));
        for ii = 1:nEpochs
            s(:,ii)= w'*squeeze(X(:,:,ii));
        end
    end
end
