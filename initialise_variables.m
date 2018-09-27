% These scripts reproduce the analysis in the paper: van Kempen et al.,
% (2018) 'Behavioural and neural signatures of perceptual evidence
% accumulation are modulated by pupil-linked arousal'. 
% 
% Many of these scripts are based on the original scripts for the paper
% Newman et al. (2017), Journal of Neuroscience.
% https://github.com/gerontium/big_dots 
%
% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject to
% the following conditions:
% 
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
% LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
% OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
% WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
% Jochem van Kempen, 2018
% Jochemvankempen@gmail.com
% https://github.com/jochemvankempen/2018_Monash
%
% -------------------------------------------------------------------------
%
% Initialise variables, this is necessary when muliple files (e.g.
% different blocks/sessions) are loaded and need to be stored in the same
% variable. Also, collect_data.m expects consistent variable sizes across
% blocks/subjects etc.

nTrFile              = length(allTrig); % n Trial per File
nTrTotal             = length(allTrig) * length(loadfilenames); % n total Trial for all Files

trialIdx            = (1:nTrTotal)';
tmpblockTrialIdx=[];
tmp = unique(blockIdx);
for iblockIdx = 1:length(unique(blockIdx))
    tmpblockTrialIdx       = [tmpblockTrialIdx 1:length(find(blockIdx==tmp(iblockIdx)))];
end
clear idx_trFile_l
idx_trFile_l = zeros(nTrFile,1);
idx_trFile_l((1:nTrFile)' + (nTrFile*(ifile-1)))=1; % get logical trial index

idx_trFile = (1:nTrFile)' + (nTrFile*(ifile-1)); % get trial idx


% trial indices
sideStimtr(idx_trFile,1)          = zeros(nTrFile,1); %side of the screen the target was presented
motiontr(idx_trFile,1)            = zeros(nTrFile,1); %side of the screen the target was presented
ititr(idx_trFile,1)               = zeros(nTrFile,1); %inter trial interval
hits(idx_trFile,1)                = zeros(nTrFile,1); %whether trial was correctly answered
misses(idx_trFile,1)              = zeros(nTrFile,1); %whether target was missed
validrlock(idx_trFile,1)          = false(nTrFile,1); %if RT was valid
corrected_drugtr(idx_trFile,1)    = zeros(nTrFile,1); %drug index, corrected for placebo-drug order (double-blind paradigm)
fileIdx(idx_trFile,1)             = zeros(nTrFile,1); %file index, used for CD paradigm, up down
subject(idx_trFile,1)             = zeros(nTrFile,1); %index of subject
blockTrialIdx(idx_trFile,1)       = tmpblockTrialIdx';

% RT
subRT(idx_trFile,1)               = zeros(nTrFile,1); %RT for each trial
subRT_log(idx_trFile,1)           = zeros(nTrFile,1); %logarithm of RT
subRT_zscore(idx_trFile,1)        = nan(nTrFile,1); %zscore of logarithm of RT
validRT(idx_trFile,1)             = false(nTrFile,1); %whether RT was valid

% subRT_window(idx_trFile,:)        = NaN(nTrFile,RTwinsize*2+1); %RT for each trial, up to 5 back
% subRTz_window(idx_trFile,:)       = NaN(nTrFile,RTwinsize*2+1); %RT for each trial, up to 5 back
% validRT_window(idx_trFile,:)      = false(nTrFile,RTwinsize*2+1); %whether RT was valid

% ERP
rawERP=0;% file gets really big if you do compute this
if rawERP
    ERP(nChan,length(t),idx_trFile)         = zeros(nChan,length(t),nTrFile); %nChan, ntimepoints, ntrial, nSub
    ERPr(nChan,length(t),idx_trFile)        = zeros(nChan,length(tr),nTrFile); %
    ERP_csd(nChan,length(t),idx_trFile)     = zeros(nChan,length(t),nTrFile); %nChan, ntimepoints, ntrial, nSub
    ERPr_csd(nChan,length(t),idx_trFile)    = zeros(nChan,length(tr),nTrFile); %
end

% N2
N2c(1:length(t),idx_trFile)                 = zeros(length(t),nTrFile); %ntimepoints, ntrial  nSub
N2i(1:length(t),idx_trFile)                 = zeros(length(t),nTrFile); %
N2c_topo                                    = zeros(nChan, nTrFile);
N2i_topo                                    = zeros(nChan, nTrFile);

% CPP
CPP(1:length(t),idx_trFile)                 = zeros(length(t),nTrFile); %ntimepoints, ntrial  nSub
CPPr(1:length(tr),idx_trFile)               = zeros(length(tr),nTrFile); %
CPP_csd(1:length(t),idx_trFile)             = zeros(length(t),nTrFile); %ntimepoints, ntrial  nSub
CPPr_csd(1:length(tr),idx_trFile)           = zeros(length(tr),nTrFile); %
CPP_topo                                    = zeros(nChan, nTrFile);

% spectral - alpha
alphaRh_preTarget(idx_trFile,1)             = zeros(nTrFile,1); %ntrial  nSub
alphaLh_preTarget(idx_trFile,1)             = zeros(nTrFile,1); %ntrial  nSub
alpha_preTarget(idx_trFile,1)               = zeros(nTrFile,1); %ntrial  nSub
alphaAsym_preTarget(idx_trFile,1)           = zeros(nTrFile,1); %ntrial  nSub

alpha_preTarget_topo(1:nChan,idx_trFile)    = zeros(nChan,nTrFile); %ntrial nChan nSub
alphaAsym_preTarget_topo(1:nChan,idx_trFile)= zeros(nChan,nTrFile); %ntrial nChan nSub

% spectral - beta
beta_postTarget_topo(1:nChan,idx_trFile)        = zeros(nChan,nTrFile); %ntrial nChan nSub
beta_preResponse_topo(1:nChan,idx_trFile)       = zeros(nChan,nTrFile); %ntrial nChan nSub
beta_base_preResponse_topo(1:nChan,idx_trFile)  = zeros(nChan,nTrFile); %ntrial nChan nSub

%pupil
lpPupil.bp.resp_locked_neg500_200(idx_trFile,1) = NaN(nTrFile,1);
lpPupil.lp.resp_locked_neg500_200(idx_trFile,1) = NaN(nTrFile,1);

% trial window
trialWindow_idx(idx_trFile,:)                   = NaN(nTrFile,RTwinsize*2+1); % for each trial, up to RTwinsize backward and forward
trialWindow_valid_neg100_RT_200(idx_trFile,:)   = false(nTrFile,RTwinsize*2+1); % whether trial was valid
trialWindow_valid_neg500_0(idx_trFile,:)        = false(nTrFile,RTwinsize*2+1); % whether trial was valid
trialWindow_valid_RT(idx_trFile,:)              = false(nTrFile,RTwinsize*2+1); % whether trial was valid



