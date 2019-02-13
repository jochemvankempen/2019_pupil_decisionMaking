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
% setAnalysisSettings_bigDots_st
% this script sets the settings that are used both for preprocessing of the
% data and the further analysis. Not only used for bigDots dataset, also
% clonidine.
 
doCSD=0;% Use/compute Current Source Density transformed erp?

targcodes = [101:106]; % triggers in the EEG data, indicating start of trial
side_tags = {'Left','Right'}; %left and right targets
motion_tags = {'Down','Down'};

%% Time

% set limit to reaction time
rtlim       = [0.150 1.700]; % newman 2017. In original script default_response_time   = 1.700;

fs      = 500; % new sample rate
trs     = [-1.000*fs:fs*1.000];% resp-locked erps
ts      = -0.800*fs:(rtlim(2)+0.1+trs(end)/fs)*fs; % stim-locked erps, samples

% take longer epoch for pupil
trs_pupil     = [-0.900*fs:fs*3.000];% resp-locked pupil
ts_pupil      = -0.800*fs:(rtlim(2)+0.1+trs_pupil(end)/fs)*fs; % stim-locked pupil, samples

if abs(trs(1) - ts(1)) > (rtlim(1) * 1000)
    error(['The extraction of response-locked EEG cannot extract trials with RT shorter than ' num2str(abs(trs(1) - ts(1)))]);
end

t       = ts*1000/fs; %timepoints
tr      = trs*1000/fs;

t_pupil       = ts_pupil*1000/fs; %timepoints
tr_pupil      = trs_pupil*1000/fs;

BL_time     = [-100 0]; % baseline interval in ms
BL_spectrum = [-100 0]; % baseline SPG

% times used for artifact rejection (in samples).
% Note that default is that 2nd argument is added to RT. If times are added
% for which this is not desired, this needs to be specifically mentioned in
% batch_preprocess.m
tt.neg500_0         = [-0.5 0]   * fs;
tt.neg100_0         = [-0.1 0]   * fs;
tt.neg100_RT_200    = [-0.1 0.2] * fs; %2nd argument will be added to the RT
tt.neg100_RT_1000   = [-0.1 1.0] * fs; %2nd argument will be added to the RT
tt.neg100_RT_1500   = [-0.1 1.5] * fs; %2nd argument will be added to the RT
tt.neg100_RT_2000   = [-0.1 2.0] * fs; %2nd argument will be added to the RT
tt.neg100_RT_2500   = [-0.1 2.5] * fs; %2nd argument will be added to the RT
tt.neg100_RT_3000   = [-0.1 3.0] * fs; %2nd argument will be added to the RT
tt.neg500_RT_200    = [-0.5 0.2] * fs; %2nd argument will be added to the RT

artifact_times = fields(tt);

%% Filters, preprocessing

LPFcutoff_35Hz  = 35;      % Low Pass Filter cutoff
LPFcutoff_8Hz   = 8;       % Low Pass Filter cutoff

HPFcutoff = 0.1;       % High Pass Filter cutoff

LPF = 1; % 1 = low-pass filter the data, 0=don't.
HPF = 1;

%% CPP

% Define CPP onset search window
CPP_search_t = [50,1000];

% Size of sliding window. This isub in fact 1/4 of the search window in ms.
% So 25 isub 100ms. (25 samples x 2ms either side of a particular sample).
max_search_window = 25;

% CPP window in samples
CPP_search_ts  = [find(t==CPP_search_t(1)),find(t==CPP_search_t(2))];
ch_CPP_str = {'Pz'};
nChanCPP = length(ch_CPP_str);

%% N2
ch_N2_str = {'P8' ; 'P7'}; %right hemisphere first, if left targets are first (isideStim==1)
nChanN2 = 1;

%% spectral

tSpectral.window    = 50; % in samples.
tSpectral.skip_step = tSpectral.window/2;
tSpectral.t_crop    = 200; %cuts time with x ms

% for chronux toolbox
SPG.TW                      = [2^8]; %nSamples
SPG.TW_256                  = [2^7]; %nSamples. Check whether shorter time window makes a difference for visual signals.
SPG.stepSize                = 25; %samples
SPG.params.tapers           = [1 1];
SPG.params.Fs               = fs;
SPG.params.fpass            = [0 35]; % online filter during recordings
movingwin                   = [SPG.TW/fs SPG.stepSize/fs];
movingwin_256               = [SPG.TW_256/fs SPG.stepSize/fs];

%% alpha
t_preTarget_alpha   = [-500 0];

bandlimits_alpha    = [8, 13]; %alpha

mean_ch_alpha = {...
    'O1','PO3','PO7';...
    'O2','PO4','PO8'};
nSideAlpha = 2;

nch_alpha = length(mean_ch_alpha);
            
%% beta
t_postTarget_beta   = [300 600];

beta_power = 2; % 1 = TSE, 2 = STFT

t_preResponse_beta      = [-100 0];
% t_preResponse_beta2     = [-130 -70];
t_preTarget_beta        = [-100 0];
bandlimits_beta         = [20, 35]; %beta

ch_beta_str = {'C3'};

% STFT parameters for beta
STFT_time=[];
no_of_cycles = 8;
% % get a broad range, e.g. beta, 20 to 35Hz
fs_STFT = [20:35]; % stftlen_STFT = 140; % 8 cycles = 280 ms (at mid frequency, i.e. 28 Hz for beta), = 140 samples. 
stftlen_STFT = round((1000/round(median(fs_STFT))*no_of_cycles)/2);

% % get a specific frequency for SSVEP
% fs_STFT = [25]; % or if you want a particular SSVEP frequency
% stftlen_STFT = round((1000/fs_STFT*no_of_cycles)/2);
% for SSVEP frequency make sure it's EXACTLY a particular number of cycles of the frequency. 
% check freq_temp_STFT to make sure SSVEP frequency falls on the range

skip_step = 20;
cc=1;
for tt_beta = 1:skip_step:length(ts)-(stftlen_STFT)
    tf = tt_beta:tt_beta+stftlen_STFT-1;
    nfft = length(tf);
    freq_temp_STFT = (0:ceil((nfft+1)/2)-1)*fs/nfft;
    STFT_time(cc) = mean(t(tf));
    cc=cc+1;
end

%% pupil

pupilSet.BL     = [-100 0]; %pupil baseline
pupilSet.t      = t_pupil;
pupilSet.ts     = ts_pupil;
pupilSet.tr     = tr_pupil;
pupilSet.trs    = trs_pupil;
pupilSet.fs     = fs;

% GLM
pupilSet.prestim_range = [-800 0];  % pre-stimulus range (in ms), here set to first time of extracted epoch
pupilSet.zscore = 0;
pupilSet.meanCenter_predictors = 0;
pupilSet.orthogonalise_predictors = 0;
pupilSet.normalise_pupil = 0;% either 1 == normalise (essentially zscore), or 2 == demean

pupilSet.bin2use = 'RT';
% pupilSet.bin2use = 'baselinePupil';
pupilSet.nbin = [2 3 5];

