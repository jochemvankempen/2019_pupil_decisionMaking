function [pupil, pupilr] = processPupilDiameter(inp_pupil, set, validtr, RT)
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
% processes pupil diameter.
% Normalizes, baseline, get response locked pupil data, get average pupil
% around RT
%
% INPUT
% inp_pupil: single trial pupil data.
%     (1,:,:) = 
%     (2,:,:) = 
%     (3,:,:) = 
%     (4,:,:) = raw pupil diameter, AREA
%     (5,:,:) = blink interpolated pupil diameter
%     (6,:,:) = interpolated and low-pass filtered pupil diameter < 4 Hz
%     (7,:,:) = interpolated and band-pass filtered pupil diameter 0.05 - 4
%     (8,:,:) = instantaneous phase
%     (9,:,:) = interpolated and low-pass filtered pupil diameter < 1 Hz
%     (10,:,:) = interpolated and band-pass filtered pupil diameter 0.05 - 1 Hz
% set: struct with fields:
%     t
%     tr
%     ts
%     trs
%     fs:
%     filter
%     BL
% validtr: idices indicating wheter the trial contains an artifact, used
%     for the normalisation
% RT: Reaction times (ms) matching the pupil trials, used for the computation of
%     the peak response and peak time.
%
% OUTPUT
% pupil:       structure with 
%                     - lp: filtered and normalised low-pass filtered pupil diameter
%                     - bp: filtered and normalised band-pass filtered pupil diameter
%                     - scaleFactor: normalisation factor
%                     - RT: average pupil diameter around RT
% pupilr:      filtered and normalised response-locked pupil diameter
%
nTrial  = size(inp_pupil,3);
tt      = set.t;
ttr     = set.tr;
ttrs    = set.trs;
fs      = set.fs;

%% pick the data you want to use
% if set.filter
    pupil.lp = squeeze(inp_pupil(6,:,:)); %low pass filtered data 4Hz
    pupil.bp = squeeze(inp_pupil(7,:,:)); %band pass filtered data 0.1-4Hz
    
    pupil.lp_1Hz = squeeze(inp_pupil(9,:,:)); %low pass filtered data 1Hz

% else
%     pupil = squeeze(inp_pupil(5,:,:)); %blink interpolated data
% end

%% normalise

pupil.scaleFactor.lp = 0;
pupil.scaleFactor.bp = 0;
pupil.scaleFactor.lp_1Hz = 0;
for itrial = 1:length(validtr)
    if validtr(itrial) == 0
        continue
    end
    
    tmp_scaleFactor_lp = max(pupil.lp(find(tt>=set.BL(1) & tt<=(RT(itrial) * 1000/fs + 200)),itrial));%scaling factor, scale by maximum value measured in any trial
    tmp_scaleFactor_bp = max(pupil.bp(find(tt>=set.BL(1) & tt<=(RT(itrial) * 1000/fs + 200)),itrial));%scaling factor, scale by maximum value measured in any trial
    tmp_scaleFactorlp_1Hz = max(pupil.lp_1Hz(find(tt>=set.BL(1) & tt<=(RT(itrial) * 1000/fs + 200)),itrial));%scaling factor, scale by maximum value measured in any trial
    if tmp_scaleFactor_lp > pupil.scaleFactor.lp
        pupil.scaleFactor.lp = tmp_scaleFactor_lp;
    end
    if tmp_scaleFactor_bp > pupil.scaleFactor.bp
        pupil.scaleFactor.bp = tmp_scaleFactor_bp;
    end
    if tmp_scaleFactorlp_1Hz > pupil.scaleFactor.lp_1Hz
        pupil.scaleFactor.lp_1Hz = tmp_scaleFactorlp_1Hz;
    end
end

pupil.lp  = pupil.lp/pupil.scaleFactor.lp; %normalise 
pupil.bp  = pupil.bp/pupil.scaleFactor.bp; %normalise 
pupil.lp_1Hz  = pupil.lp_1Hz/pupil.scaleFactor.lp_1Hz; %normalise 



%% get response locked pupil
RTs = RT/1000*fs; %RT in samples

pupilr.lp  = zeros(length(ttr),nTrial);
pupilr.bp  = zeros(length(ttr),nTrial);
for n=1:nTrial
    [~,RTsamp] = min(abs(tt*fs/1000-RTs(n))); % get the sample point of the RT.
    if RTsamp+ttrs(1) >0 & RTsamp+ttrs(end)<=length(tt) & RTs(n)>0 % the RT larger than 1st stim RT point, smaller than last RT point.
        pupilr.lp(:,n)   = pupil.lp(RTsamp+ttrs,n);
        pupilr.bp(:,n)   = pupil.bp(RTsamp+ttrs,n);
    end
end

%% pre target pupil diameter
tIdx = tt>set.BL(1) & tt<set.BL(2);
pupil.baseline.lp = squeeze(mean(pupil.lp(tIdx,:),1))'; %use non-baselined pupil diameter
pupil.baseline.bp = squeeze(mean(pupil.bp(tIdx,:),1))'; %use non-baselined pupil diameter
pupil.baseline.lp_1Hz = squeeze(mean(pupil.lp_1Hz(tIdx,:),1))'; %use non-baselined pupil diameter

clear tIdx

%% pupil dilation response
%  Pupil dilation response was computed as the difference between the peak
%  diameter recorded during the 4 s following trial onset and the preceding
%  baseline diameter (Eldar et al. 2013, nat neuro,
%  http://www.nature.com/neuro/journal/v16/n8/full/nn.3428.html#methods)

baseline_bp             = mean(pupil.bp(find(tt>=set.BL(1) & tt<=set.BL(2)),:),1);
baseline_lp             = mean(pupil.lp(find(tt>=set.BL(1) & tt<=set.BL(2)),:),1);
pupil.RT.bp.neg200_200  = zeros(nTrial,1);
pupil.RT.lp.neg200_200  = zeros(nTrial,1);

for itrial = 1:nTrial 

    RTwindow = (tt > (RT(itrial)-200) & tt < (RT(itrial)+200));
    pupil.RT.bp.neg200_200(itrial,1)   = mean(pupil.bp(RTwindow, itrial)) - baseline_bp(itrial);
    pupil.RT.lp.neg200_200(itrial,1)   = mean(pupil.lp(RTwindow, itrial)) - baseline_lp(itrial);
 
end

%% Baseline
baseline     = mean(pupil.lp(find(tt>=set.BL(1) & tt<=set.BL(2)),:),1);
pupil.lp     = pupil.lp  - repmat(baseline,[size(pupil.lp,1)],1); % baseline full erp
pupilr.lp    = pupilr.lp - repmat(baseline,[size(pupilr.lp,1)],1); % baseline full erp

baseline     = mean(pupil.bp(find(tt>=set.BL(1) & tt<=set.BL(2)),:),1);
pupil.bp     = pupil.bp  - repmat(baseline,[size(pupil.bp,1)],1); % baseline full erp
pupilr.bp    = pupilr.bp - repmat(baseline,[size(pupilr.bp,1)],1); % baseline full erp

baseline     = mean(pupil.lp_1Hz(find(tt>=set.BL(1) & tt<=set.BL(2)),:),1);
pupil.lp_1Hz = pupil.lp_1Hz  - repmat(baseline,[size(pupil.lp_1Hz,1)],1); % baseline full erp

%% zscore 

pupil.baseline_zscore.lp = NaN(size(pupil.baseline.lp));
pupil.baseline_zscore.bp = NaN(size(pupil.baseline.lp));
pupil.baseline_zscore.lp_1Hz = NaN(size(pupil.baseline.lp));
pupil.RT_zscore.lp.neg200_200 = NaN(size(pupil.baseline.lp));
pupil.RT_zscore.bp.neg200_200 = NaN(size(pupil.baseline.lp));

pupil.baseline_zscore.lp(validtr) = zscore(pupil.baseline.lp(validtr));
pupil.baseline_zscore.bp(validtr) = zscore(pupil.baseline.bp(validtr));
pupil.baseline_zscore.lp_1Hz(validtr) = zscore(pupil.baseline.lp_1Hz(validtr));
pupil.RT_zscore.lp.neg200_200(validtr) = zscore(pupil.RT.lp.neg200_200(validtr));
pupil.RT_zscore.bp.neg200_200(validtr) = zscore(pupil.RT.bp.neg200_200(validtr));

%% get phase of Pupil

pupil.baselinephase.bp  = squeeze(angle(mean(exp(1i*(inp_pupil(8,find(tt>=set.BL(1) & tt<=set.BL(2)),:)))))); %phase pupil

