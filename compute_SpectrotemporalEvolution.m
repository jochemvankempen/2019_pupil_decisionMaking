function [spectral_t, t_crop, spectrum, spectrum_TSE_base, spectrum_asym, phase, mat_filt] = compute_SpectrotemporalEvolution(erp,bandlimits,ch,tOriginal,tt,fs,BL)
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
% Spectrotemporal Evolution a la Thut 2006

if tOriginal(1)<0
    t2crop(1) = tOriginal(1) + tt.t_crop;
else
    t2crop(1) = tOriginal(1) - tt.t_crop;
end
if tOriginal(end)<0
    t2crop(2) = tOriginal(end) + tt.t_crop;
else
    t2crop(2) = tOriginal(end) - tt.t_crop;
end

ts_crop = t2crop(1)/1000*fs:t2crop(end)/1000*fs;
t_crop = ts_crop*1000/fs;

% time
spectral_t=[]; cca=1;
for it = 1:tt.skip_step:length(t_crop)-tt.window
    spectral_t(:,cca) = mean(t_crop(it:it+tt.window-1));
    cca=cca+1;
end
convert2spectrum=0;
% filter
if sum(bandlimits == [8 13])==2 || sum(bandlimits == [20 35])==2
    % alpha_bandlimits = [8,13]; % defining the filter for alpha bandpass.
    [H,G]=butter(4,[2*(bandlimits(1)/fs) 2*(bandlimits(2)/fs)]); % alpha bandpass for 500Hz
    if 0 % visualise filter
        [A,B,C,D] = butter(4,[2*(bandlimits(1)/fs) 2*(bandlimits(2)/fs)]); % alpha bandpass for 500Hz
        d = designfilt('bandpassiir','FilterOrder',8, ...
            'HalfPowerFrequency1',bandlimits(1),'HalfPowerFrequency2',bandlimits(2), ...
            'SampleRate',fs);
        sos = ss2sos(A,B,C,D);
        fvtool(sos,d,'fs',fs);
    end
    spectrum        = zeros(size(erp,1),length(spectral_t),size(erp,3));
    phase           = zeros(size(erp));
    mat_filt        = zeros(size(erp));
    spectrum_asym   = zeros(size(erp,1),length(spectral_t),size(erp,3));
    convert2spectrum=1;
    nTrial = size(erp,3);
elseif sum(bandlimits == [0.05 4])==2
    [H,G]=butter(3,[2*(bandlimits(1)/fs) 2*(bandlimits(2)/fs)]); % pupil bandpass for 500Hz
    if 0 % visualise filter
        [A,B,C,D] = butter(3,[2*(bandlimits(1)/fs) 2*(bandlimits(2)/fs)]); % alpha bandpass for 500Hz
        d = designfilt('bandpassiir','FilterOrder',8, ...
            'HalfPowerFrequency1',bandlimits(1),'HalfPowerFrequency2',bandlimits(2), ...
            'SampleRate',fs);
        sos = ss2sos(A,B,C,D);
        fvtool(sos,d,'fs',fs);
    end
    nTrial = size(erp,2);
    spectrum        = zeros(length(t_crop),size(erp,2));
end


% Spectrotemporal Evolution a la Thut
for itrial = 1:nTrial
    % filtering
    
    if convert2spectrum
        ep_filt = filtfilt(H,G,squeeze(erp(:,:,itrial))')';
        phase(:,:,itrial) = angle(hilbert(ep_filt'))';
        mat_filt(:,:,itrial) = ep_filt;
        
        % chop off ends and rectify
        ep_filt = abs(ep_filt(:,find(tOriginal>=t_crop(1) & tOriginal<=t_crop(end))));
        % smooth
        cca=1;
        for it = 1:tt.skip_step:size(ep_filt,2)-tt.window
            spectrum(:,cca,itrial)  = mean(ep_filt(:,it:it+tt.window-1),2);
            cca=cca+1;
        end
        % asymmetry: right minus left. more positive = more right hemi alph
        spectrum_asym(ch.elec_pairs(:,2),:,itrial) = (spectrum(ch.elec_pairs(:,2),:,itrial)-spectrum(ch.elec_pairs(:,1),:,itrial))./...
            ((spectrum(ch.elec_pairs(:,2),:,itrial)+spectrum(ch.elec_pairs(:,1),:,itrial))/2);
    else
        ep_filt = filtfilt(H,G,squeeze(erp(:,itrial))')';
        spectrum(:,itrial) = ep_filt(find(tOriginal>=t_crop(1) & tOriginal<=t_crop(end)),:);
%         subplot(2,1,1), plot(tOriginal,erp(:,itrial)), subplot(2,1,2), plot(tOriginal, ep_filt(:,itrial))
%         pause
    end
end

% Baseline
if convert2spectrum
    if length(BL)==1
        baseline_spectrum = mean(spectrum(:,find(spectral_t<=BL),:),2);
    else
        baseline_spectrum = mean(spectrum(:,find(spectral_t>=BL(1) & spectral_t<=BL(2)),:),2);
    end
    spectrum_TSE_base = spectrum-repmat(baseline_spectrum,[1,size(spectrum,2),1]); % baseline full erp
else
    baseline = mean(spectrum(find(t_crop>=BL(1) & t_crop<=BL(2)),:),1);
    spectrum_TSE_base = spectrum-repmat(baseline,[size(spectrum,1)],1); % baseline full erp
end

