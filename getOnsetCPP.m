function [CPP_side_onsets, starttime] = getOnsetCPP(CPP, subject_folder, t, CPP_search_t, max_search_window, win_mean_change, consecutive_windows, plot_onset,type)
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
% Compute CPP onset latency
%
% This script is almost a direct copy of
% https://github.com/gerontium/big_dots/blob/master/BigDots_ERP_analysis.m,
% line 319-376. 

% CPP : time x trial, mean across channels if there are multiple

% constrain the search window according to parameters set in setAnalysisSettings_bigDots.
CPP = squeeze(CPP(find(t>=CPP_search_t(1) & t<=CPP_search_t(2)),:));
prestim_temp = find(t<CPP_search_t(1)); % so we can add it on after getting max peak.

% we want sliding windows for each trial, create smoothed waveform.
clear win_mean win_mean_inds tstats ps
for trial = 1:size(CPP,2)
    counter = 1;
    for j = max_search_window:2:size(CPP,1)-max_search_window
        win_mean(counter,trial) = mean(CPP([j-max_search_window+1:j+max_search_window-1],trial));
        win_mean_inds(counter) = j;
        counter = counter+1;
    end
end
starttime = t(win_mean_inds(1)+length(prestim_temp));

% do t-test to zero across the smoothed trials.
for tt = 1:size(win_mean,1)
    [~,P,~,STATS] = ttest(win_mean(tt,:) + win_mean_change);
    tstats(tt) = STATS.tstat;
    ps(tt) = P;
end

%%

%DN: added this in to explicitly make sure the "consecutive_windows" number of following p-values from onset are also lower than 0.05.
allp05 = find(ps<0.05 & tstats>0);
interuptions = [1 find(diff(allp05)>1) length(allp05)];

onsetp05=[];
correction = 0;
for iInterupt = 1:length(interuptions)-1
    
    tmpIdx = (interuptions(iInterupt)+ correction):interuptions(iInterupt+1);
    tmpAllp05 = allp05(tmpIdx);
    for i = 1:length(tmpAllp05)
        if  (i+consecutive_windows-1)<=length(tmpAllp05)
            if tmpAllp05(i+consecutive_windows-1)-tmpAllp05(i)==consecutive_windows-1 %if there is at least 15 consecutive windows where p<.05
                onsetp05=tmpAllp05(i);
                break
            end
        else
%             onsetp05=allp05(i);
            break
        end
    end
    if ~isempty(onsetp05)
        break
    end
    
    correction = 1;
end

% get timepoint of min index.
if ~isempty(onsetp05)
    onset_ind       = win_mean_inds(onsetp05);
    CPP_onset_ind   = onset_ind + length(prestim_temp); % see above, this needs to be added to get the overall time with respect to t.
    CPP_side_onsets = t(CPP_onset_ind);
else % onsetp05 is empty, no significant CPP.
    disp([subject_folder,': bugger']) %AD48C has no CPP onset
    CPP_side_onsets = NaN;
end

if plot_onset
    % %     plot the smoothed waveforms, the corresponding t-tests and p-values.
    % %     Make sure the 10 (DN:30) following p-values from onset are also lower than
    % %     0.05.
    figure(1),clf
    h = suptitle([subject_folder]);
    h.Interpreter = 'none';
    subplot(3,1,1)
    title('CPP')
    hold on
    plot(win_mean_inds,mean(win_mean,2))
    
    subplot(3,1,2)
    title('t-stat')
    hold on
    plot(win_mean_inds,tstats)
    subplot(3,1,3)
    title('p-value')
    hold on
    plot(win_mean_inds,ps), hold on
    line(xlim,[0.05,0.05],'Color','k','LineWidth',1);
    if ~isempty(onsetp05)
        line([onset_ind,onset_ind],ylim,'Color','g','LineWidth',1);
    else
        line([0,0],ylim,'Color','r','LineWidth',1);
    end
    pause(0.2)
end

