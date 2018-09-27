function [CPPr_slope] = getSlopeCPP(CPPr, nside, subject_folder, tr, side_tags, plot_slope,twin,type)
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
% Calculate slope of CPP

%% Extract Response Locked CPP slope"
%CPP build-up defined as the slope of a straight line fitted to the
%response-locked waveform at during "slope_timeframe_index" defined for
%each participant here:
clear slope_timeframe_index 

if size(CPPr,1)>nside %switch to row vector if there is only 1 side
    CPPr = CPPr';
end

if type==1
    % as in Newman (2017) J Neuro
    if nside>1
        slope_timeframe_index(2)=find(mean(CPPr,1)==max(mean(CPPr(:,find(tr<0 & tr>-400)),1)),1,'last');%max amplitude index
    else
        slope_timeframe_index(2)=find(CPPr==max(CPPr(find(tr<0 & tr>-400))),1,'last');%max amplitude index
    end
    slope_timeframe_index(1)=slope_timeframe_index(2)-50;%subtract 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    
elseif type==2
    % get time index
    slope_timeframe_index(1)=find(tr==twin(1));%
    slope_timeframe_index(2)=find(tr==twin(2));%
%     slope_timeframe_index(1)=slope_timeframe_index(2)-(diff(twin)/2);%subtract n samples from slope_timeframe_index(2) index to form slope_timeframe window, divide by 2 because fs = 500
end

%Now find and save CPPr slope
for side = 1:nside
    coef = polyfit(tr(slope_timeframe_index(1):slope_timeframe_index(2)),(CPPr(side,slope_timeframe_index(1):slope_timeframe_index(2))),1);% coef gives 2 coefficients fitting r = slope * x + intercept
    CPPr_slope(side)=coef(1);
end

%%%Plot each individual participant's CPPr_slope with time-window varying
%%%per participant
if plot_slope
    clear h
    figure(1),clf
    for side = 1:nside
        h(side) = plot(tr,CPPr(side,:),'LineWidth',3,'LineStyle','-');hold on
        coef = polyfit(tr(slope_timeframe_index(1):slope_timeframe_index(2)),(CPPr(side,slope_timeframe_index(1):slope_timeframe_index(2))),1);% coef gives 2 coefficients fitting r = slope * x + intercept
        CPP_slope(side)=coef(1);
        r = coef(1) .* tr(slope_timeframe_index(1):slope_timeframe_index(2)) + coef(2); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
        plot(tr(slope_timeframe_index(1):slope_timeframe_index(2)), r,'Linewidth',2, 'LineStyle', ':');
    end
    
    set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title([subject_folder, ' CPP (resp-locked) by Hemifield'],'interpret','none')
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    shg
%     legend(h,side_tags,'FontSize',16,'Location','NorthWest');
    keyboard
end