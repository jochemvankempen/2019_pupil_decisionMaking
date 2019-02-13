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

%% plot N2c versus N2c aligned to RT modulated target onset

% The reviewers said that with this level of coherently moving dots, we
% cannot assume there is any evidence accumulation. Rather, it could be
% that target onset times are just delayed on trials with larger RT. To
% control for this, we align data to a modified target onset time, which is
% the difference between a subject's fastest RT and the RT on any given
% trial. This essentially aligns the data to response. If a signal's
% amplitude would be higher, or less variable, this would indicate better
% alignment.
% The mean of N2c is lower, and the std is higer when aligned to the new
% 'target onset' event. This indicates that target selection is not
% necessarilly delayed on trials with slower RT, or at least that there is
% no fixed amount of time between target selection (N2c) and detection
% (CPP)

clear n2c_amplitude n2cr_amplitude

nrow = 1;
ncol = 3;
subplotgap = [0.06 0.08];
subplotmargin = [0.1 0.1];

ylim2use = [-3 1];

[figHandle, fSet] = figInit('fig',1, {'height', 8; 'width', 22});



data2use = allN2c;
% data2use = allN2i;
% data2use = allCPP;
% data2use = allPupil_bp;

tt2use = t;

measure = 'mean';
% measure = 'std';

validtr2use = allvalidtr_neg100_RT_200;

tWin = [-100 800];

tIdx = find(tt2use>=tWin(1) & tt2use<=tWin(2));

allData_extract = NaN(nSub,length(tIdx));
allData_extract_mod = NaN(nSub,length(tIdx));
for isub = 1:nSub
    
    trIdx = find(allSubject==isub & validtr2use);
    
    tmpData = NaN(length(find(trIdx)), length(tIdx));
    tmpData_mod = NaN(length(find(trIdx)), length(tIdx));
    
    RTdiff = allRT(trIdx) - min(allRT(trIdx));
    for itrial = 1:length(find(trIdx))
        tmpData(itrial, :) = data2use(tIdx, trIdx(itrial));
        tmpData_mod(itrial, :) = data2use(tIdx + ceil(RTdiff(itrial)/1000*fs), trIdx(itrial));
    end 
    
    switch measure
        case 'mean'
            allData_extract(isub, :) = squeeze(mean(tmpData,1));
            allData_extract_mod(isub, :) = squeeze(mean(tmpData_mod,1));
        case 'std'
            allData_extract(isub, :) = squeeze(std(tmpData,[],1));
            allData_extract_mod(isub, :) = squeeze(std(tmpData_mod,[],1));
    end
end

axH = subtightplot(1,3,1, subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 1)
hold on
% plot(tt2use(tIdx),mean(allData_extract,1))
% plot(tt2use(tIdx),mean(allData_extract_mod,1))
h = boundedline(tt2use(tIdx),mean(allData_extract,1), std(allData_extract,[],1)/sqrt(nSub));
set(h,'LineWidth', fSet.LineWidth);

xlabel('Time from target onset (ms)', 'Fontsize', fSet.Fontsize_text)
ylabel({'Amplitude (\muV)'},'fontsize',fSet.Fontsize_text)
xlim([-100 800])
ylim(ylim2use)
[~,tIdx_amp] = min(mean(allData_extract,1));

plot( tt2use(tIdx(tIdx_amp)) + [-25 25] * (1000/fs),[-2.8 -2.8], 'k', 'LineWidth', fSet.LineWidth)
n2c_amplitude = mean(allData_extract(:, tIdx_amp + [-25:25]),2);


%% allN2cr
data2use = allN2cr;

tt2use = tr;

tWin = [-1000 100];

tIdx = find(tt2use>=tWin(1) & tt2use<=tWin(2));

allData_extract = NaN(nSub,length(tIdx));
for isub = 1:nSub
    
    trIdx = find(allSubject==isub & validtr2use);
    
    tmpData = NaN(length(find(trIdx)), length(tIdx));
    
    for itrial = 1:length(find(trIdx))
        tmpData(itrial, :) = data2use(tIdx, trIdx(itrial));
    end 
    
    switch measure
        case 'mean'
            allData_extract(isub, :) = squeeze(mean(tmpData,1));
            allData_extract_mod(isub, :) = squeeze(mean(tmpData_mod,1));
        case 'std'
            allData_extract(isub, :) = squeeze(std(tmpData,[],1));
            allData_extract_mod(isub, :) = squeeze(std(tmpData_mod,[],1));
    end
end

axH = subtightplot(1,3,2, subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 2)
hold on
% plot(tt2use(tIdx),mean(allData_extract,1))
h = boundedline(tt2use(tIdx),mean(allData_extract,1), std(allData_extract,[],1)/sqrt(nSub));
set(h,'LineWidth', fSet.LineWidth);

xlabel('Time from response (ms)', 'Fontsize', fSet.Fontsize_text)
ylabel({'Amplitude (\muV)'},'fontsize',fSet.Fontsize_text)

xlim([-1000 0])
ylim(ylim2use)

[~,tIdx_amp] = min(mean(allData_extract,1));
plot( tt2use(tIdx(tIdx_amp)) + [-25 25] * (1000/fs),[-2.8 -2.8], 'k', 'LineWidth', fSet.LineWidth)

n2cr_amplitude = mean(allData_extract(:,  tIdx_amp + [-25:25]),2);

%% 
axH = subtightplot(1,3,3, subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 3)
hold on
hbar = bar(mean([n2c_amplitude n2cr_amplitude]));
hbar.FaceColor = [0.5 0.5 0.5];

plot([1 1], mean(n2c_amplitude) + [-std(n2c_amplitude) +std(n2c_amplitude)]/sqrt(nSub), 'k', 'linewidth', fSet.LineWidth)
plot([2 2], mean(n2cr_amplitude) + [-std(n2cr_amplitude) +std(n2cr_amplitude)]/sqrt(nSub), 'k', 'linewidth', fSet.LineWidth)
% clf

set(gca,'xtick', 1:2, 'xticklabel',{'Onset', 'Response'}, 'xticklabelrotation', 25, 'Fontsize', fSet.Fontsize_text)
ylabel({'Peak amplitude (\muV)'},'fontsize',fSet.Fontsize_text)

[P,H,STATS] = signrank(n2c_amplitude, n2cr_amplitude);

sigstar([1 2], P)


saveFigName = ['N2c_onset_response'];

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], {'png','svg'})

