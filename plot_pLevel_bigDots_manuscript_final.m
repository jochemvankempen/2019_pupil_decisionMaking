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
% produce the plots from manuscript

%% select which data to plot


x = 1:nbin2use;

idx_BL_lp       = strcmpi(allBin2use, 'pupil_lp_baseline_regress_iti_side');
% idx_BL_lp       = strcmpi(allBin2use, 'pupil_lp_baselineSlope_regress_iti_side');
% idx_BL_lp       = strcmpi(allBin2use, 'pupil_lp_baselineDiff_regress_iti_side');

% idx_BL_lp_1Hz   = strcmpi(allBin2use, 'pupil_lp_1Hz_baseline_regress_iti_side');

idx_BL_bp       = strcmpi(allBin2use, 'pupil_bp_baseline_regress_iti_side');
idx_BL_bp       = strcmpi(allBin2use, 'pupil_bp_baselinePhase_regress_iti_side');
% idx_BL_bp       = strcmpi(allBin2use, 'pupil_bp_baselinePhase_regress_iti_side_fix');

% idx_resp        = strcmpi(allBin2use, 'pupil_bp_RT_neg200_200_regress_bl_iti_side');
% idx_resp        = strcmpi(allBin2use, 'pupil_bp_RT_neg200_200_regress_bl_iti_side_RT');
% 
% idx_resp        = strcmpi(allBin2use, 'pupil_bp_average_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side');
% idx_resp        = strcmpi(allBin2use, 'pupil_bp_slope_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side');
% idx_resp        = strcmpi(allBin2use, 'pupil_bp_diff_maxDiff_pupilIRF_200_200_regress_bl_iti_side');

% idx_resp        = strcmpi(allBin2use, 'pupil_bp_linearProjection_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side');
% idx_resp        = strcmpi(allBin2use, 'pupil_bp_linearProjection_maxDiff_pupilIRF_0_200_regress_bl_iti_side');
% 

% idx_resp        = strcmpi(allBin2use, 'GLM_pupil_StimResp_stim_regress_bl_iti_side');
% idx_resp        = strcmpi(allBin2use, 'GLM_pupil_Ramp_stim_regress_iti_side');
idx_resp        = strcmpi(allBin2use, 'GLM_pupil_Ramp_stim_regress_bl_iti_side');
idx_resp        = strcmpi(allBin2use, 'GLM_pupil_Ramp_stim_regress_blPhase_iti_side');
idx_resp        = strcmpi(allBin2use, 'GLM_pupil_Ramp_stim_regress_bl_blPhase_iti_side');
% idx_resp        = strcmpi(allBin2use, 'GLM_pupil_Ramp_sust_regress_iti_side');
% idx_resp        = strcmpi(allBin2use, 'GLM_pupil_Ramp_sust_regress_bl_iti_side');

% idx_resp        = strcmpi(allBin2use, 'GLM_pupil_Ramp_stim_VIF5_regress_bl_iti_side');


idx_alpha       = strcmpi(allBin2use, 'pretarget_alpha');

idx2plot   = idx_BL_bp;
idx2plot   = idx_BL_lp;
% idx2plot   = idx_BL_lp_1Hz;
idx2plot   = idx_resp;

% plot Settings
chanlocs = readlocs('cap64.loc'); %biosemi
[~,plot_chans, exclude_chans] = getChannelName;

figureFileType = {'png','svg'};
clear cfg
cfg.t = t;
cfg.tr = tr;
cfg.t_pupil = t_pupil;
cfg.tr_pupil = tr_pupil;
plotMarker = {'o','s','d','^','h'};


[~,idx] = grep(allBin2use(idx2plot),'GLM_pupil');
if idx == 0
    filename_R = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{idx2plot} '_' bintype fileExt_preprocess  fileExt_CDT   ];
else
    filename_R = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{idx2plot} '_' bintype fileExt_preprocess  fileExt_CDT fileExt_GLM  ];
end


% betaString = [char(946) '_{1} = ' num2str(stats.B, '%1.2f')];
% ylabel({['Amplitude (' char(956) 'V)']},'fontsize',fSet.Fontsize_text)



%% Figure 1.

nrow = 1;
ncol = 3;
subplotgap = [0.06 0.08];
subplotmargin = [];

[figHandle, fSet] = figInit('fig',1, {'height', 10; 'width', 22});

set(gca,'color','none','XColor','k','YColor','k')
set(gcf,'PaperPositionMode','auto')
% set(gca,'fontsize',fSet.plFontsize)

switch allBin2use{idx2plot}
    case {'pupil_lp_baseline_regress_iti_side','pupil_bp_baseline_regress_iti_side'}
        text_y = 0.86;
        panelLabel = 'B';
    case {'GLM_pupil_Ramp_stim_regress_bl_blPhase_iti_side'}
        text_y = 0.835;
        panelLabel = 'C';
    otherwise
        text_y = 0.855;
%         text_y = 0.55;
        panelLabel = 'C';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% paradigm figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subtightplot(nrow, ncol, 1, subplotgap, subplotmargin, subplotmargin);
% import paradimg png
fig2import = ['C:\Jochem\Dropbox\Monash\bigDots_st\figures\paradigm_test.png'];
paradigm = imread(fig2import);

image(paradigm)
axis off
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pupil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, 2, subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, panelLabel)

[~,idx1] = grep(allBin2use(idx2plot),'pupil_lp_1Hz');
[~,idx2] = grep(allBin2use(idx2plot),'pupil_lp');
if idx2 == 0
    fprintf('Plotting band-pass filtered pupil diameter\n')
    pupil2plot = pupil_bp;
    cfg.ylim.Pupil = [-0.2 0.3];
elseif idx2==1 && idx1==0
    fprintf('Plotting low-pass filtered pupil diameter\n')
    pupil2plot = pupil_lp;
    cfg.ylim.Pupil = [-0.05 0.09];
elseif idx2==1 && idx1==1
    fprintf('Plotting low-pass filtered (1Hz) pupil diameter\n')
    pupil2plot = pupil_lp_1Hz;
    cfg.ylim.Pupil = [-0.05 0.09];
end

cfg.xlim.Pupil     = [-150 1500];

hold on

clear LEGEND
for ibin = 1:nbin2use
    plot(cfg.t_pupil(1:size(pupil2plot{1},4)),squeeze(mean(squeeze(pupil2plot{idx2plot}(:,ibin,:,:)),1)),'linewidth',fSet.LineWidth_in,'Color',fSet.colors(ibin,:));
    LEGEND{ibin} = ['Bin ' num2str(ibin)];
end
[hLeg,icons,plots] = legend(LEGEND,'location','northwest','fontsize',fSet.Fontsize_text , 'box','off');
% hLeg.Position = hLeg.Position .* [0.95 0.87 1 1];
hLeg.Position = hLeg.Position .* [0.95 1.1 1 1];


XDAT = icons(6).XData;
icons(6).LineStyle = '-';
icons(6).LineWidth = fSet.LineWidth;
icons(6).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
icons(8).LineStyle = '-';
icons(8).LineWidth = fSet.LineWidth;
icons(8).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
icons(10).LineStyle = '-';
icons(10).LineWidth = fSet.LineWidth;
icons(10).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
icons(12).LineStyle = '-';
icons(12).LineWidth = fSet.LineWidth;
icons(12).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
icons(14).LineStyle = '-';
icons(14).LineWidth = fSet.LineWidth;
icons(14).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];


xlim(cfg.xlim.Pupil)
ylim(cfg.ylim.Pupil)

for ibin = 1:nbin2use
    A = boundedline(cfg.t_pupil(1:size(pupil2plot{1},4)),squeeze(mean(pupil2plot{idx2plot}(:,ibin,:,:),1)),squeeze(std(pupil2plot{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(A,'linewidth',fSet.LineWidth_in)
end
plot([0 0], cfg.ylim.Pupil ,'k','linewidth',1)

% ylim([-0.2 0.4])
xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text);
ylabel('Normalized pupil diameter','fontsize',fSet.Fontsize_text);
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, 3, subplotgap, subplotmargin, subplotmargin);
figInit('ax');
% plot_subplot_label(axH, fSet, 'E')

RT = RT;

cfg.ylim_RT     = [470 610];
% cfg.ylim_RT     = [550 680];

% cfg.ylim_RT     = [520 595];
% cfg.ylim_RT     = [440 630];
% cfg.ylim_RT     = [200 700];
% cfg.ylim_RT     = [800 1200];
cfg.ylim_RT_CV  = [0.18 0.310];
cfg.ylim_RT_CV  = [0.175 0.29];
% cfg.ylim_RT_CV  = [0.12 0.26];
% cfg.ylim_RT_CV  = [0.07 0.30];
% cfg.ylim_RT_CV  = [0.13 0.32];
cfg.ylabel.RT   = 'RT (ms)';

bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'RT',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');

switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        %         model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
%         if stats.U
%             betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
%         else
%             betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
%         end
%         %         model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end


ax1 = gca;
hold on
clear hPlotMark

if ~isempty(mfit)
    % plot model fit
    h = plot(mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', fSet.colors(1,:));
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
% plot se
for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT{idx2plot}(:,ibin))-std(RT{idx2plot}(:,ibin))/sqrt(nSub) ...
        mean(RT{idx2plot}(:,ibin))+std(RT{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.LineWidth_in);
end
% plot plotMarker
hPlotMark(ibinType) = plot(x, mean(RT{idx2plot}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(1,:),'MarkerFaceColor', fSet.colors(1,:), 'MarkerEdgeColor', 'k');

ylim([cfg.ylim_RT])
h = ylabel(cfg.ylabel.RT,'fontsize',fSet.Fontsize_text);
h.Color = fSet.colors(1,:);
set(gca,'xtick',x)
xlabel('Pupil Bin','fontsize',fSet.Fontsize_text);
axis square

h = text(0.7, mean(RT{idx2plot}(:,1))*text_y, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in, 'color', fSet.colors(1,:));

% text for RT_CV
[stats, mfit] = loadStatsR(filename_R, 'RT_CV',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        %         model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
%         if stats.U
%             betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
%         else
%             betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
%         end
    otherwise
        betaString = [];
end

text(3.3, mean(RT{idx2plot}(:,1))*text_y, {betaString, pstring_chi }, 'fontsize',fSet.Fontsize_text_in, 'color', fSet.colors(2,:))

% plot RT_CV
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k','XTickLabel',[]);
figInit('ax');
ax2.YLim = [cfg.ylim_RT_CV];
ax2.XColor = 'none';
linkaxes([ax1 ax2],'x');
hold on

if ~isempty(mfit)
    h = plot(mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', fSet.colors(2,:));
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT_CV{idx2plot}(:,ibin))-std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ...
        mean(RT_CV{idx2plot}(:,ibin))+std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.LineWidth_in);
end

hPlotMark(ibinType) = plot(x, mean(RT_CV{idx2plot}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(2,:),'MarkerFaceColor', fSet.colors(2,:), 'MarkerEdgeColor', 'k');
h = ylabel('RT CV','fontsize',fSet.Fontsize_text);
h.Color = fSet.colors(2,:);
axis square
subplot_ax1 = get(gca);

xlim([0.5 nbin2use+0.5])

saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_RT'];

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CPP ALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setAlpha = 1;
subplotgap = [0.08 0.07];
subplotmargin = [0.1 0.1];

[figHandle, fSet] = figInit('fig', 2, {'width',22;'height', 28});

nsubplot = 3; % onset, response, ITPC
nrow = 3;
ncol = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% stim locked/ onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, [1:6], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'A')
hold on

CPP2plot        = squeeze(CPP{idx2plot});
CPP_onset2plot  = squeeze(CPP_onset{idx2plot});
y_onset = [6:-0.45:4.2];
y_onset = [6.5:-0.6:3.8];
% y_onset = [23:-1.5:17]+0.7;


% clear h2
tmpCPP = nanmean(CPP_onset2plot);
for ibin = 1:nbin2use
    if setAlpha
        plot([tmpCPP(ibin) tmpCPP(ibin)], [-2 25],'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
    else
        plot([tmpCPP(ibin) tmpCPP(ibin)], [-2 25],'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) ]);
    end
end

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'CPP_onset',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

% plot model fit
if ~isempty(mfit)
    h = plot(mfit.mean,y_onset, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], [y_onset flip(y_onset)] , h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
clear h h2 LEGEND
for ibin = 1:nbin2use
    % plot CPP over time
    if setAlpha
        hLine(ibin) = boundedline(cfg.t, squeeze(mean(CPP2plot(:,ibin, :),1)), squeeze(std(CPP2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    else
        hLine(ibin) = boundedline(cfg.t, squeeze(mean(CPP2plot(:,ibin, :),1)), squeeze(std(CPP2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:));
    end
    set(hLine(ibin),'linewidth',fSet.LineWidth)
    
    % plot SE lines on plotMarker
    tmpstd = nanstd(CPP_onset2plot(:,ibin))/sqrt(nSub);
    plot([tmpCPP(ibin)-tmpstd tmpCPP(ibin)+tmpstd], [y_onset(ibin) y_onset(ibin)],'linewidth',fSet.LineWidth_in,'color','k');
    
    % plot plotMarker
    h2(ibin) = plot(tmpCPP(ibin), y_onset(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    
end
[hLeg,icons,plots] = legend([h2], LEGEND,'location','southeast','fontsize',fSet.Fontsize_text, 'box','off');

Pos = hLeg.Position;
hLeg.Position(1) = Pos(1)*1.05;
hLeg.Position(2) = Pos(2)*1.01;

icons(1).FontSize = fSet.Fontsize_text;
icons(2).FontSize = fSet.Fontsize_text;
icons(3).FontSize = fSet.Fontsize_text;
icons(4).FontSize = fSet.Fontsize_text;
icons(5).FontSize = fSet.Fontsize_text;

icons(6).LineStyle = '-';
icons(6).LineWidth = fSet.LineWidth;
icons(8).LineStyle = '-';
icons(8).LineWidth = fSet.LineWidth;
icons(10).LineStyle = '-';
icons(10).LineWidth = fSet.LineWidth;
icons(12).LineStyle = '-';
icons(12).LineWidth = fSet.LineWidth;
icons(14).LineStyle = '-';
icons(14).LineWidth = fSet.LineWidth;

plot([0 0], [-5 30] ,'-k','linewidth',1)
plot([-100 650], [0 0] ,'-k','linewidth',1)
xlim([-50 650])
ylim([-0.6 7])
% ylim([-2 25])

xBox = [min(tmpCPP)-min(tmpCPP)*0.05 max(tmpCPP)+min(tmpCPP)*0.05];
yBox = [min(y_onset)-min(y_onset)*0.08 max(y_onset)+min(y_onset)*0.08];
h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

text(max(xBox)*1.02, max(yBox)*0.95, [{betaString, pstring_chi}], 'fontsize',fSet.Fontsize_text_in)

xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text)
% ylabel({'Amplitude (\muV)'},'fontsize',fSet.Fontsize_text)
ylabel({['Amplitude (' char(956) 'V)']},'fontsize',fSet.Fontsize_text)


%%%% inset with topoplot
ax1 = axes('position',[0.1 0.76 0.4 0.13]) ; % inset
maplimits = [min(min(squeeze(mean(mean(CPP_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(CPP_topo{idx2plot},1),2))))];
topoplot(squeeze(mean(mean(CPP_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
colormap(ax1,'Jet')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% resp locked/ threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, [7:12], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'B')
hold on

CPP2plot        = squeeze(CPPr_csd{idx2plot});
CPP_amp2plot    = squeeze(CPPr_csd_amplitude{idx2plot});
CPP_slope2plot  = squeeze(CPPr_csd_slope{idx2plot});
modelLabel = {'CPPr_csd_amplitude'; 'CPP_csd_slope2'};
YLIM = [-2 30];
YLABEL = {'Amplitude (\muV/m^{2})'};
y_lines = [1.3 0.6];

% 
% CPP2plot        = squeeze(CPPr{idx2plot});
% CPP_amp2plot    = squeeze(CPPr_amplitude{idx2plot});
% CPP_slope2plot  = squeeze(CPPr_slope{idx2plot});
% modelLabel = 'CPPr_amplitude';
% YLIM = [-0.5 8];
% YLABEL = {'Amplitude (\muV)'};
% y_lines = [.1 .2];


% location to plot plotMarker CPP threshold
x_amp = [-320:15:-250]+140;

t2test_amp      = [-50 50]; % time threshold is computed on
t2test_slope    = [-250 -50]; % time slope is computed on

% plot lines indicating threshold
tmpCPP = mean(CPP_amp2plot);
for ibin = 1:nbin2use
    if setAlpha
        plot([max(x_amp)+min(x_amp)*0.38 100], [tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
    else
        plot([max(x_amp)+min(x_amp)*0.38 100], [tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:)]);
    end
end

% get stats
[stats, mfit] = loadStatsR(filename_R, modelLabel{1},nSub,nbin2use);

[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    h = plot(x_amp, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([x_amp flip(x_amp)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

% plot times amplitude and slope are computed over
plot([t2test_amp(1) t2test_amp(2)], [y_lines(1) y_lines(1)], 'k', 'linewidth', 4)
plot([t2test_slope(1) t2test_slope(2)], [y_lines(2) y_lines(2)], 'color',[0.5 0.5 0.5], 'linewidth', 4)

clear h h2 LEGEND
for ibin = 1:nbin2use
    % plot CPP timecourse
    if setAlpha
        h(ibin) = boundedline(cfg.tr, squeeze(mean(CPP2plot(:,ibin, :),1)), squeeze(std(CPP2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
        %     h(ibin) = plot(cfg.tr,squeeze(mean(CPP2plot(:,ibin, :),1)),'color',fSet.colors(ibin,:));
    else
        h(ibin) = boundedline(cfg.tr, squeeze(mean(CPP2plot(:,ibin, :),1)), squeeze(std(CPP2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:));
    end
    set(h(ibin),'linewidth',fSet.LineWidth)
    
    % plot plotMarker and se
    tmpstd = std(CPP_amp2plot(:,ibin))/sqrt(nSub);
    plot([x_amp(ibin) x_amp(ibin)], [tmpCPP(ibin)-tmpstd tmpCPP(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    h2(ibin) = plot(x_amp(ibin), tmpCPP(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end

plot([0 0], [YLIM] ,'-k','linewidth',1)
xlim([-350 60])
ylim(YLIM)

yBox = [min(tmpCPP)-min(tmpCPP)*0.1 max(tmpCPP)+min(tmpCPP)*0.1];
xBox = [min(x_amp)-min(x_amp)*0.38 max(x_amp)+min(x_amp)*0.38];

h = rectangle('Position',[xBox(2) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';


text(min(xBox)*0.99, max(yBox)*1.12, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

xlabel('Time from response (ms)','fontsize',fSet.Fontsize_text)
ylabel(YLABEL,'fontsize',fSet.Fontsize_text)
%     ylabel({'CSD ()'},'fontsize',fSet.Fontsize_text)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% resp locked/ slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes('position',[0.17 0.48 0.17 0.09]) ; % inset
figInit('ax',[],{'fontsize',10});
set(gca,'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4])
hold on

tmpCPP_slope = mean(CPP_slope2plot);
% plot(1:5, tmpCPP_slope,'linewidth',4,'color',[0.5 0.5 0.5])

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, modelLabel{2},nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
text(min(xBox)*0.95, max(yBox)*1.11, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

% plot model fit
if ~isempty(mfit)
    h = plot(1:5, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
for ibin = 1:nbin2use
    % se
    tmpCPP_slope_std = squeeze(std(CPP_slope2plot(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpCPP_slope(ibin)-tmpCPP_slope_std tmpCPP_slope(ibin)+tmpCPP_slope_std],'linewidth',fSet.LineWidth_in,'color','k');
    % plotMarker
    h2(ibin) = plot((ibin), tmpCPP_slope(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize, 'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0.5 5.5])
ylim([min(tmpCPP_slope)*0.75 max(tmpCPP_slope)*1.15])

text(1, max(tmpCPP_slope)*1.3, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)


set(gca,'xtick',[])
xlabel('Pupil bin','color',[0.4 0.4 0.4],'fontsize',fSet.Fontsize_text)
ylabel('Build-up rate','color',[0.4 0.4 0.4],'fontsize',fSet.Fontsize_text)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CPP2plot = 'CPP';
CPP2plot = 'CPPr';

axH = subtightplot(nrow, ncol, [13 14], subplotgap, subplotmargin, subplotmargin);
% hsub3 = subplot(nrow, ncol, [5]);
plot_subplot_label(axH, fSet, 'C')


switch CPP2plot
    case 'CPP'
        switch dataset
            case 'bigDots'
                t2test      = [300 550];
                XLIM = [-50 700];
            case 'CD'
                t2test      = [750 1200];
                XLIM = [-50 1500];
        end
        ttSPG = allSPG_times-abs(t(1)/1000);
        ttSPG = ttSPG * 1000;
        
        ITPC2plot = CPP_ITPC;
        ITPC_bar2plot = CPP_ITPC_bar;
        ITPC_band2plot = CPP_ITPC_band;
        
        ITPC2plot = CPP_csd_ITPC;
        ITPC_bar2plot = CPP_csd_ITPC_bar;
        ITPC_band2plot = CPP_csd_ITPC_band;
        
        x_ITPC = [45:35:185];
        modelString = 'CPP_ITPC';
        modelString = 'CPP_csd_ITPC';
    case 'CPPr'
        
        t2test = [-300 -50];
        ttSPG = allSPG_timesr-abs(tr(1)/1000);
        ttSPG = ttSPG * 1000;
        
        ITPC2plot = CPPr_ITPC;
        ITPC_bar2plot = CPPr_ITPC_bar;
        ITPC_band2plot = CPPr_ITPC_band;
        
        ITPC2plot = CPPr_csd_ITPC;
        ITPC_bar2plot = CPPr_csd_ITPC_bar;
        ITPC_band2plot = CPPr_csd_ITPC_band;
        
        XLIM = [-500 100];
        x_ITPC = [-450:35:-310];
        modelString = 'CPPr_ITPC';
        modelString = 'CPPr_csd_ITPC';
end
f2test      = [0.1 4]; % in reality 0 4, but this translates to this in the plot

maplimits = [min(min(min(squeeze(mean(ITPC2plot{idx2plot},1))))) max(max(max(squeeze(mean(ITPC2plot{idx2plot},1)))))];
maplimits = [0 max(max(max(squeeze(mean(ITPC2plot{idx2plot},1)))))];

contourf(ttSPG, allSPG_freq, squeeze(mean(mean(ITPC2plot{idx2plot},1),2))',100,'LineColor','none');
set(gca,'clim',maplimits)
ylim([0 18])
axis xy
xlim(XLIM)
% ylim([-0.1 34])
hold on
figInit('ax');
set(axH,'color','white')
colormap(axH,'Default')


h = rectangle('Position',[t2test(1) f2test(1) diff(t2test) diff(f2test)]);%, '-w','linewidth',3)
h.LineWidth = fSet.LineWidth;
h.EdgeColor = [1 1 1];
% set(gca,'TickDir','out')

ylabel('Frequency (Hz)','fontsize',fSet.Fontsize_text)
xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text)

B=colorbar;
xpos = axH.Position(3) + axH.Position(1) - .035 ;
ypos = axH.Position(4) + axH.Position(2) - .08 ;

set(B, 'Position', [xpos ypos .01981 * 0.7 .1 * 0.7], 'Limits', maplimits, 'FontSize', fSet.Fontsize_text_in)
B.AxisLocation = 'in';
B.FontWeight = 'bold';

%%%% ITPC band

axH = subtightplot(nrow, ncol, [15 16 17 18], subplotgap, subplotmargin, subplotmargin);
% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');
plot_subplot_label(axH, fSet, 'D')

clear LEGEND
hold on

tmpCPP = nanmean(ITPC_bar2plot{idx2plot});
for ibin = 1:nbin2use
    if setAlpha
        plot([-1500 700],[tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
    else
        plot([-1500 700],[tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:)]);
    end
end
plot([t2test(1) t2test(2)], [0.12 0.12], 'k', 'linewidth', 4)

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, modelString, nSub, nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');

switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end


% plot model fit
if ~isempty(mfit)
    h = plot(x_ITPC, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([x_ITPC flip(x_ITPC)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

clear h h2 LEGEND
for ibin = 1:nbin2use
    if setAlpha
        h(ibin) = boundedline(ttSPG, squeeze(mean(ITPC_band2plot{idx2plot}(:,ibin, :, :),1)), squeeze(std(ITPC_band2plot{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    else
        h(ibin) = boundedline(ttSPG, squeeze(mean(ITPC_band2plot{idx2plot}(:,ibin, :, :),1)), squeeze(std(ITPC_band2plot{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:));
    end
    set(h(ibin),'linewidth',fSet.LineWidth)
    
    % se
    tmpstd = squeeze(std(ITPC_band2plot{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpCPP(ibin)-tmpstd tmpCPP(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    % plotMarker
    h2(ibin) = plot(x_ITPC(ibin), tmpCPP(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end

yBox = [min(tmpCPP)-min(tmpCPP)*0.07 max(tmpCPP)+min(tmpCPP)*0.07];
xBox = [min(x_ITPC)-25 max(x_ITPC)+25];
text(min(xBox)*0.95, max(yBox)*1.09, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([0 0], [0 0.5] ,'k','linewidth',1)
xlim(XLIM)
        
ylim([0.1 0.5])
% axis square
xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text)
ylabel({'ITPC'},'fontsize',fSet.Fontsize_text)


set(gcf,'renderer','painters')

saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_CPP'];
figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% supp fig . CPP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setAlpha = 1;
subplotgap = [0.1 0.08];
subplotmargin = [0.1 0.1];
[figHandle, fSet] = figInit('fig', 11, {'width',22 ;'height', 7});

nsubplot = 3; % onset, response, ITPC
nrow = 1;
ncol = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CPP2plot = 'CPP';
CPP2plot = 'CPPr';

axH = subtightplot(nrow, ncol, [1], subplotgap, subplotmargin, subplotmargin);
% hsub3 = subplot(nrow, ncol, [5]);
colormap(axH,'Default')
plot_subplot_label(axH, fSet, 'A')


switch CPP2plot
    case 'CPP'
        switch dataset
            case 'bigDots'
                t2test      = [300 550];
                XLIM = [-50 700];
            case 'CD'
                t2test      = [750 1200];
                XLIM = [-50 1500];
        end
        
        ttSPG = allSPG_times-abs(t(1)/1000);
        ttSPG = ttSPG * 1000;
        ffSPG = allSPG_freq;
        
        ITPC2plot = CPP_ITPC;
        ITPC_bar2plot = CPP_ITPC_bar;
        ITPC_band2plot = CPP_ITPC_band;
        modelString = 'CPP_ITPC';
        
%         ITPC2plot = CPP_csd_ITPC;
%         ITPC_bar2plot = CPP_csd_ITPC_bar;
%         ITPC_band2plot = CPP_csd_ITPC_band;
%         modelString = 'CPP_csd_ITPC';
%         
        XLABEL = 'Time from motion onset (ms)';
        x_ITPC = [45:35:185];
    case 'CPPr'
        
        t2test = [-300 -50];

        ttSPG = allSPG_timesr-abs(tr(1)/1000);
        ttSPG = ttSPG * 1000;
        ffSPG = allSPG_freq;
        
        ITPC2plot = CPPr_ITPC;
        ITPC_bar2plot = CPPr_ITPC_bar;
        ITPC_band2plot = CPPr_ITPC_band;
        
        modelString = 'CPPr_ITPC';
        
%         ITPC2plot = CPPr_csd_ITPC;
%         ITPC_bar2plot = CPPr_csd_ITPC_bar;
%         ITPC_band2plot = CPPr_csd_ITPC_band;
%         
%         modelString = 'CPPr_csd_ITPC';
        
        XLIM = [-500 100];
        x_ITPC = [-450:45:-270];
        XLABEL = 'Time from response (ms)';
end
f2test      = [0.1 4]; % in reality 0 4, but this translates to this in the plot

maplimits = [min(min(min(squeeze(mean(ITPC2plot{idx2plot},1))))) max(max(max(squeeze(mean(ITPC2plot{idx2plot},1)))))];
maplimits = [0 max(max(max(squeeze(mean(ITPC2plot{idx2plot},1)))))];

contourf(ttSPG, ffSPG, squeeze(mean(mean(ITPC2plot{idx2plot},1),2))',20,'LineColor','none');
set(gca,'clim',maplimits)
ylim([0 18])
axis xy
xlim(XLIM)
% ylim([-0.1 34])
hold on
figInit('ax');
set(axH,'color','white')


h = rectangle('Position',[t2test(1) f2test(1) diff(t2test) diff(f2test)]);%, '-w','linewidth',3)
h.LineWidth = fSet.LineWidth;
h.EdgeColor = [1 1 1];
% set(gca,'TickDir','out')

ylabel('Frequency (Hz)','fontsize',fSet.Fontsize_text)
xlabel(XLABEL,'fontsize',fSet.Fontsize_text)

B=colorbar;
% set(B, 'Position', [.345 .055 .01981 .1], 'Limits', maplimits)
% xpos = hsub3.Position(3) + hsub3.Position(1) - .01981*figsize - 0.01;
% ypos = hsub3.Position(4) + hsub3.Position(2) - .1*figsize - 0.01;
xpos = axH.Position(3) + axH.Position(1) - .01981 - 0.015;
ypos = axH.Position(4) + axH.Position(2) - .2 ;

set(B, 'Position', [xpos ypos .01981 .15], 'Limits', maplimits, 'FontSize', fSet.Fontsize_text_in)
B.AxisLocation = 'in';
B.FontWeight = 'bold';
% title(B,'ITPC','fontsize',fSet.axFontsize)
% axis square

%%%% ITPC band


for iplot = 1:2
    
    switch iplot
        case 1
            pupilIdx = idx_resp;
            panelLabel = 'B';
        case 2
            pupilIdx = idx_BL_lp;
            panelLabel = 'C';
    end
    
    
    axH = subtightplot(nrow, ncol, [iplot + 1], subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
    plot_subplot_label(axH, fSet, panelLabel)
    
    clear LEGEND
    hold on
    
    tmpCPP = nanmean(ITPC_bar2plot{pupilIdx});
    for ibin = 1:nbin2use
        if setAlpha
            plot([-1500 700],[tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
        else
            plot([-1500 700],[tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:)]);
        end
    end
    plot([t2test(1) t2test(2)], [0.12 0.12], 'k', 'linewidth', 4)
    
    
    
    [~,idx] = grep(allBin2use(pupilIdx),'GLM_pupil');
    if idx == 0
        filename_R_tmp = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{pupilIdx} '_' bintype fileExt_preprocess  fileExt_CDT   ];
    else
        filename_R_tmp = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{pupilIdx} '_' bintype fileExt_preprocess  fileExt_CDT fileExt_GLM  ];
    end
    
    % get stats
    [stats, mfit] = loadStatsR(filename_R_tmp, modelString, nSub, nbin2use);
    [pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
    
    switch stats.model
        case 'Linear'
            betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
            model_color = fSet.colors(2,:);
        case 'Quadratic'
            if stats.U
                betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
            else
                betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
            end
            model_color = fSet.colors(1,:);
        otherwise
            betaString = [];
    end
    
    
    % plot model fit
    if ~isempty(mfit)
        h = plot(x_ITPC, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
        h = patch([x_ITPC flip(x_ITPC)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
        h.FaceAlpha = 0.3;
        h.EdgeAlpha = 0.3;
        h.EdgeColor = [h.FaceColor];
    end
    
    clear h h2 LEGEND
    for ibin = 1:nbin2use
        if setAlpha
            h(ibin) = boundedline(ttSPG, squeeze(mean(ITPC_band2plot{pupilIdx}(:,ibin, :, :),1)), squeeze(std(ITPC_band2plot{pupilIdx}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
        else
            h(ibin) = boundedline(ttSPG, squeeze(mean(ITPC_band2plot{pupilIdx}(:,ibin, :, :),1)), squeeze(std(ITPC_band2plot{pupilIdx}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:));
        end
        set(h(ibin),'linewidth',fSet.LineWidth)
        
        % se
        tmpstd = squeeze(std(ITPC_band2plot{pupilIdx}(:,ibin, :, :)))/sqrt(nSub);
        plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpCPP(ibin)-tmpstd tmpCPP(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
        
        % plotMarker
        h2(ibin) = plot(x_ITPC(ibin), tmpCPP(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
        LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    end
    
    yBox = [min(tmpCPP)-min(tmpCPP)*0.07 max(tmpCPP)+min(tmpCPP)*0.07];
    xBox = [min(x_ITPC)-25 max(x_ITPC)+25];
    text(min(xBox)*0.95, max(yBox)*1.12, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)
    
    h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
    h.LineWidth = 1;
    h.LineStyle = '--';
    
    plot([0 0], [0 0.6] ,'k','linewidth',1)
    xlim(XLIM)
    
    ylim([0.1 0.55])
    % axis square
    xlabel(XLABEL,'fontsize',fSet.Fontsize_text)
    ylabel({'ITPC'},'fontsize',fSet.Fontsize_text)
    
end
saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_CPPsupp'];
figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig 3. Alpha, beta, N2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[fH, fSet] = figInit('fig', 3, {'height', 22; 'width',22});

nrow = 3;
ncol = 8;
subplotgap = [0.1 0.07];
subplotmargin = [0.1 0.1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, [1 2], subplotgap + [0 0.03], subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'A')


%%% alpha - RT
cfg.ylim_RT     = [510 600];
cfg.ylim_RT_CV  = [0.21 0.29];
cfg.ylabel.RT   = 'RT (ms)';

bin2use = allBin2use{idx_alpha};
filename_R_alpha = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_pretarget_alpha_' bintype fileExt_preprocess  fileExt_CDT   ];
[stats, mfit] = loadStatsR(filename_R_alpha, 'RT',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end


ax1 = gca;
hold on
clear hPlotMark

% plot model fit
if ~isempty(mfit)
    h = plot(mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', fSet.colors(1,:));
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
% plot se
for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT{idx_alpha}(:,ibin))-std(RT{idx_alpha}(:,ibin))/sqrt(nSub) ...
        mean(RT{idx_alpha}(:,ibin))+std(RT{idx_alpha}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.LineWidth_in);
end
% plot plotMarker
hPlotMark(ibinType) = plot(x, mean(RT{idx_alpha}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(1,:),'MarkerFaceColor', fSet.colors(1,:), 'MarkerEdgeColor', 'k');

ylim([cfg.ylim_RT])
h = ylabel(cfg.ylabel.RT,'fontsize',fSet.Fontsize_text);
h.Color = fSet.colors(1,:);
set(gca,'xtick',x)
xlabel('\alpha Power Bin','fontsize',fSet.Fontsize_text);
% axis square

h = text(0.5, max(mean(RT{idx_alpha})) * 1.022, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in, 'color', fSet.colors(1,:));

% text for RT_CV
[stats, mfit] = loadStatsR(filename_R_alpha, 'RT_CV',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
text(0.5, mean(RT{idx_alpha}(:,1))*0.92, {betaString,pstring_chi }, 'fontsize',fSet.Fontsize_text_in, 'color', fSet.colors(2,:))

% plot model fit
if ~isempty(mfit)
    h = plot(mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', fSet.colors(2,:));
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

% plot RT_CV
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k','XTickLabel',[]);
figInit('ax');
ax2.YLim = [cfg.ylim_RT_CV];
ax2.XColor = 'none';
linkaxes([ax1 ax2],'x');
hold on

% h = plot(mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', fSet.colors(2,:));
% h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
% h.FaceAlpha = 0.3;
% h.EdgeAlpha = 0.3;
% h.EdgeColor = [h.FaceColor];

for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT_CV{idx_alpha}(:,ibin))-std(RT_CV{idx_alpha}(:,ibin))/sqrt(nSub) ...
        mean(RT_CV{idx_alpha}(:,ibin))+std(RT_CV{idx_alpha}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.LineWidth_in);
end

hPlotMark(ibinType) = plot(x, mean(RT_CV{idx_alpha}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(2,:),'MarkerFaceColor', fSet.colors(2,:), 'MarkerEdgeColor', 'k');
h = ylabel('RT CV','fontsize',fSet.Fontsize_text);
h.Color = fSet.colors(2,:);
% axis square
subplot_ax1 = get(gca);
xlim([0.5 nbin2use+0.5])
xlim([0 nbin2use+1])

%%% aplha-pupil
axH = subtightplot(nrow, ncol, [3 4], subplotgap + [0 0.03], subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'B')

cfg.ylim.alpha = [1.9 2.9];
hold on

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'alpha',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
% plot model fit
if ~isempty(mfit)
    h = plot(1:5, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(alpha{idx2plot}(:,ibin))-std(alpha{idx2plot}(:,ibin))/sqrt(nSub) mean(alpha{idx2plot}(:,ibin))+std(alpha{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'linewidth',fSet.LineWidth_in,'color','k');
    
    plot(x(ibin), mean(alpha{idx2plot}(:,ibin)),plotMarker{ibin}, ...
        'color', fSet.colors(ibin,:),'markersize', fSet.MarkerSize,...
        'MarkerFaceColor', fSet.colors(ibin,:),'MarkerEdgeColor','k')
end

text(1, 2.1, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)


ylim([cfg.ylim.alpha])
h = ylabel('\alpha Power','fontsize',fSet.Fontsize_text);
% xlim([0.5 5.5])
xlim([0 nbin2use+1])

set(gca,'xtick',x)
xlabel('Pupil Bin','fontsize',fSet.Fontsize_text);

% axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, [5:8], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'C')

% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');

cfg.xlim.beta = [-530 80];
cfg.ylim.beta = [0.52 0.76];
cfg.ylim.beta = [-0.2 0.01];


x_beta = [-440:30:-320];

spectral_t = allStft_timesr;

hold on

% beta2use_mean = beta_base_response_amplitude; modelString = {'preRespBeta_base', 'preRespBeta_slope'};
% beta2use = beta_base_response;
% beta2use_mean = beta_response_mean; modelString = {'preRespBeta', 'preRespBeta_slope'};
% beta2use = beta_response;
beta2use_mean = beta_baseAT_response_amplitude; modelString = {'preRespBeta_baseAT', 'preRespBeta_slope'};
beta2use = beta_baseAT_response;


tmpBeta = mean(beta2use_mean{idx2plot});
tmpBeta_slope = mean(beta_pre_response_slope{idx2plot});
for ibin = 1:nbin2use
    plot([-600 100],[tmpBeta(ibin) tmpBeta(ibin)],'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
%
% plot(x_beta, tmpBeta,'linewidth',4,'color',[0.5 0.5 0.5 0.5])
plot([0 0], [cfg.ylim.beta] ,'k','linewidth',1)

twin_bar    = [-130 -70];
t2test_slope= [-300 0];

plot([twin_bar(1) twin_bar(2)], cfg.ylim.beta(2) - [0.017 0.017], 'k', 'linewidth', 4)
plot([t2test_slope(1) t2test_slope(2)], cfg.ylim.beta(2) - [0.009 0.009], 'color',[0.5 0.5 0.5], 'linewidth', 4)

% yBox = [min(tmpBeta)-min(tmpBeta)*0.045 max(tmpBeta)+min(tmpBeta)*0.045];
yBox = [min(tmpBeta)+min(tmpBeta)*0.13 max(tmpBeta)-min(tmpBeta)*0.13];
xBox = [min(x_beta)+min(x_beta)*0.045 max(x_beta)-min(x_beta)*0.045];

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, modelString{1},nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    %     plot model fit
    h = plot(x_beta, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([x_beta flip(x_beta)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

clear h h2 LEGEND
for ibin = 1:nbin2use
    % plot line
    h(ibin) = boundedline(spectral_t, squeeze(mean(beta2use{idx2plot}(:,ibin, :, :),1)), squeeze(std(beta2use{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(h(ibin),'linewidth',fSet.LineWidth)
    
    % amplitude se
    tmpstd = squeeze(std(beta2use_mean{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_beta(ibin) x_beta(ibin)],[tmpBeta(ibin)-tmpstd tmpBeta(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot(x_beta(ibin), tmpBeta(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
end

xlim([cfg.xlim.beta])
ylim([cfg.ylim.beta])
% set(gca,'fontsize',fSet.plFontsize,'ydir','reverse')
set(gca,'ydir','reverse')
xlabel('Time from response (ms)','fontsize',fSet.Fontsize_text)
ylabel({'Power (dB)'},'fontsize',fSet.Fontsize_text)
% ydir('reverse')

% text(xBox(2)*0.99, yBox(1)*1.02, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)
text(xBox(1)*0.99, yBox(2)*0.9, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% beta slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2 = axes('position',[0.61 0.84 0.14 0.07]) ; % inset
figInit('ax',[],{'fontsize',10});
set(gca,'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4])
hold on

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, modelString{2},nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(beta_pre_response_slope{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpBeta_slope(ibin)-tmpstd tmpBeta_slope(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpBeta_slope(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpBeta_slope)*1.15 max(tmpBeta_slope)*0.75])

YTICK = ax2.YTick;
ax2.YTick = linspace(min(YTICK), max(YTICK), 3);

text(3, min(mean(beta_pre_response_slope{idx2plot})) * 1.35, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

set(gca,'xtick',1:nbin2use,'xticklabel',[],'ydir','reverse')
% xlabel('Pupil bin', 'fontsize',fSet.Fontsize_text)
ylabel('LHB Slope', 'fontsize',fSet.Fontsize_text)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2 combined

%%% general settings

% x and y values where to plot N2 latency and amplitude
x_amplitude = [30:25:130];
y_latency   = [1.5:-0.25:0.5];

n2_xlim = [-20 420];
n2_ylim = [-3 1.9];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2, timecourse and latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axH = subtightplot(nrow, ncol, [9:14 17:22 ], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'D')
hold on

%%% N2i
N2_2plot = squeeze(N2i{idx2plot});
N2_latency2plot = squeeze(N2i_latency{idx2plot});
N2_amplitude2plot = squeeze(N2i_amplitude{idx2plot});

% times N2 values have been computed on
t2test_amp      = 340 + [-50 50];
t2test_lat      = [250 450];

plot([0 0], [-5 5] ,' k','linewidth',1)
plot(t2test_amp, [-2.9 -2.9] ,'color',[0.5 0.5 0.5],'linewidth',4)

% plot lines to indicate N2c lat and amp
tmpN2_lat = mean(N2_latency2plot);
for ibin = 1:nbin2use
    plot([tmpN2_lat(ibin) tmpN2_lat(ibin)], n2_ylim,'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
tmpN2_amp = mean(N2_amplitude2plot);
for ibin = 1:nbin2use
    plot(n2_xlim,[tmpN2_amp(ibin) tmpN2_amp(ibin)],'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end

% get stats for latency
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'N2i_latency',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
% plot model fit, latency
if ~isempty(mfit)
    h = plot(mfit.mean, y_latency, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], [y_latency flip(y_latency)] , h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

% draw box and print stats
xBox = [min(tmpN2_lat)-min(tmpN2_lat)*0.05 max(tmpN2_lat)+min(tmpN2_lat)*0.05];
yBox = [min(y_latency)-min(y_latency)*0.35 max(y_latency)+min(y_latency)*0.35];

text(max(xBox)*0.97, max(yBox)*1.2, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([-200 xBox(1)], [0 0] ,' k','linewidth',1)
plot([xBox(2) 800], [0 0] ,' k','linewidth',1)

% get stats, amplitude
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'N2i_amplitude',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(x_amplitude, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([x_amplitude flip(x_amplitude)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

yBox = [min(tmpN2_amp)+min(tmpN2_amp)*0.18 max(tmpN2_amp)-min(tmpN2_amp)*0.18];
xBox = [min(x_amplitude)-min(x_amplitude)*0.5 max(x_amplitude)+min(x_amplitude)*0.5];

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

text(min(xBox)*1.03, min(yBox)*1.1, {[betaString ' ' pstring_chi]}, 'fontsize',fSet.Fontsize_text_in)

clear h h_n2i h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2
    h_n2i(ibin) = boundedline(cfg.t, squeeze(mean(N2_2plot(:,ibin, :),1)), squeeze(std(N2_2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
    set(h_n2i(ibin),'linewidth',fSet.LineWidth,'LineStyle','--')
    
    % se
    tmpN2_lat_std = std(N2_latency2plot(:,ibin))/sqrt(nSub);
    plot([tmpN2_lat(ibin)-tmpN2_lat_std tmpN2_lat(ibin)+tmpN2_lat_std], [y_latency(ibin) y_latency(ibin)],'linewidth',fSet.LineWidth_in,'Color','k')
    h2(ibin) = plot(tmpN2_lat(ibin), y_latency(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
    tmpN2_amp_std = std(N2_amplitude2plot(:,ibin))/sqrt(nSub);
    plot([x_amplitude(ibin) x_amplitude(ibin)], [tmpN2_amp(ibin)-tmpN2_amp_std tmpN2_amp(ibin)+tmpN2_amp_std], 'linewidth',fSet.LineWidth_in,'Color','k')
    h2(ibin) = plot(x_amplitude(ibin),tmpN2_amp(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
end


%%%%%%%%%%%%%%%%%%%
%%% N2c
%%%%%%%%%%%%%%%%%%%

N2_2plot = squeeze(N2c{idx2plot});
N2_latency2plot = squeeze(N2c_latency{idx2plot});
N2_amplitude2plot = squeeze(N2c_amplitude{idx2plot});

% times N2 values have been computed on
t2test_amp      = 266 + [-50 50];
plot(t2test_amp, [-2.8 -2.8] ,'color','k','linewidth',4)

% plot lines to indicate N2c lat and amp
tmpN2_lat = mean(N2_latency2plot);
for ibin = 1:nbin2use
    plot([tmpN2_lat(ibin) tmpN2_lat(ibin)], n2_ylim,'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
tmpN2_amp = mean(N2_amplitude2plot);
for ibin = 1:nbin2use
    plot(n2_xlim,[tmpN2_amp(ibin) tmpN2_amp(ibin)],'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end

% get stats for latency
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'N2c_latency',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
% plot model fit, latency
if ~isempty(mfit)
    h = plot(mfit.mean, y_latency, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], [y_latency flip(y_latency)] , h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

% draw box and print stats
xBox = [min(tmpN2_lat)-min(tmpN2_lat)*0.05 max(tmpN2_lat)+min(tmpN2_lat)*0.05];
yBox = [min(y_latency)-min(y_latency)*0.35 max(y_latency)+min(y_latency)*0.35];

text(max(xBox)*0.97, max(yBox)*1.2, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([-200 xBox(1)], [0 0] ,' k','linewidth',1)
plot([xBox(2) 800], [0 0] ,' k','linewidth',1)

% get stats, amplitude
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'N2c_amplitude',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(x_amplitude, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([x_amplitude flip(x_amplitude)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

yBox = [min(tmpN2_amp)+min(tmpN2_amp)*0.12 max(tmpN2_amp)-min(tmpN2_amp)*0.12];
xBox = [min(x_amplitude)-min(x_amplitude)*0.5 max(x_amplitude)+min(x_amplitude)*0.5];

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

text(min(xBox)*1.03, min(yBox)*1.03, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

clear h h_n2c h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2
    h_n2c(ibin) = boundedline(cfg.t, squeeze(mean(N2_2plot(:,ibin, :),1)), squeeze(std(N2_2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
    set(h_n2c(ibin),'linewidth',fSet.LineWidth)
    
    % se
    tmpN2_lat_std = std(N2_latency2plot(:,ibin))/sqrt(nSub);
    plot([tmpN2_lat(ibin)-tmpN2_lat_std tmpN2_lat(ibin)+tmpN2_lat_std], [y_latency(ibin) y_latency(ibin)],'linewidth',fSet.LineWidth_in,'Color','k')
    h2(ibin) = plot(tmpN2_lat(ibin), y_latency(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
    tmpN2_amp_std = std(N2_amplitude2plot(:,ibin))/sqrt(nSub);
    plot([x_amplitude(ibin) x_amplitude(ibin)], [tmpN2_amp(ibin)-tmpN2_amp_std tmpN2_amp(ibin)+tmpN2_amp_std], 'linewidth',fSet.LineWidth_in,'Color','k')
    h2(ibin) = plot(x_amplitude(ibin),tmpN2_amp(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
end

h_n2 = [h_n2i(1) h_n2c(1)];
[lgd,icons,plots,txt] = legend(h_n2, {'N2i', 'N2c'}, 'Fontsize', fSet.Fontsize_text, 'Box', 'off', 'Location', 'NorthWest');
icons(3).Color = [0 0 0];
icons(5).Color = [0 0 0];
lgd.Position = lgd.Position + [0.03 0 0 0];
%%% general settings
xlim(n2_xlim)
ylim(n2_ylim)

xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text);
ylabel({'Amplitude (\muV)'},'fontsize',fSet.Fontsize_text);

%%%%%%%%%%%%%%%%%%%
%%% ITPC, N2c
%%%%%%%%%%%%%%%%%%%

ITPC2use = N2c_ITPC_bar; modelLabel = 'N2c_ITPC';
% ITPC2use = N2c_256_ITPC_bar; modelLabel = 'N2c_256_ITPC';

axH = subtightplot(nrow, ncol, [15 16], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'E')
hold on
% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, modelLabel,nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(ITPC2use{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(ITPC2use{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.8 1.2])

text(1.5, min(tmpmean) * 0.88, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.Fontsize_text)
ylabel('N2c ITPC', 'fontsize',fSet.Fontsize_text)




%%%%%%%%%%%%%%%%%%%
%%% ITPC, N2i
%%%%%%%%%%%%%%%%%%%

ITPC2use = N2i_ITPC_bar; modelLabel = 'N2i_ITPC';
% ITPC2use = N2i_256_ITPC_bar; modelLabel = 'N2i_256_ITPC';


axH = subtightplot(nrow, ncol, [23 24], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'F')
hold on
% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, modelLabel,nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(ITPC2use{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(ITPC2use{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.8 1.2])

text(1.5, min(tmpmean) * 0.88, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.Fontsize_text)
ylabel('N2i ITPC', 'fontsize',fSet.Fontsize_text)

set(gcf,'renderer','painters')

saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_otherVar'];
figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary figure N2 ITPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ITPC2use = '256';
ITPC2use = '512';


[fH, fSet] = figInit('fig', 3, {'height', 12; 'width',22});

nrow = 2;
ncol = 8;
subplotgap = [0.12 0.1];
subplotmargin = [0.1 0.1];

%%%% ITPC
%%% imagesc plot
plotOrder = [1 2; 9 10];
for iplot = 1:2 % N2c and N2i
    
    switch iplot
        case 1
            %%% N2c
            switch ITPC2use
                case '512'
                    N2_2use = N2c_ITPC;
                    N2_2use_band = N2c_ITPC_band;
                    N2_2use_bar = N2c_ITPC_bar;
                    t2test = [200 400];
                case '256'
                    N2_2use = N2c_256_ITPC;
                    N2_2use_band = N2c_256_ITPC_band;
                    N2_2use_bar = N2c_256_ITPC_bar;
                    t2test = [200 400];
            end
            panelLabel= 'A';
       case 2
           %%% N2i
           switch ITPC2use
               case '512'
                   N2_2use = N2i_ITPC;
                   N2_2use_band = N2i_ITPC_band;
                   N2_2use_bar = N2i_ITPC_bar;
                   t2test      = [250 450];
               case '256'
                   N2_2use = N2i_256_ITPC;
                   N2_2use_band = N2i_256_ITPC_band;
                   N2_2use_bar = N2i_256_ITPC_bar;
                   t2test      = [250 450];
           end
            panelLabel= 'D';
    end
    
    switch ITPC2use
        case '512'
            ttSPG = allSPG_times-abs(t(1)/1000);
            ttSPG = ttSPG * 1000;
            freq2use = allSPG_freq;
            f2test      = [0.1 4]; %
        case '256'
            ttSPG = allSPG_times256-abs(t(1)/1000);
            ttSPG = ttSPG * 1000;
            freq2use = allSPG_freq256;
            f2test      = [0.1 4]; %
    end
    
    axH = subtightplot(nrow, ncol, plotOrder(iplot, :), subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
    plot_subplot_label(axH, fSet, panelLabel)

    hold on
    
    colormap(axH,'Default')
    
    
    maplimits = [min(min(min(squeeze(mean(N2_2use{idx2plot},1))))) max(max(max(squeeze(mean(N2_2use{idx2plot},1)))))];
    maplimits = [0 max(max(max(squeeze(mean(N2_2use{idx2plot},1)))))];
    
    contourf(ttSPG, freq2use, squeeze(mean(mean(N2_2use{idx2plot},1),2))',20,'LineColor','none');
    set(gca,'clim',maplimits)
    ylim([0 18])
    axis xy
    xlim([-50 600])
    hold on
    figInit('ax');

    h = rectangle('Position',[t2test(1) f2test(1) diff(t2test) diff(f2test)]);%, '-w','linewidth',3)
    h.LineWidth = fSet.LineWidth;
    h.EdgeColor = [1 1 1];
    
    ylabel('Frequency (Hz)','fontsize',fSet.Fontsize_text)
%     if iplot==2
    xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text)
%     end
    figInit('ax');
    set(axH,'color','white','linewidth',2)
    B=colorbar;
    % set(B, 'Position', [.345 .375 .01981 .1], 'Limits', maplimits)
    xpos = axH.Position(3) + axH.Position(1) - .01981 - 0.01;
    ypos = axH.Position(4) + axH.Position(2) - .1 - 0.02;
    set(B, 'Position', [xpos ypos .01981 .1], 'Limits', maplimits)
    B.AxisLocation = 'in';
    B.FontWeight = 'bold';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOrder = [...
    3 4 5;
    11 12 13; 
    6 7 8;
    14 15 16;    
    ];
    
for iplot = 1:4 %N2c/Nci, baseline and pupil response
    
    switch iplot
        case 1
            %%% N2c, pupil resp
            pupilIdx = idx_resp;
            switch ITPC2use
                case '512'
                    N2_2use = N2c_ITPC;
                    N2_2use_band = N2c_ITPC_band;
                    N2_2use_bar = N2c_ITPC_bar;
                    stats_label = 'N2c_ITPC';
                    t2test      = [200 400];
                case '256'
                    N2_2use = N2c_256_ITPC;
                    N2_2use_band = N2c_256_ITPC_band;
                    N2_2use_bar = N2c_256_ITPC_bar;
                    stats_label = 'N2c_ITPC';
                    t2test      = [150 200];
            end
            YLABEL = 'N2c ITPC';
            panelLabel = 'B';
        case 2
            %%% N2i, pupil resp
            N2_2use = N2i_ITPC;
            N2_2use_band = N2i_ITPC_band;
            N2_2use_bar = N2i_ITPC_bar;
            stats_label = 'N2i_ITPC';
            pupilIdx = idx_resp;
            t2test      = [250 450];
            YLABEL = 'N2i ITPC';
            panelLabel = 'E';
        case 3
            %%% N2c, pupil baseline
            N2_2use = N2c_ITPC;
            N2_2use_band = N2c_ITPC_band;
            N2_2use_bar = N2c_ITPC_bar;
            
            stats_label = 'N2c_ITPC';
            pupilIdx = idx_BL_lp;
            t2test      = [200 400];
            YLABEL = 'N2c ITPC';
            panelLabel = 'C';
        case 4
            %%% N2i, pupil baseline
            N2_2use = N2i_ITPC;
            N2_2use_band = N2i_ITPC_band;
            N2_2use_bar = N2i_ITPC_bar;
            stats_label = 'N2i_ITPC';
            pupilIdx = idx_BL_lp;
            t2test      = [250 450];
            YLABEL = 'N2i ITPC';
            panelLabel = 'F';
    end
    
    
    [~,idx] = grep(allBin2use(pupilIdx),'GLM_pupil');
    if idx == 0
        filename_R_tmp = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{pupilIdx} '_' bintype fileExt_preprocess  fileExt_CDT   ];
    else
        filename_R_tmp = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{pupilIdx} '_' bintype fileExt_preprocess  fileExt_CDT fileExt_GLM  ];
    end
    
    
    switch ITPC2use
        case '512'
            ttSPG = allSPG_times-abs(t(1)/1000);
            ttSPG = ttSPG * 1000;
            freq2use = allSPG_freq;
            f2test      = [0.1 4]; %
        case '256'
            ttSPG = allSPG_times256-abs(t(1)/1000);
            ttSPG = ttSPG * 1000;
            freq2use = allSPG_freq256;
            f2test      = [0.1 4]; %
    end

    axH = subtightplot(nrow, ncol, plotOrder(iplot, :), subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
    plot_subplot_label(axH, fSet, panelLabel)
    
    x_ITPC = [-80:35:65]+15;
    
    clear LEGEND
    hold on
    
    tmpN2 = nanmean(N2_2use_bar{pupilIdx});
    for ibin = 1:nbin2use
        plot([-200 700],[tmpN2(ibin) tmpN2(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
    end
    plot([t2test(1) t2test(2)], [0.11 0.11], 'k', 'linewidth', 4)
    
    % get stats
    bin2use = allBin2use{pupilIdx};
    [stats, mfit] = loadStatsR(filename_R_tmp, stats_label, nSub, nbin2use);
    [pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
    switch stats.model
        case 'Linear'
            betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
            model_color = fSet.colors(2,:);
        case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
            model_color = fSet.colors(1,:);
        otherwise
            betaString = [];
    end
    if ~isempty(mfit)
        %     plot model fit
        h = plot(x_ITPC, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
        h = patch([x_ITPC flip(x_ITPC)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
        h.FaceAlpha = 0.3;
        h.EdgeAlpha = 0.3;
        h.EdgeColor = [h.FaceColor];
    end
    yBox = [min(tmpN2)-min(tmpN2)*0.07 max(tmpN2)+min(tmpN2)*0.07];
    xBox = [min(x_ITPC)-35 max(x_ITPC)+35];
    text(30, max(yBox)*1.07, {pstring_chi}, 'fontsize',fSet.Fontsize_text)
    
    h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
    h.LineWidth = 1;
    h.LineStyle = '--';
    
    
    clear h h2 LEGEND
    for ibin = 1:nbin2use
        h(ibin) = boundedline(ttSPG, squeeze(mean(N2_2use_band{pupilIdx}(:,ibin, :, :),1)), squeeze(std(N2_2use_band{pupilIdx}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
        
        set(h(ibin),'linewidth',fSet.LineWidth)
        
        % se
        tmpstd = squeeze(std(N2_2use_band{pupilIdx}(:,ibin, :, :)))/sqrt(nSub);
        plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpN2(ibin)-tmpstd tmpN2(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
        
        % plotMarker
        h2(ibin) = plot(x_ITPC(ibin), tmpN2(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
        LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    end
        
    plot([0 0], [0 yBox(1)] ,' k','linewidth',1)
    plot([0 0], [yBox(2) 0.5] ,' k','linewidth',1)
    
    xlim([-130 550])
    ylim([0.1 0.33])
    % axis square
    xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text)
    ylabel(YLABEL, 'fontsize', fSet.Fontsize_text, 'interpret','none')
    
end

saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_N2_ITPC'];
figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4. reduced - CPP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [fH, fSet] = figInit('fig', 5, {'height', 2.5/5; 'width',6/5});
[fH, fSet] = figInit('fig', 5, {'height', 6; 'width',23});

nrow = 1;
ncol = 5;
subplotgap = [0.1 0.05];
subplotmargin = [0.1 0.1];

subplotLabelDisplacement = [0.06 0.01];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CPP onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, [1], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'A', subplotLabelDisplacement)

hold on

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'CPP_onset',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(nanmean(CPP_onset{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(nanstd(CPP_onset{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.95 1.05])

text(1, min(tmpmean) * 0.97, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

set(gca,'xtick',1:nbin2use,'xticklabel',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.Fontsize_text)
ylabel('CPP onset latency (ms)', 'fontsize',fSet.Fontsize_text)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CPP slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CPPr2plot = CPPr_csd_slope; CPP_label = 'CPP_csd_slope2';
% CPPr2plot = CPPr_slope; CPP_label = 'CPP_slope2';

axH = subtightplot(nrow, ncol, [2], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'B', subplotLabelDisplacement)

hold on

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, CPP_label, nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(CPPr2plot{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(CPPr2plot{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.8 1.2])

text(1, min(tmpmean) * 0.87, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

set(gca,'xtick',1:nbin2use,'xticklabel',1:nbin2use)
set(gca,'ytick',get(gca,'ytick'),'yticklabel',get(gca,'ytick'))
xlabel('Pupil bin', 'fontsize', fSet.Fontsize_text)
ylabel('CPP build-up rate', 'fontsize',fSet.Fontsize_text)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CPP amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CPPr2plot = CPPr_csd_amplitude; CPP_label = 'CPPr_csd_amplitude'; YLABEL = 'CPP amplitude (\muV/m^{2})'; %YLIM = 
% CPPr2plot = CPPr_amplitude; CPP_label = 'CPPr_amplitude';  YLABEL = 'CPP amplitude (\muV)';

axH = subtightplot(nrow, ncol, [3], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'C', subplotLabelDisplacement)

hold on

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, CPP_label,nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(CPPr2plot{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(CPPr2plot{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)].* [0.75 1.2] )%

text(1, min(tmpmean) * 0.83, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

set(gca,'xtick',1:nbin2use,'xticklabel',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.Fontsize_text)
ylabel(YLABEL, 'fontsize',fSet.Fontsize_text)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ITPC band



CPP2plot = 'CPP';
% CPP2plot = 'CPPr';

switch CPP2plot
    case 'CPP'
        switch dataset
            case 'bigDots'
                t2test      = [300 550];
                XLIM = [-50 700];
            case 'CD'
                t2test      = [750 1200];
                XLIM = [-50 1500];
        end
        ttSPG = allSPG_times-abs(t(1)/1000);
        ttSPG = ttSPG * 1000;
        
        ITPC2plot = CPP_ITPC;
        ITPC_bar2plot = CPP_ITPC_bar;
        ITPC_band2plot = CPP_ITPC_band;
        
        x_ITPC = [45:45:225];
        modelString = 'CPP_ITPC';
    case 'CPPr'
        
        t2test = [-300 -50];
        ttSPG = allSPG_timesr-abs(tr(1)/1000);
        ttSPG = ttSPG * 1000;
        
        ITPC2plot = CPPr_ITPC;
        ITPC_bar2plot = CPPr_ITPC_bar;
        ITPC_band2plot = CPPr_ITPC_band;
        
        XLIM = [-500 100];
        x_ITPC = [-450:35:-310];
        modelString = 'CPPr_ITPC';
end



axH = subtightplot(nrow, ncol, [4 5], subplotgap, subplotmargin, subplotmargin);
% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');
plot_subplot_label(axH, fSet, 'D', subplotLabelDisplacement)

clear LEGEND
hold on

tmpCPP = nanmean(ITPC_bar2plot{idx2plot});
for ibin = 1:nbin2use
    plot(XLIM,[tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
end
plot([t2test(1) t2test(2)], [0.12 0.12], 'k', 'linewidth', 4)

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, modelString,nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(x_ITPC, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([x_ITPC flip(x_ITPC)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
clear h h2 LEGEND
for ibin = 1:nbin2use
    h(ibin) = boundedline(ttSPG, squeeze(mean(ITPC_band2plot{idx2plot}(:,ibin, :, :),1)), squeeze(std(ITPC_band2plot{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(h(ibin),'linewidth',fSet.LineWidth)
    
    % se
    tmpstd = squeeze(std(ITPC_band2plot{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpCPP(ibin)-tmpstd tmpCPP(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    % plotMarker
    h2(ibin) = plot(x_ITPC(ibin), tmpCPP(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end

yBox = [min(tmpCPP)-min(tmpCPP)*0.07 max(tmpCPP)+min(tmpCPP)*0.07];
xBox = [min(x_ITPC)-min(x_ITPC)*0.6 max(x_ITPC)+min(x_ITPC)*0.6];
% xBox = [min(x_ITPC)-abs(min(x_ITPC))*0.05 max(x_ITPC)+abs(min(x_ITPC))*0.05];
text(min(xBox)*1.1, max(yBox)*1.12, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([0 0], [0 1] ,'k','linewidth',1)
xlim(XLIM)
ylim([0.1 0.5])
% axis square
xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text)
ylabel({'CPP ITPC'},'fontsize',fSet.Fontsize_text)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_redCPP'];

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5. reduced - other var
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fH, fSet] = figInit('fig', 7, {'height', 15; 'width', 22});

nrow = 2;
ncol = 6;
subplotgap = [0.12 0.1];
subplotmargin = [0.1 0.1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, [1 2], subplotgap + [0 0.02], subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'A')

hold on

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'alpha',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(alpha{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(alpha{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.85 1.15])

text(1.5, min(tmpmean) * 0.92, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.Fontsize_text)
ylabel('\alpha power', 'fontsize',fSet.Fontsize_text)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, [3 4 5 6], subplotgap - [0 0.02], subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'B')

cfg.xlim = [-530 80];
cfg.ylim = [0.52 0.76];
cfg.ylim = [-0.2 0.01];


x_beta = [-440:30:-320];

spectral_t = allStft_timesr;

hold on

% beta2use_mean = beta_base_response_amplitude;  modelLabel='preRespBeta_base';
% beta2use = beta_base_response;
beta2use_mean = beta_baseAT_response_amplitude; modelLabel='preRespBeta_baseAT';
beta2use = beta_baseAT_response;
% beta2use_mean = beta_response_mean;
% beta2use = beta_response;

tmpBeta = mean(beta2use_mean{idx2plot});
tmpBeta_slope = mean(beta_pre_response_slope{idx2plot});
for ibin = 1:nbin2use
    plot([-600 100],[tmpBeta(ibin) tmpBeta(ibin)],'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
%
% plot(x_beta, tmpBeta,'linewidth',4,'color',[0.5 0.5 0.5 0.5])
plot([0 0], [cfg.ylim] ,'k','linewidth',1)

twin_bar    = [-130 -70];
% twin_bar    = [-50 50];
t2test_slope= [-300 0];

plot([twin_bar(1) twin_bar(2)], cfg.ylim(2) - [0.017 0.017], 'k', 'linewidth', 4)
plot([t2test_slope(1) t2test_slope(2)], cfg.ylim(2) - [0.009 0.009], 'color',[0.5 0.5 0.5], 'linewidth', 4)

yBox = [min(tmpBeta)+min(tmpBeta)*0.13 max(tmpBeta)-min(tmpBeta)*0.13];
xBox = [min(x_beta)+min(x_beta)*0.045 max(x_beta)-min(x_beta)*0.045];

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, modelLabel,nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    %     plot model fit
    h = plot(x_beta, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([x_beta flip(x_beta)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

clear h h2 LEGEND
for ibin = 1:nbin2use
    % plot line
    h(ibin) = boundedline(spectral_t, squeeze(mean(beta2use{idx2plot}(:,ibin, :, :),1)), squeeze(std(beta2use{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(h(ibin),'linewidth',fSet.LineWidth)
    
    % amplitude se
    tmpstd = squeeze(std(beta2use_mean{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_beta(ibin) x_beta(ibin)],[tmpBeta(ibin)-tmpstd tmpBeta(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot(x_beta(ibin), tmpBeta(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
end


xlim([cfg.xlim])
ylim([cfg.ylim])
% set(gca,'fontsize',fSet.plFontsize,'ydir','reverse')
set(gca,'ydir','reverse')
xlabel('Time from response (ms)','fontsize',fSet.Fontsize_text)
ylabel({'Power (dB)'},'fontsize',fSet.Fontsize_text)
% ydir('reverse')

text(xBox(1)*0.99, yBox(2)*0.9, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% beta slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2 = axes('position',[0.48 0.80 0.14 0.09]) ; % inset
figInit('ax',[],{'fontsize',10});
set(gca,'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4])
hold on

hold on

% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'preRespBeta_slope',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(beta_pre_response_slope{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpBeta_slope(ibin)-tmpstd tmpBeta_slope(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpBeta_slope(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpBeta_slope)*1.15 max(tmpBeta_slope)*0.75])

YTICK = ax2.YTick;
ax2.YTick = linspace(min(YTICK), max(YTICK), 3);
% text(5, max(tmpBeta_slope)*0.85, '*','fontsize',35)

text(3, max(mean(beta_pre_response_slope{idx2plot})) * 1.8, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

set(gca,'xtick',1:nbin2use,'xticklabel',[],'ydir','reverse')
% xlabel('Pupil bin', 'fontsize',fSet.Fontsize_text)
ylabel('LHB slope', 'fontsize',fSet.Fontsize_text)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2, timecourse and latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, [8 9 10], subplotgap, subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'C')
hold on

N2_2plot = squeeze(N2c{idx2plot});
N2_latency2plot = squeeze(N2c_latency{idx2plot});
N2_amplitude2plot = squeeze(N2c_amplitude{idx2plot});

% x and y values where to plot N2 latency and amplitude
y_latency = [0.3:-0.3:-0.9];
x_amplitude = [40:30:160];

% times N2 values have been computed on
t2test_amp      = 266 + [-50 50];
t2test_lat      = [150 400];

plot([0 0], [-4 1] ,' k','linewidth',1)
plot(t2test_amp, [-3 -3] ,' k','linewidth',4)

% plot lines to indicate N2c lat and amp
tmpN2_lat = mean(N2_latency2plot);
for ibin = 1:nbin2use
    plot([tmpN2_lat(ibin) tmpN2_lat(ibin)], [-20 1],'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
tmpN2_amp = mean(N2_amplitude2plot);
for ibin = 1:nbin2use
    plot([0 400],[tmpN2_amp(ibin) tmpN2_amp(ibin)],'linewidth',fSet.LineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end


% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'N2c_latency',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

% plot model fit
if ~isempty(mfit)
    clear h
    h = plot(mfit.mean, y_latency, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    h = patch([mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], [y_latency flip(y_latency)] , h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

% plot box around plotMarker
xBox = [min(tmpN2_lat)-min(tmpN2_lat)*0.05 max(tmpN2_lat)+min(tmpN2_lat)*0.05];
yBox = [min(y_latency)+min(y_latency)*0.22 max(y_latency)-min(y_latency)*0.22];

clear h
h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

text(max(xBox)*1.02, max(yBox)*0.85, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

% zero lines
plot([-200 xBox(1)], [0 0] ,' k','linewidth',1)
plot([xBox(2) 800], [0 0] ,' k','linewidth',1)


% get stats
bin2use = allBin2use{idx2plot};
[stats, mfit] = loadStatsR(filename_R, 'N2c_amplitude',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

% plot model fit
if ~isempty(mfit)
    clear hM
    hM = plot(x_amplitude, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
    hM = patch([x_amplitude flip(x_amplitude)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], [y_latency flip(y_latency)] , hM.Color);
    hM.FaceAlpha = 0.3;
    hM.EdgeAlpha = 0.3;
    hM.EdgeColor = [hM.FaceColor];
end

% plot box around plotMarker
yBox = [min(tmpN2_amp)+min(tmpN2_amp)*0.12 max(tmpN2_amp)-min(tmpN2_amp)*0.12];
xBox = [min(x_amplitude)-min(x_amplitude)*0.5 max(x_amplitude)+min(x_amplitude)*0.5];

clear h
h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

text(min(xBox)*1.02, max(yBox)*0.7, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in)

clear h1 h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2c
    h1 = boundedline(cfg.t, squeeze(mean(N2_2plot(:,ibin, :),1)), squeeze(std(N2_2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
    set(h1,'linewidth',fSet.LineWidth)
    
    % se
    tmpN2_lat_std = std(N2_latency2plot(:,ibin))/sqrt(nSub);
    plot([tmpN2_lat(ibin)-tmpN2_lat_std tmpN2_lat(ibin)+tmpN2_lat_std], [y_latency(ibin) y_latency(ibin)],'linewidth',fSet.LineWidth_in,'Color','k')
    tmpN2_amp_std = std(N2_amplitude2plot(:,ibin))/sqrt(nSub);
    plot([x_amplitude(ibin) x_amplitude(ibin)], [tmpN2_amp(ibin)-tmpN2_amp_std tmpN2_amp(ibin)+tmpN2_amp_std], 'linewidth',fSet.LineWidth_in,'Color','k')
    
    
    plot(tmpN2_lat(ibin), y_latency(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    plot(x_amplitude(ibin),tmpN2_amp(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
end

xlim([-20 350])
ylim([-3.2 0.5])

xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text);
ylabel({'Amplitude (\muV)'},'fontsize',fSet.Fontsize_text);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow, ncol, [11 12], subplotgap + [0 0.02], subplotmargin, subplotmargin);
figInit('ax');
plot_subplot_label(axH, fSet, 'D')

hold on

for itype = 1:2
    switch itype
        case 1
            stats2load = 'N2c_ITPC';
            data2plot = N2c_ITPC_bar;
            LME_text = {'N2c, ', 0.92};
        case 2
            stats2load = 'N2i_ITPC';
            data2plot = N2i_ITPC_bar;
            LME_text = {'N2i, ', 0.92};
            
    end
    [stats, mfit] = loadStatsR(filename_R, stats2load, nSub, nbin2use);
    [pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
    switch stats.model
        case 'Linear'
            betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
            model_color = fSet.colors(2,:);
        case 'Quadratic'
        if stats.U
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
        else
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
        end
            model_color = fSet.colors(1,:);
        otherwise
            betaString = [];
    end
    if ~isempty(mfit)
        % plot model fit
        h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', model_color);
        h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
        h.FaceAlpha = 0.3;
        h.EdgeAlpha = 0.3;
        h.EdgeColor = [h.FaceColor];
    end
    tmpmean = squeeze(mean(data2plot{idx2plot}));
    for ibin = 1:nbin2use
        %     slope
        tmpstd = squeeze(std(data2plot{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
        plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.LineWidth_in,'color','k');
        
        h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
        LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    end
    
    text(1, min(tmpmean) * LME_text{2}, {[LME_text{1}, betaString, pstring_chi]}, 'fontsize',fSet.Fontsize_text_in)
    clear data2plot
end
xlim([0 6])
YLIM = get(gca,'ylim');
ylim([YLIM .* [0.9 1.1]]);

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.Fontsize_text)
ylabel('N2 ITPC', 'fontsize',fSet.Fontsize_text)


saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_redOtherVar'];
figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)



%% fitted pupil diameter

nrow = 1;
ncol = 3;
subplotgap = [0.06 0.08];
subplotmargin = [0.1 0.1];

[figHandle, fSet] = figInit('fig',10, {'height', 8; 'width', 22});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pupil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.t_pupilFit = cfg.t_pupil(cfg.t_pupil>=pupilSet.prestim_range(1));

for iplot = 1:3
    axH = subtightplot(nrow, ncol, iplot, subplotgap, subplotmargin, subplotmargin);
    figInit('ax');

    if iplot == 1
        pupil2plot = pupil_bp;
        tt2use = cfg.t_pupil;
        panelLabel = 'A';
    elseif iplot == 2
        pupil2plot = pupil_conc_GLM_Ramp_yhat;
        tt2use = cfg.t_pupilFit;
        panelLabel = 'B';
    elseif iplot == 3
        pupil2plot = pupil_GLM_Ramp_yhat;
        tt2use = cfg.t_pupilFit;
        panelLabel = 'C';
    end
    cfg.ylim = [-0.2 0.4];
    %     cfg.ylim.Pupil = [-1.2 1.2];
    
    cfg.xlim     = [-250 3000];
    
    plot_subplot_label(axH, fSet, panelLabel)
    hold on
    
    clear LEGEND
    for ibin = 1:nbin2use
        plot(tt2use(1:size(pupil2plot{1},4)),squeeze(nanmean(squeeze(pupil2plot{idx2plot}(:,ibin,:,:)),1)),'linewidth',fSet.LineWidth_in,'Color',fSet.colors(ibin,:));
        LEGEND{ibin} = ['Bin ' num2str(ibin)];
    end
    
    if iplot == 1
        
        [hLeg,icons,plots] = legend(LEGEND,'location','northwest','fontsize',fSet.Fontsize_text , 'box','off');
        % hLeg.Position = hLeg.Position .* [0.95 0.87 1 1];
        hLeg.Position = hLeg.Position .* [0.95 1.1 1 1];
        
        
        XDAT = icons(6).XData;
        icons(6).LineStyle = '-';
        icons(6).LineWidth = fSet.LineWidth;
        icons(6).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
        icons(8).LineStyle = '-';
        icons(8).LineWidth = fSet.LineWidth;
        icons(8).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
        icons(10).LineStyle = '-';
        icons(10).LineWidth = fSet.LineWidth;
        icons(10).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
        icons(12).LineStyle = '-';
        icons(12).LineWidth = fSet.LineWidth;
        icons(12).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
        icons(14).LineStyle = '-';
        icons(14).LineWidth = fSet.LineWidth;
        icons(14).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
        
    end
    xlim(cfg.xlim)
    ylim(cfg.ylim)
    
    for ibin = 1:nbin2use
        A = boundedline(tt2use(1:size(pupil2plot{1},4)),squeeze(nanmean(pupil2plot{idx2plot}(:,ibin,:,:),1)),squeeze(nanstd(pupil2plot{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
        set(A,'linewidth',fSet.LineWidth_in)
    end
    plot([0 0], cfg.ylim ,'k','linewidth',1)
    
    % ylim([-0.2 0.4])
    xlabel('Time from motion onset (ms)','fontsize',fSet.Fontsize_text);
    ylabel('Normalized pupil diameter','fontsize',fSet.Fontsize_text);
    axis square
    
end


saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_pupilFit'];

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)


%% replot all topoplots in higher res
idx2plot    = [find(idx_resp)];

%CPP
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(CPP_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(CPP_topo{idx2plot},1),2))))];
topoplot(squeeze(mean(mean(CPP_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
saveFigName = [bin2use '_' bintype fileExt_preprocess '_CPP_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
% set(gca,'color','none')

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)

% export_fig([paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.png'],'-transparent')

%N2 - left
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(N2c_topo{idx2plot}(:,:,1,:),1),2)))) 1];
topoplot(squeeze(mean(mean(N2c_topo{idx2plot}(:,:,1,:),1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_N2l_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
set(gca,'color','none')

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)


%N2 - right
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(N2c_topo{idx2plot}(:,:,2,:),1),2)))) 1];
topoplot(squeeze(mean(mean(N2c_topo{idx2plot}(:,:,2,:),1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_N2r_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
set(gca,'color','none')

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)

%alpha
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(alpha_preTarget_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(alpha_preTarget_topo{idx2plot},1),2))))];
topoplot(squeeze(mean(mean(alpha_preTarget_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_Alpha_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
set(gca,'color','none')

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)

%beta
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(beta_base_pre_response_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(beta_base_pre_response_topo{idx2plot},1),2))))];
topoplot(squeeze(mean(mean(beta_base_pre_response_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);

saveFigName = [bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM '_Beta_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
set(gca,'color','none')

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)


%% Supplementary figure 2. plot pupil diameter and RT for various pupil measures


subplotgap = [0.08 0.09];
subplotmargin = [0.1 0.1];

[fH, fSet] = figInit('fig', 7, {'height', 18; 'width', 22});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pupil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrow=3;
ncol=4;
subplotOrder = [1 5 9 3 7];

for ibin2use = 1:length(subplotOrder)
    
    switch ibin2use
        case 1
            idx2plot        = strcmpi(allBin2use, 'pupil_bp_RT_neg200_200_regress_bl_iti_side');
            cfg.ylim_RT     = [500 620];
            cfg.ylim_RT_CV  = [0.20 0.3];
        case 2
            idx2plot        = strcmpi(allBin2use, 'pupil_bp_RT_neg200_200_regress_bl_iti_side_RT');
            cfg.ylim_RT     = [500 620];
            cfg.ylim_RT_CV  = [0.20 0.3];
        case 3
            idx2plot        = strcmpi(allBin2use, 'pupil_bp_average_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side');
            cfg.ylim_RT     = [500 620];
            cfg.ylim_RT_CV  = [0.20 0.3];
        case 4
            idx2plot        = strcmpi(allBin2use, 'pupil_bp_slope_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side');
            cfg.ylim_RT     = [500 620];
            cfg.ylim_RT_CV  = [0.20 0.3];
        case 5
            idx2plot        = strcmpi(allBin2use, 'pupil_bp_linearProjection_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side');
            cfg.ylim_RT     = [500 620];
            cfg.ylim_RT_CV  = [0.20 0.3];
    end
    
    [~,idx] = grep(allBin2use(idx2plot),'GLM_pupil');
    if idx == 0
        tmp_filename_R = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{idx2plot} '_' bintype fileExt_preprocess  fileExt_CDT   ];
    else
        tmp_filename_R = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{idx2plot} '_' bintype fileExt_preprocess  fileExt_CDT fileExt_GLM  ];
    end

    
    axH = subtightplot(nrow, ncol, subplotOrder(ibin2use), subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
    plot_subplot_label(axH, fSet, ibin2use)
    
    pupil2plot = pupil_bp;
    cfg.ylim = [-0.2 0.4];
    cfg.xlim = [-400 3000];
    
    hold on
    
    for ibin = 1:nbin2use
        A = boundedline(cfg.t_pupil(1:size(pupil2plot{1},4)),squeeze(mean(pupil2plot{idx2plot}(:,ibin,:,:),1)),squeeze(std(pupil2plot{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
        set(A,'linewidth',fSet.LineWidth_in)
    end
    plot([0 0], cfg.ylim ,'k','linewidth',1)
    
    xlim(cfg.xlim)
    ylim(cfg.ylim)
    xlabel('Time (ms)','fontsize',fSet.Fontsize_text);
    ylabel('Pupil diameter (a.u.)','fontsize',fSet.Fontsize_text);
    axis square
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% RT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    axH = subtightplot(nrow, ncol,subplotOrder(ibin2use)+1, subplotgap, subplotmargin, subplotmargin);
    axH.Position = axH.Position + [-0.02 0 0 0];
    figInit('ax');
    % plot_subplot_label(axH, fSet, 'E')
    
    [stats, mfit] = loadStatsR(tmp_filename_R, 'RT',nSub,nbin2use);
    [pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
    
    switch stats.model
        case 'Linear'
            betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        case 'Quadratic'
            if stats.U
                betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
            else
                betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
            end
        otherwise
            betaString = [];
    end
    

    ax1 = gca;
    hold on
    clear hPlotMark
    
    if ~isempty(mfit)
        % plot model fit
        h = plot(mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', fSet.colors(1,:));
        h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
        h.FaceAlpha = 0.3;
        h.EdgeAlpha = 0.3;
        h.EdgeColor = [h.FaceColor];
    end
    % plot se
    for ibin = 1:nbin2use
        plot([x(ibin) x(ibin)],...
            [mean(RT{idx2plot}(:,ibin))-std(RT{idx2plot}(:,ibin))/sqrt(nSub) ...
            mean(RT{idx2plot}(:,ibin))+std(RT{idx2plot}(:,ibin))/sqrt(nSub) ],...
            'color', 'k','linewidth',fSet.LineWidth_in);
    end
    % plot plotMarker
    hPlotMark(ibinType) = plot(x, mean(RT{idx2plot}),plotMarker{1},...
        'markersize', fSet.MarkerSize,...
        'color', fSet.colors(1,:),'MarkerFaceColor', fSet.colors(1,:), 'MarkerEdgeColor', 'k');
    
    ylim([cfg.ylim_RT])
    h = ylabel('RT (ms)','fontsize',fSet.Fontsize_text);
    h.Color = fSet.colors(1,:);
    set(gca,'xtick',x)
    xlabel('Pupil Bin','fontsize',fSet.Fontsize_text);
    axis square
    
%     h = text(0.7, mean(RT{idx2plot}(:,1))*text_y, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in, 'color', fSet.colors(1,:));
    
    % text for RT_CV
    [stats, mfit] = loadStatsR(tmp_filename_R, 'RT_CV',nSub,nbin2use);
    [pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
    switch stats.model
        case 'Linear'
            betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
            %         model_color = fSet.colors(2,:);
        case 'Quadratic'
            if stats.U
                betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{+})'];
            else
                betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f') ' (U^{-})'];
            end
        otherwise
            betaString = [];
    end
    
%     text(3.3, mean(RT{idx2plot}(:,1))*text_y, {betaString, pstring_chi }, 'fontsize',fSet.Fontsize_text_in, 'color', fSet.colors(2,:))
    
    % plot RT_CV
    ax2 = axes('Position',get(ax1,'Position'),...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k','YColor','k','XTickLabel',[]);
    figInit('ax');
    ax2.YLim = [cfg.ylim_RT_CV];
    ax2.XColor = 'none';
    linkaxes([ax1 ax2],'x');
    hold on
    
    if ~isempty(mfit)
        h = plot(mfit.mean, 'linewidth', fSet.LineWidth_in, 'color', fSet.colors(2,:));
        h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
        h.FaceAlpha = 0.3;
        h.EdgeAlpha = 0.3;
        h.EdgeColor = [h.FaceColor];
    end
    for ibin = 1:nbin2use
        plot([x(ibin) x(ibin)],...
            [mean(RT_CV{idx2plot}(:,ibin))-std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ...
            mean(RT_CV{idx2plot}(:,ibin))+std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ],...
            'color', 'k','linewidth',fSet.LineWidth_in);
    end
    
    hPlotMark(ibinType) = plot(x, mean(RT_CV{idx2plot}),plotMarker{1},...
        'markersize', fSet.MarkerSize,...
        'color', fSet.colors(2,:),'MarkerFaceColor', fSet.colors(2,:), 'MarkerEdgeColor', 'k');
    h = ylabel('RT CV','fontsize',fSet.Fontsize_text);
    h.Color = fSet.colors(2,:);
    axis square
    subplot_ax1 = get(gca);
    
    xlim([0.5 nbin2use+0.5])
end
saveFigName = ['SuppFig_allRT'];

figSave(saveFigName, [paths.pop 'fig' filesep 'manuscript' filesep], figureFileType)




