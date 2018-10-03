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

idx_BL_bp       = strcmpi(allBin2use, 'pupil_bp_baseline_regress_iti_side');
idx_BL_lp       = strcmpi(allBin2use, 'pupil_lp_baseline_regress_iti_side');
idx_resp        = strcmpi(allBin2use, 'pupil_lp_RT_neg200_200_regress_bl_iti_side');
% idx_resp        = strcmpi(allBin2use, 'pupil_bp_RT_neg200_200_regress_bl_iti_side');

idx_alpha       = strcmpi(allBin2use, 'pretarget_alpha');
idx_alpha_asym  = strcmpi(allBin2use, 'pretarget_alpha_asym');

idx_N2i         = strcmpi(allBin2use, 'N2i_amplitude_regress_iti_side');

% idx2plot   = idx_BL_bp;
idx2plot   = idx_resp;
% idx2plot   = idx_BL_lp;

% plot Settings
chanlocs = readlocs('cap64.loc'); %biosemi
[~,plot_chans, exclude_chans] = getChannelName;

figureFileType = {'png','eps'};
clear cfg
cfg.t = t;
cfg.tr = tr;
plotMarker = {'o','s','d','^','h'};


%% Figure 1

nrow = 1;
ncol = 3;
subplotgap = [0.06 0.08];

[figHandle, fSet] = figInit('fig',1, {'height', 1/3; 'width', 6/5});
set(gca,'linewidth',2)

set(gca,'color','none','XColor','k','YColor','k')
set(gcf,'PaperPositionMode','auto')
set(gca,'fontsize',fSet.plFontsize)

switch allBin2use{idx2plot}
    case {'pupil_lp_baseline_regress_iti_side','pupil_bp_baseline_regress_iti_side'}
        text_y = 0.855;
    case {'pupil_lp_RT_neg200_200_regress_bl_iti_side','pupil_bp_RT_neg200_200_regress_bl_iti_side'}
        text_y = 0.845;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% paradigm figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subtightplot(nrow, ncol, 1, subplotgap)

% import paradimg png
fig2import = ['C:\Jochem\Dropbox\Monash\bigDots_st\figures\paradigm_test.png'];
paradigm = imread(fig2import);

image(paradigm)
axis off
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pupil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subtightplot(nrow, ncol, 2, subplotgap)
figInit('ax');

pupil2plot = pupil_lp;
% pupil2plot = pupil_bp; 

cfg.xlim.Pupil     = [-150 1200];
cfg.ylim.Pupil      = [-0.2 0.45];
cfg.ylim.Pupil      = [-0.05 0.09];

hold on

clear LEGEND
for ibin = 1:nbin2use
    plot(cfg.t,squeeze(mean(squeeze(pupil2plot{idx2plot}(:,ibin,:,:)),1)),'linewidth',fSet.plLineWidth_in,'Color',fSet.colors(ibin,:));
    LEGEND{ibin} = ['Bin ' num2str(ibin)];
end
[hLeg,icons,plots] = legend(LEGEND,'location','northwest','fontsize',fSet.axLabFontsize , 'box','off');
hLeg.Position = hLeg.Position .* [0.95 0.87 1 1];


XDAT = icons(6).XData;
icons(6).LineStyle = '-';
icons(6).LineWidth = fSet.plLineWidth;
icons(6).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
icons(8).LineStyle = '-';
icons(8).LineWidth = fSet.plLineWidth;
icons(8).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
icons(10).LineStyle = '-';
icons(10).LineWidth = fSet.plLineWidth;
icons(10).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
icons(12).LineStyle = '-';
icons(12).LineWidth = fSet.plLineWidth;
icons(12).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];
icons(14).LineStyle = '-';
icons(14).LineWidth = fSet.plLineWidth;
icons(14).XData = [XDAT(1)+XDAT(2)/2 XDAT(2)];


xlim(cfg.xlim.Pupil)
ylim(cfg.ylim.Pupil)

for ibin = 1:nbin2use
    A = boundedline(cfg.t,squeeze(mean(pupil2plot{idx2plot}(:,ibin,:,:),1)),squeeze(std(pupil2plot{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(A,'linewidth',fSet.plLineWidth_in)
end
plot([0 0], [-0.2 0.4] ,'k','linewidth',1)

% ylim([-0.2 0.4])
xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
ylabel('Normalized pupil diameter','fontsize',fSet.axLabFontsize);
axis square


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subtightplot(nrow, ncol, 3, subplotgap)
figInit('ax');

RT = RT;

cfg.ylim.RT     = [470 620];
cfg.ylim.RT_CV  = [0.195 0.28];
cfg.ylabel.RT   = 'RT (ms)';

bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'RT',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');

switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        %         model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        %         model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end


ax1 = gca;
hold on
clear hPlotMark

if ~isempty(mfit)
    % plot model fit
    h = plot(mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', fSet.colors(1,:));
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
        'color', 'k','linewidth',fSet.plLineWidth_in);
end
% plot plotMarker
hPlotMark(ibinType) = plot(x, mean(RT{idx2plot}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(1,:),'MarkerFaceColor', fSet.colors(1,:), 'MarkerEdgeColor', 'k');

ylim([cfg.ylim.RT])
h = ylabel(cfg.ylabel.RT,'fontsize',fSet.axLabFontsize);
h.Color = fSet.colors(1,:);
set(gca,'xtick',x)
xlabel('Pupil Bin','fontsize',fSet.axLabFontsize);
axis square

h = text(0.7, mean(RT{idx2plot}(:,1))*text_y, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize, 'color', fSet.colors(1,:));

% text for RT_CV
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'RT_CV',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        %         model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        %         model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

text(3.7, mean(RT{idx2plot}(:,1))*text_y, {betaString, pstring_chi }, 'fontsize',fSet.plFontsize, 'color', fSet.colors(2,:))

% plot RT_CV
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k','XTickLabel',[]);
figInit('ax');
ax2.YLim = [cfg.ylim.RT_CV];
ax2.XColor = 'none';
linkaxes([ax1 ax2],'x');
hold on

if ~isempty(mfit)
    h = plot(mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', fSet.colors(2,:));
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT_CV{idx2plot}(:,ibin))-std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ...
        mean(RT_CV{idx2plot}(:,ibin))+std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.plLineWidth_in);
end

hPlotMark(ibinType) = plot(x, mean(RT_CV{idx2plot}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(2,:),'MarkerFaceColor', fSet.colors(2,:), 'MarkerEdgeColor', 'k');
h = ylabel('RT CV','fontsize',fSet.axLabFontsize);
h.Color = fSet.colors(2,:);
axis square
subplot_ax1 = get(gca);

xlim([0.5 nbin2use+0.5])

saveFigName = [bin2use '_' bintype fileExt '_RT'];

for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            %             print(gcf,['-d' figureFileType{ifiletype} 'c'],[paths.pop 'fig' filesep 'p_level' filesep saveFigName '.' figureFileType{ifiletype}])
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
            
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CPP ALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setAlpha = 1;
subplotgap = [0.08 0.07];

figsize = 4/5;

[figHandle, fSet] = figInit('fig', 2, {'width',6/5 ;'height', figsize});

nsubplot = 3; % onset, response, ITPC
nrow = 3;
ncol = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% stim locked/ onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub1 = subtightplot(nrow, ncol, [1:6], subplotgap);
figInit('ax');

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
        plot([tmpCPP(ibin) tmpCPP(ibin)], [-2 25],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
    else
        plot([tmpCPP(ibin) tmpCPP(ibin)], [-2 25],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) ]);
    end
end

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'CPP_onset',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

% plot model fit
if ~isempty(mfit)
    h = plot(mfit.mean,y_onset, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
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
    set(hLine(ibin),'linewidth',fSet.plLineWidth)
    
    % plot SE lines on plotMarker
    tmpstd = nanstd(CPP_onset2plot(:,ibin))/sqrt(nSub);
    plot([tmpCPP(ibin)-tmpstd tmpCPP(ibin)+tmpstd], [y_onset(ibin) y_onset(ibin)],'linewidth',fSet.plLineWidth_in,'color','k');
    
    % plot plotMarker
    h2(ibin) = plot(tmpCPP(ibin), y_onset(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    
end
[hLeg,icons,plots] = legend([h2], LEGEND,'location','southeast','fontsize',fSet.axLabFontsize, 'box','off');

Pos = hLeg.Position;
% hLeg.Position(2) = Pos(2)+Pos(1)*0.03;
hLeg.Position(1) = Pos(1)*1.05;
hLeg.Position(2) = Pos(2)*1.01;

icons(1).FontSize = fSet.axLabFontsize;
icons(2).FontSize = fSet.axLabFontsize;
icons(3).FontSize = fSet.axLabFontsize;
icons(4).FontSize = fSet.axLabFontsize;
icons(5).FontSize = fSet.axLabFontsize;
% POS = icons(1).Position;
% icons(1).Position = POS + [0.1 0 0];
% POS = icons(2).Position;
% icons(2).Position = POS + [0.1 0 0];
% POS = icons(3).Position;
% icons(3).Position = POS + [0.1 0 0];
% POS = icons(4).Position;
% icons(4).Position = POS + [0.1 0 0];
% POS = icons(5).Position;
% icons(5).Position = POS + [0.1 0 0];

% XDAT = icons(6).XData;
icons(6).LineStyle = '-';
icons(6).LineWidth = fSet.plLineWidth;
% icons(6).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];
icons(8).LineStyle = '-';
icons(8).LineWidth = fSet.plLineWidth;
% icons(8).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];
icons(10).LineStyle = '-';
icons(10).LineWidth = fSet.plLineWidth;
% icons(10).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];
icons(12).LineStyle = '-';
icons(12).LineWidth = fSet.plLineWidth;
% icons(12).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];
icons(14).LineStyle = '-';
icons(14).LineWidth = fSet.plLineWidth;
% icons(14).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];

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

text(max(xBox)*1.02, max(yBox)*0.95, [{betaString, pstring_chi}], 'fontsize',fSet.plFontsize)
% text(max(xBox)*1.02, max(yBox)*0.95, [pstring_Q], 'fontsize',fSet.plFontsize)

xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize)
% ylabel({'CPP \muV'},'fontsize',fSet.axLabFontsize)
ylabel({'Amplitude (\muV)'},'fontsize',fSet.plFontsize)


%%%% inset with topoplot
ax1 = axes('position',[0.1 0.76 0.4 0.15]) ; % inset
maplimits = [min(min(squeeze(mean(mean(CPP_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(CPP_topo{idx2plot},1),2))))];
topoplot(squeeze(mean(mean(CPP_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
colormap(ax1,'Jet')
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% resp locked/ threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub2 = subtightplot(nrow, ncol, [7:12], subplotgap);
figInit('ax');
hold on

CPP2plot        = squeeze(CPPr_csd{idx2plot});
CPP_amp2plot    = squeeze(CPPr_amplitude{idx2plot});
CPP_slope2plot  = squeeze(CPPr_csd_slope{idx2plot});

% location to plot plotMarker CPP threshold
x_amp = [-320:15:-250]+140;

t2test_amp      = [-100 0]; % time threshold is computed on
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

% plot times amplitude and slope are computed over
plot([t2test_amp(1) t2test_amp(2)], [1.3 1.3], 'k', 'linewidth', 4)
plot([t2test_slope(1) t2test_slope(2)], [0.6 0.6], 'color',[0.5 0.5 0.5], 'linewidth', 4)

clear h h2 LEGEND
for ibin = 1:nbin2use
    % plot CPP timecourse
    if setAlpha
        h(ibin) = boundedline(cfg.tr, squeeze(mean(CPP2plot(:,ibin, :),1)), squeeze(std(CPP2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
        %     h(ibin) = plot(cfg.tr,squeeze(mean(CPP2plot(:,ibin, :),1)),'color',fSet.colors(ibin,:));
    else
        h(ibin) = boundedline(cfg.tr, squeeze(mean(CPP2plot(:,ibin, :),1)), squeeze(std(CPP2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:));
    end
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % plot plotMarker and se
    tmpstd = std(CPP_amp2plot(:,ibin))/sqrt(nSub);
    plot([x_amp(ibin) x_amp(ibin)], [tmpCPP(ibin)-tmpstd tmpCPP(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    h2(ibin) = plot(x_amp(ibin), tmpCPP(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end

plot([0 0], [-5 40] ,'-k','linewidth',1)
xlim([-350 50])
ylim([-2 30])

yBox = [min(tmpCPP)-min(tmpCPP)*0.1 max(tmpCPP)+min(tmpCPP)*0.1];
xBox = [min(x_amp)-min(x_amp)*0.38 max(x_amp)+min(x_amp)*0.38];

h = rectangle('Position',[xBox(2) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'CPPr_amplitude',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
%     betaString = ['\beta2 = ' num2str(stats.B, '%1.3f')];
%     betaString = ['B2 = ' num2str(stats.B, '%1.2f')];
text(min(xBox)*0.99, max(yBox)*1.07, {pstring_chi}, 'fontsize',fSet.plFontsize)

xlabel('Time from response (ms)','fontsize',fSet.axLabFontsize)
ylabel({'Amplitude (\muV/m^{2})'},'fontsize',fSet.axLabFontsize)
%     ylabel({'CSD ()'},'fontsize',fSet.axLabFontsize)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% resp locked/ slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes('position',[0.155 0.48 0.22 0.09]) ; % inset
figInit('ax',[],{'fontsize',10});
set(gca,'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4])
hold on

tmpCPP_slope = mean(CPP_slope2plot);
% plot(1:5, tmpCPP_slope,'linewidth',4,'color',[0.5 0.5 0.5])

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'CPP_slope2',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
text(min(xBox)*0.95, max(yBox)*1.11, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

% plot model fit
if ~isempty(mfit)
    h = plot(1:5, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
for ibin = 1:nbin2use
    % se
    tmpCPP_slope_std = squeeze(std(CPP_slope2plot(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpCPP_slope(ibin)-tmpCPP_slope_std tmpCPP_slope(ibin)+tmpCPP_slope_std],'linewidth',fSet.plLineWidth_in,'color','k');
    % plotMarker
    h2(ibin) = plot((ibin), tmpCPP_slope(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize, 'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0.5 5.5])
ylim([min(tmpCPP_slope)*0.75 max(tmpCPP_slope)*1.15])

text(1, max(tmpCPP_slope)*1.3, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)


set(gca,'xtick',[])
xlabel('Pupil bin','color',[0.4 0.4 0.4],'fontsize',fSet.plFontsize)
ylabel('Slope','color',[0.4 0.4 0.4],'fontsize',fSet.plFontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub3 = subtightplot(nrow, ncol, [13 14], subplotgap);
% hsub3 = subplot(nrow, ncol, [5]);
colormap(hsub3,'Default')

t2test      = [300 550];
f2test      = [0.1 4]; % in reality 0 4, but this translates to this in the plot

ttSPG = SPG_times-abs(t(1)/1000);
ttSPG = ttSPG * 1000;
maplimits = [min(min(min(squeeze(mean(CPP_ITPC{idx2plot},1))))) max(max(max(squeeze(mean(CPP_ITPC{idx2plot},1)))))];
maplimits = [0 max(max(max(squeeze(mean(CPP_ITPC{idx2plot},1)))))];

contourf(ttSPG, SPG_freq, squeeze(mean(mean(CPP_ITPC{idx2plot},1),2))',20,'LineColor','none');
set(gca,'clim',maplimits)
ylim([0 18])
axis xy
xlim([-50 700])
% ylim([-0.1 34])
hold on
figInit('ax');
set(hsub3,'color','white')


h = rectangle('Position',[t2test(1) f2test(1) diff(t2test) diff(f2test)]);%, '-w','linewidth',3)
h.LineWidth = fSet.plLineWidth;
h.EdgeColor = [1 1 1];
% set(gca,'TickDir','out')

ylabel('Frequency (Hz)','fontsize',fSet.axLabFontsize)
xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize)

B=colorbar;
% set(B, 'Position', [.345 .055 .01981 .1], 'Limits', maplimits)
% xpos = hsub3.Position(3) + hsub3.Position(1) - .01981*figsize - 0.01;
% ypos = hsub3.Position(4) + hsub3.Position(2) - .1*figsize - 0.01;
xpos = hsub3.Position(3) + hsub3.Position(1) - .01981 - 0.01;
ypos = hsub3.Position(4) + hsub3.Position(2) - .1 ;

set(B, 'Position', [xpos ypos .01981*figsize .1*figsize], 'Limits', maplimits, 'FontSize', fSet.plFontsize)
B.AxisLocation = 'in';
B.FontWeight = 'bold';
% title(B,'ITPC','fontsize',fSet.axFontsize)
% axis square

%%%% ITPC band

hsubplot = subtightplot(nrow, ncol, [15 16 17 18], subplotgap);
% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');


x_ITPC = [45:35:185];

ttSPG = SPG_times-abs(t(1)/1000);
ttSPG = ttSPG*1000;
clear LEGEND
hold on

t2test      = [300 550];

tmpCPP = nanmean(CPP_ITPC_bar{idx2plot});
for ibin = 1:nbin2use
    if setAlpha
        plot([-200 700],[tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
    else
        plot([-200 700],[tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:)]);
    end
end
plot([t2test(1) t2test(2)], [0.12 0.12], 'k', 'linewidth', 4)

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'CPP_ITPC',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
% betaString = ['\beta2 = ' num2str(stats.B, '%1.3f')];
% betaString = ['B2 = ' num2str(stats.B, '%1.2f')];

% plot model fit
h = plot(x_ITPC, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', fSet.colors(1,:));
h = patch([x_ITPC flip(x_ITPC)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
h.FaceAlpha = 0.3;
h.EdgeAlpha = 0.3;
h.EdgeColor = [h.FaceColor];

clear h h2 LEGEND
for ibin = 1:nbin2use
    if setAlpha
        h(ibin) = boundedline(ttSPG, squeeze(mean(CPP_ITPC_band{idx2plot}(:,ibin, :, :),1)), squeeze(std(CPP_ITPC_band{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    else
        h(ibin) = boundedline(ttSPG, squeeze(mean(CPP_ITPC_band{idx2plot}(:,ibin, :, :),1)), squeeze(std(CPP_ITPC_band{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:));
    end
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % se
    tmpstd = squeeze(std(CPP_ITPC_band{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpCPP(ibin)-tmpstd tmpCPP(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    % plotMarker
    h2(ibin) = plot(x_ITPC(ibin), tmpCPP(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end

yBox = [min(tmpCPP)-min(tmpCPP)*0.07 max(tmpCPP)+min(tmpCPP)*0.07];
xBox = [min(x_ITPC)-min(x_ITPC)*0.6 max(x_ITPC)+min(x_ITPC)*0.6];
text(min(xBox)*0.95, max(yBox)*1.09, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([0 0], [0 0.5] ,'k','linewidth',1)
xlim([-50 600])
ylim([0.1 0.5])
% axis square
xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize)
ylabel({'ITPC'},'fontsize',fSet.axLabFontsize)


set(gcf,'renderer','painters')

saveFigName = [bin2use '_' bintype fileExt '_CPP'];
for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig 3. Alpha, beta, N2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[fH, fSet] = figInit('fig', 3, {'height', 3.4/5; 'width',6/5});

nrow = 3;
ncol = 8;
subplotgap = [0.09 0.058];
subplotgap2 = [0.09 0.09];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [1 2], subplotgap2);
figInit('ax');


%%% alpha - RT
cfg.ylim.RT     = [490 600];
cfg.ylim.RT_CV  = [0.20 0.28];
cfg.ylabel.RT   = 'RT (ms)';

bin2use = allBin2use{idx_alpha};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'RT',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end


ax1 = gca;
hold on
clear hPlotMark

% plot model fit
if ~isempty(mfit)
    h = plot(mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', fSet.colors(1,:));
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
        'color', 'k','linewidth',fSet.plLineWidth_in);
end
% plot plotMarker
hPlotMark(ibinType) = plot(x, mean(RT{idx_alpha}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(1,:),'MarkerFaceColor', fSet.colors(1,:), 'MarkerEdgeColor', 'k');

ylim([cfg.ylim.RT])
h = ylabel(cfg.ylabel.RT,'fontsize',fSet.axLabFontsize);
h.Color = fSet.colors(1,:);
set(gca,'xtick',x)
xlabel('\alpha Power Bin','fontsize',fSet.axLabFontsize);
% axis square

h = text(0.5, mean(RT{idx_alpha}(:,1))*0.91, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize, 'color', fSet.colors(1,:));

% text for RT_CV
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'RT_CV',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
text(0.5, mean(RT{idx_alpha}(:,1))*0.885, {betaString,pstring_chi }, 'fontsize',fSet.plFontsize, 'color', fSet.colors(2,:))

% plot RT_CV
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k','XTickLabel',[]);
figInit('ax');
ax2.YLim = [cfg.ylim.RT_CV];
ax2.XColor = 'none';
linkaxes([ax1 ax2],'x');
hold on

% h = plot(mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', fSet.colors(2,:));
% h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
% h.FaceAlpha = 0.3;
% h.EdgeAlpha = 0.3;
% h.EdgeColor = [h.FaceColor];

for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT_CV{idx_alpha}(:,ibin))-std(RT_CV{idx_alpha}(:,ibin))/sqrt(nSub) ...
        mean(RT_CV{idx_alpha}(:,ibin))+std(RT_CV{idx_alpha}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.plLineWidth_in);
end

hPlotMark(ibinType) = plot(x, mean(RT_CV{idx_alpha}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(2,:),'MarkerFaceColor', fSet.colors(2,:), 'MarkerEdgeColor', 'k');
h = ylabel('RT CV','fontsize',fSet.axLabFontsize);
h.Color = fSet.colors(2,:);
% axis square
subplot_ax1 = get(gca);
xlim([0.5 nbin2use+0.5])
xlim([0 nbin2use+1])

%%% aplha-pupil
hsub = subtightplot(nrow, ncol, [3 4], subplotgap2);

% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');

cfg.ylim.alpha = [1.72 2.3];
cfg.ylim.alpha = [2.1 3];
hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'alpha',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
% plot model fit
if ~isempty(mfit)
    h = plot(1:5, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(alpha{idx2plot}(:,ibin))-std(alpha{idx2plot}(:,ibin))/sqrt(nSub) mean(alpha{idx2plot}(:,ibin))+std(alpha{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'linewidth',fSet.plLineWidth_in,'color','k');
    
    plot(x(ibin), mean(alpha{idx2plot}(:,ibin)),plotMarker{ibin}, ...
        'color', fSet.colors(ibin,:),'markersize', fSet.MarkerSize,...
        'MarkerFaceColor', fSet.colors(ibin,:),'MarkerEdgeColor','k')
end

text(1, 2.2, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)


ylim([cfg.ylim.alpha])
h = ylabel('\alpha Power','fontsize',fSet.axLabFontsize);
% xlim([0.5 5.5])
xlim([0 nbin2use+1])

set(gca,'xtick',x)
xlabel('Pupil Bin','fontsize',fSet.axLabFontsize);

% axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub1 = subtightplot(nrow, ncol, [5:8], subplotgap);
figInit('ax');

% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');

cfg.xlim.beta = [-530 80];
cfg.ylim.beta = [0.52 0.76];
cfg.ylim.beta = [-0.2 0.01];


x_beta = [-440:30:-320];

spectral_t = stft_timesr;

hold on

beta2use_mean = beta_base_response_amplitude;
beta2use = beta_base_response;
% beta2use_mean = beta_response_mean;
% beta2use = beta_response;

tmpBeta = mean(beta2use_mean{idx2plot});
tmpBeta_slope = mean(beta_pre_response_slope{idx2plot});
for ibin = 1:nbin2use
    plot([-600 100],[tmpBeta(ibin) tmpBeta(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
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
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
% [stats, mfit] = loadStatsR(filename, 'preRespBeta',nSub,nbin2use);
[stats, mfit] = loadStatsR(filename, 'preRespBeta_base',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    %     plot model fit
    h = plot(x_beta, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([x_beta flip(x_beta)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

clear h h2 LEGEND
for ibin = 1:nbin2use
    % plot line
    h(ibin) = boundedline(spectral_t, squeeze(mean(beta2use{idx2plot}(:,ibin, :, :),1)), squeeze(std(beta2use{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % amplitude se
    tmpstd = squeeze(std(beta2use_mean{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_beta(ibin) x_beta(ibin)],[tmpBeta(ibin)-tmpstd tmpBeta(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot(x_beta(ibin), tmpBeta(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
end


xlim([cfg.xlim.beta])
ylim([cfg.ylim.beta])
set(gca,'fontsize',fSet.plFontsize,'ydir','reverse')
xlabel('Time from response (ms)','fontsize',fSet.axLabFontsize)
ylabel({'Power (dB)'},'fontsize',fSet.axLabFontsize)
% ydir('reverse')

% text(xBox(2)*0.99, yBox(1)*1.02, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)
text(xBox(1)*0.99, yBox(2)*0.7, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% beta slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hsub4 = subtightplot(nrow, ncol, [8 ], subplotgap2);
% figInit('ax');
axes('position',[0.61 0.88 0.14 0.07]) ; % inset
% axes('position',[0.81 0.12 0.12 0.09]) ; % inset
figInit('ax',[],{'fontsize',10});
set(gca,'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4])
hold on

hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'preRespBeta_slope',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(beta_pre_response_slope{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpBeta_slope(ibin)-tmpstd tmpBeta_slope(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpBeta_slope(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpBeta_slope)*1.15 max(tmpBeta_slope)*0.75])

YTICK = ax2.YTick;
ax2.YTick = linspace(min(YTICK), max(YTICK), 3);
% text(5, max(tmpBeta_slope)*0.85, '*','fontsize',35)

text(5, min(mean(beta_pre_response_slope{idx2plot})) * 1.1, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use,'xticklabel',[],'ydir','reverse')
% xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('LHB Slope', 'fontsize',fSet.axLabFontsize)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2 combined

%%% general settings

% x and y values where to plot N2 latency and amplitude
x_amplitude = [30:25:130];
y_latency   = [1.5:-0.25:0.5];

n2_xlim = [-20 500];
n2_ylim = [-3 1.9];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2, timecourse and latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hsub1 = subtightplot(nrow, ncol, [9:14 17:22 ], subplotgap);
figInit('ax');
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
    plot([tmpN2_lat(ibin) tmpN2_lat(ibin)], n2_ylim,'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
tmpN2_amp = mean(N2_amplitude2plot);
for ibin = 1:nbin2use
    plot(n2_xlim,[tmpN2_amp(ibin) tmpN2_amp(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end

% get stats for latency
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'N2i_latency',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
% plot model fit, latency
if ~isempty(mfit)
    h = plot(mfit.mean, y_latency, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], [y_latency flip(y_latency)] , h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

% draw box and print stats
xBox = [min(tmpN2_lat)-min(tmpN2_lat)*0.05 max(tmpN2_lat)+min(tmpN2_lat)*0.05];
yBox = [min(y_latency)-min(y_latency)*0.35 max(y_latency)+min(y_latency)*0.35];

text(max(xBox)*0.95, max(yBox)*1.2, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([-200 xBox(1)], [0 0] ,' k','linewidth',1)
plot([xBox(2) 800], [0 0] ,' k','linewidth',1)

% get stats, amplitude
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'N2i_amplitude',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(x_amplitude, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
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

text(min(xBox)*1.02, min(yBox)*1.1, {[betaString ' ' pstring_chi]}, 'fontsize',fSet.plFontsize)

clear h h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2
    h(ibin) = boundedline(cfg.t, squeeze(mean(N2_2plot(:,ibin, :),1)), squeeze(std(N2_2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
    set(h(ibin),'linewidth',fSet.plLineWidth,'LineStyle','--')
    
    % se
    tmpN2_lat_std = std(N2_latency2plot(:,ibin))/sqrt(nSub);
    plot([tmpN2_lat(ibin)-tmpN2_lat_std tmpN2_lat(ibin)+tmpN2_lat_std], [y_latency(ibin) y_latency(ibin)],'linewidth',fSet.plLineWidth_in,'Color','k')
    h2(ibin) = plot(tmpN2_lat(ibin), y_latency(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
    tmpN2_amp_std = std(N2_amplitude2plot(:,ibin))/sqrt(nSub);
    plot([x_amplitude(ibin) x_amplitude(ibin)], [tmpN2_amp(ibin)-tmpN2_amp_std tmpN2_amp(ibin)+tmpN2_amp_std], 'linewidth',fSet.plLineWidth_in,'Color','k')
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
    plot([tmpN2_lat(ibin) tmpN2_lat(ibin)], n2_ylim,'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
tmpN2_amp = mean(N2_amplitude2plot);
for ibin = 1:nbin2use
    plot(n2_xlim,[tmpN2_amp(ibin) tmpN2_amp(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end

% get stats for latency
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'N2c_latency',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
% plot model fit, latency
if ~isempty(mfit)
    h = plot(mfit.mean, y_latency, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], [y_latency flip(y_latency)] , h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

% draw box and print stats
xBox = [min(tmpN2_lat)-min(tmpN2_lat)*0.05 max(tmpN2_lat)+min(tmpN2_lat)*0.05];
yBox = [min(y_latency)-min(y_latency)*0.35 max(y_latency)+min(y_latency)*0.35];

text(max(xBox)*0.95, max(yBox)*1.2, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([-200 xBox(1)], [0 0] ,' k','linewidth',1)
plot([xBox(2) 800], [0 0] ,' k','linewidth',1)

% get stats, amplitude
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'N2c_amplitude',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(x_amplitude, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
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

text(min(xBox)*1.02, min(yBox)*1.1, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

clear h h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2
    h(ibin) = boundedline(cfg.t, squeeze(mean(N2_2plot(:,ibin, :),1)), squeeze(std(N2_2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % se
    tmpN2_lat_std = std(N2_latency2plot(:,ibin))/sqrt(nSub);
    plot([tmpN2_lat(ibin)-tmpN2_lat_std tmpN2_lat(ibin)+tmpN2_lat_std], [y_latency(ibin) y_latency(ibin)],'linewidth',fSet.plLineWidth_in,'Color','k')
    h2(ibin) = plot(tmpN2_lat(ibin), y_latency(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
    tmpN2_amp_std = std(N2_amplitude2plot(:,ibin))/sqrt(nSub);
    plot([x_amplitude(ibin) x_amplitude(ibin)], [tmpN2_amp(ibin)-tmpN2_amp_std tmpN2_amp(ibin)+tmpN2_amp_std], 'linewidth',fSet.plLineWidth_in,'Color','k')
    h2(ibin) = plot(x_amplitude(ibin),tmpN2_amp(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
end


%%% general settings
xlim(n2_xlim)
ylim(n2_ylim)

xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
ylabel({'Amplitude (\muV)'},'fontsize',fSet.axLabFontsize);

%%%%%%%%%%%%%%%%%%%
%%% ITPC, N2c
%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [15 16], subplotgap2);
figInit('ax');
hold on
% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'N2c_ITPC',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(N2c_ITPC_bar{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(N2c_ITPC_bar{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.8 1.2])

text(1.5, min(tmpmean) * 0.86, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('N2c ITPC', 'fontsize',fSet.axLabFontsize)




%%%%%%%%%%%%%%%%%%%
%%% ITPC, N2i
%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [23 24], subplotgap2);
figInit('ax');
hold on
% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'N2i_ITPC',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(N2i_ITPC_bar{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(N2i_ITPC_bar{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.8 1.2])

text(1.5, min(tmpmean) * 0.86, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('N2i ITPC', 'fontsize',fSet.axLabFontsize)

set(gcf,'renderer','painters')

saveFigName = [bin2use '_' bintype fileExt '_otherVar'];
for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary figure N2 ITPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fH, fSet] = figInit('fig', 3, {'height', 2.5/5; 'width',6/5});

nrow = 2;
ncol = 8;
subplotgap = [0.09 0.09];

%%%% ITPC
%%% imagesc plot
plotOrder = [1 2; 9 10];
for iplot = 1:2 % N2c and N2i
    
    switch iplot
        case 1
            %%% N2c
            N2_2use = N2c_ITPC;
            N2_2use_band = N2c_ITPC_band;
            N2_2use_bar = N2c_ITPC_bar;
            
            t2test      = [200 400];
        case 2
            %%% N2i
            N2_2use = N2i_ITPC;
            N2_2use_band = N2i_ITPC_band;
            N2_2use_bar = N2i_ITPC_bar;
            t2test      = [250 450];
    end
    
    hsub3 = subtightplot(nrow, ncol, plotOrder(iplot, :), subplotgap);
    figInit('ax');
    hold on
    
    colormap(hsub3,'Default')
    
    f2test      = [0.1 4]; %
    
    ttSPG = SPG_times-abs(t(1)/1000);
    ttSPG = ttSPG * 1000;
    maplimits = [min(min(min(squeeze(mean(N2_2use{idx2plot},1))))) max(max(max(squeeze(mean(N2_2use{idx2plot},1)))))];
    maplimits = [0 max(max(max(squeeze(mean(N2_2use{idx2plot},1)))))];
    
    contourf(ttSPG, SPG_freq, squeeze(mean(mean(N2_2use{idx2plot},1),2))',20,'LineColor','none');
    set(gca,'clim',maplimits)
    ylim([0 18])
    axis xy
    xlim([-50 600])
    hold on
    
    h = rectangle('Position',[t2test(1) f2test(1) diff(t2test) diff(f2test)]);%, '-w','linewidth',3)
    h.LineWidth = fSet.plLineWidth;
    h.EdgeColor = [1 1 1];
    
    ylabel('Frequency (Hz)','fontsize',fSet.axLabFontsize)
%     if iplot==2
    xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize)
%     end
    figInit('ax');
    set(gca,'color','white','linewidth',2)
    B=colorbar;
    % set(B, 'Position', [.345 .375 .01981 .1], 'Limits', maplimits)
    xpos = hsub3.Position(3) + hsub3.Position(1) - .01981 - 0.01;
    ypos = hsub3.Position(4) + hsub3.Position(2) - .1 - 0.02;
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
            N2_2use = N2c_ITPC;
            N2_2use_band = N2c_ITPC_band;
            N2_2use_bar = N2c_ITPC_bar;
            stats_label = 'N2c_ITPC';
            pupilIdx = idx_resp;
            t2test      = [200 400];
        case 2
            %%% N2i, pupil resp
            N2_2use = N2i_ITPC;
            N2_2use_band = N2i_ITPC_band;
            N2_2use_bar = N2i_ITPC_bar;
            stats_label = 'N2i_ITPC';
            pupilIdx = idx_resp;
            t2test      = [250 450];
        case 3
            %%% N2c, pupil baseline
            N2_2use = N2c_ITPC;
            N2_2use_band = N2c_ITPC_band;
            N2_2use_bar = N2c_ITPC_bar;
            
            stats_label = 'N2c_ITPC';
            pupilIdx = idx_BL_lp;
            t2test      = [200 400];
        case 4
            %%% N2i, pupil baseline
            N2_2use = N2i_ITPC;
            N2_2use_band = N2i_ITPC_band;
            N2_2use_bar = N2i_ITPC_bar;
            stats_label = 'N2i_ITPC';
            pupilIdx = idx_BL_lp;
            t2test      = [250 450];
    end
    
    hsubplot = subtightplot(nrow, ncol, plotOrder(iplot, :), subplotgap);
    figInit('ax');
    set(gca,'color','white')
    
    x_ITPC = [-80:35:65]+15;
    
    ttSPG = SPG_times-abs(t(1)/1000);
    ttSPG = ttSPG*1000;
    clear LEGEND
    hold on
    
    tmpN2 = nanmean(N2_2use_bar{pupilIdx});
    for ibin = 1:nbin2use
        plot([-200 700],[tmpN2(ibin) tmpN2(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
    end
    plot([t2test(1) t2test(2)], [0.11 0.11], 'k', 'linewidth', 4)
    
    % get stats
    bin2use = allBin2use{pupilIdx};
    filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
    [stats, mfit] = loadStatsR(filename, stats_label, nSub, nbin2use);
    [pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
    switch stats.model
        case 'Linear'
            betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
            model_color = fSet.colors(2,:);
        case 'Quadratic'
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
            model_color = fSet.colors(1,:);
        otherwise
            betaString = [];
    end
    if ~isempty(mfit)
        %     plot model fit
        h = plot(x_ITPC, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
        h = patch([x_ITPC flip(x_ITPC)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
        h.FaceAlpha = 0.3;
        h.EdgeAlpha = 0.3;
        h.EdgeColor = [h.FaceColor];
    end
    yBox = [min(tmpN2)-min(tmpN2)*0.07 max(tmpN2)+min(tmpN2)*0.07];
    xBox = [min(x_ITPC)-35 max(x_ITPC)+35];
    text(10, max(yBox)*1.04, {pstring_chi}, 'fontsize',fSet.plFontsize)
    
    h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
    h.LineWidth = 1;
    h.LineStyle = '--';
    
    
    clear h h2 LEGEND
    for ibin = 1:nbin2use
        h(ibin) = boundedline(ttSPG, squeeze(mean(N2_2use_band{pupilIdx}(:,ibin, :, :),1)), squeeze(std(N2_2use_band{pupilIdx}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
        
        set(h(ibin),'linewidth',fSet.plLineWidth)
        
        % se
        tmpstd = squeeze(std(N2_2use_band{pupilIdx}(:,ibin, :, :)))/sqrt(nSub);
        plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpN2(ibin)-tmpstd tmpN2(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
        
        % plotMarker
        h2(ibin) = plot(x_ITPC(ibin), tmpN2(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
        LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    end
        
    plot([0 0], [0 yBox(1)] ,' k','linewidth',1)
    plot([0 0], [yBox(2) 0.5] ,' k','linewidth',1)
    
    xlim([-130 550])
    ylim([0.1 0.33])
    % axis square
    xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize)
    ylabel(stats_label, 'fontsize', fSet.axLabFontsize)
    
end


% clear h h2 LEGEND
% for ibin = 1:nbin2use
%     h(ibin) = boundedline(ttSPG, squeeze(mean(N2_2use_band{idx2plot}(:,ibin, :, :),1)), squeeze(std(N2_2use_band{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
%     
%     set(h(ibin),'linewidth',fSet.plLineWidth)
%     
%     % se
%     tmpstd = squeeze(std(N2_2use_band{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
%     plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpN2(ibin)-tmpstd tmpN2(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
%     
%     % plotMarker
%     h2(ibin) = plot(x_ITPC(ibin), tmpN2(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
%     LEGEND{ibin} = ['Bin ' num2str(ibin) ];
% end
% 
% 
% % plot([0 0], [0 0.5] ,'k','linewidth',1)
% 
% plot([0 0], [0 yBox(1)] ,' k','linewidth',1)
% plot([0 0], [yBox(2) 0.5] ,' k','linewidth',1)
% 
% xlim([-130 550])
% ylim([0.1 0.27])
% % axis square
% xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize)
% ylabel('N2c ITPC','fontsize',fSet.axLabFontsize)


saveFigName = [bin2use '_' bintype fileExt '_N2_ITPC'];
for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% control, fig 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fH, fSet] = figInit('fig', 4, {'height', 3.5/5; 'width',6/5});

nrow = 3;
ncol = 4;
subplotgap = [0.1 0.065];
subplotgap2 = [0.1 0.065];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2i, timecourse and latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N2_2plot = squeeze(N2i{idx2plot});

n2_xlim = [-50 500];
n2_ylim = [-1.7 0.75];


hsub1 = subtightplot(nrow, ncol, [1 2], subplotgap);
figInit('ax');
hold on

plot([0 0], [-4 1] ,' k','linewidth',1)

% get stats
twinsize = 100;
tstep = 10;
trange = n2_xlim(1):tstep:n2_xlim(2);
p_N2i = NaN(1,length(trange));
for istep = 1:length(trange)
    
    twin = trange(istep) + [-twinsize/2 twinsize/2];
    twinidx = (twin(1) <= cfg.t & twin(2) >= cfg.t);
    
    N2_lme = squeeze(mean(N2_2plot(:,:,twinidx),3));
    
    lme_table = table(plevelSub(:), plevelBin(:), N2_lme(:), 'VariableNames',{'Subject','Bin','Y'});
    [statsreport, ~, ~, ~, STATS, model] = fitlme_pupil_singleVar(lme_table);
    p_N2i(istep) = statsreport(3);
    if strcmpi(model,'linear')
%         p_N2i(istep) = statsreport(3);
    elseif strcmpi(model,'quadratic')
        p_N2i(istep) = 999;
    end
    
end

signbin = false(length(p_N2i),1);
[pt, p_masked] = fdr(p_N2i,0.05);
signbin(p_N2i<pt) = true;


clear h h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2i
    h(ibin) = boundedline(cfg.t, squeeze(mean(N2_2plot(:,ibin, :),1)), squeeze(std(N2_2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
    set(h(ibin),'linewidth',fSet.plLineWidth)
end
h2 = plotBrokenVector(trange, -1.5, signbin);
h2.LineWidth = 3;

xlim(n2_xlim)
ylim(n2_ylim)

xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
% ylabel({'N2c \muV'},'fontsize',fSet.plFontsize)
ylabel({'Amplitude (\muV)'},'fontsize',fSet.axLabFontsize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2i, pupil response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [3], subplotgap);
figInit('ax');
hold on

tmpPupil = mean(pupil_lp_RT_neg200_200{idx_N2i});

% get stats
bin2use = allBin2use{idx_N2i};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'plevelPupil_RT200',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

if ~isempty(mfit)
    
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

text(0.5, min(tmpPupil)*0.85, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

for ibin = 1:nbin2use
    % se
    tmpPupil_std = squeeze(std(pupil_lp_RT_neg200_200{idx_N2i}(:,ibin)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpPupil(ibin)-tmpPupil_std tmpPupil(ibin)+tmpPupil_std],'linewidth',fSet.plLineWidth_in,'color','k');
    % plotMarker
    h2(ibin) = plot((ibin), tmpPupil(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize, 'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpPupil)*0.75 max(tmpPupil)*1.15])

set(gca,'xtick',1:nbin2use)
xlabel('N2i amplitude bin','fontsize',fSet.axLabFontsize)
ylabel('Pupil response','fontsize',fSet.axLabFontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% development of RT across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data2plot = RT_window;
% data2plot = RTcv_window;

data_min1   = {squeeze(data2plot{idx2plot}(:,:,:,5))};
data_0      = {squeeze(data2plot{idx2plot}(:,:,:,6))};
data_plus1  = {squeeze(data2plot{idx2plot}(:,:,:,7))};

subplotgap = [0.1 0.058];

allx = -5:5;
x2plot = -5:2;
x2plotIdx = ismember(allx, x2plot);
hsub = subtightplot(nrow, ncol, [9 10 ], subplotgap);
figInit('ax');

cfg.xlim     = [-5 2.5];
cfg.ylim     = [-0.2 0.2];

hold on

for ibin = 1:nbin2use
    A = boundedline(x2plot,squeeze(mean(data2plot{idx2plot}(:,ibin,:,x2plotIdx),1)),squeeze(std(data2plot{idx2plot}(:,ibin,:,x2plotIdx),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(A,'linewidth',fSet.plLineWidth_in)
end
xlim(cfg.xlim)
ylim(cfg.ylim)
set(gca,'xtick',-5:3)
xlabel('Trial relative to pupil response','fontsize',fSet.axLabFontsize);
ylabel('RT (Z-score)','fontsize',fSet.axLabFontsize);

clear P
bins2compare = [1 2 3 4];
RT2compare = [];
for itrial = 1:length(find(x2plotIdx))
    RT2compare = [mean(data2plot{idx2plot}(:,bins2compare,:,itrial),2), data2plot{idx2plot}(:,[5],:,itrial)];
    
    [H,P(itrial,1),CI,STATS] = ttest(RT2compare(:,1), RT2compare(:,2));
end

if 1
    [pt,p_masked] = fdr(P,0.05);
else
    pt = 0.05;
end
P(~p_masked) = 1;
[pstring_chi,starstring] = getSignificanceStrings(P, 0, 1);

YLIM = get(gca,'ylim');
for iP = 1:length(P)
    text(x2plot(iP), YLIM(2)*0.95, starstring{iP},'fontsize',fSet.axLabFontsize,'HorizontalAlignment','center')
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT across trials, barplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hsub = subtightplot(nrow, ncol, [11 12 ], subplotgap);

figInit('ax');
hold on

RT2compare = [];
x2plot = 1:3:nbin2use*3;
clear P
for ibin = 1:nbin2use
    RT2compare = [data_min1{1}(:,ibin), data_0{1}(:,ibin), data_plus1{1}(:,ibin)];
    
    [H,P(ibin,1),CI,STATS] = ttest(RT2compare(:,1), RT2compare(:,2));
    [H,P(ibin,2),CI,STATS] = ttest(RT2compare(:,2), RT2compare(:,3));
    
    h = bar([x2plot(ibin):x2plot(ibin)+2], mean(RT2compare));
    h.FaceColor = fSet.colors(ibin,:);
    
    plot([x2plot(ibin) x2plot(ibin)], mean(RT2compare(:,1)) + [-std(RT2compare(:,1))/sqrt(nSub) std(RT2compare(:,1))/sqrt(nSub)],'linewidth',1,'color','k')
    plot([x2plot(ibin)+1 x2plot(ibin)+1], mean(RT2compare(:,2)) + [-std(RT2compare(:,2))/sqrt(nSub) std(RT2compare(:,2))/sqrt(nSub)],'linewidth',1,'color','k')
    plot([x2plot(ibin)+2 x2plot(ibin)+2], mean(RT2compare(:,3)) + [-std(RT2compare(:,3))/sqrt(nSub) std(RT2compare(:,3))/sqrt(nSub)],'linewidth',1,'color','k')   
end

if 1
    [pt,p_masked] = fdr(P,0.05);
else
    pt = 0.05;
end
P(~p_masked) = 1;

clear groups
for ibin = 1:nbin2use
    
    for iC = 1:2
        if p_masked(ibin,iC)
            groups{ibin,iC} = [x2plot(ibin) x2plot(ibin)+1] + (iC-1);
        end
    end
    
end
groups = groups(p_masked);
sigstar(groups(:), P(p_masked))

set(gca,'xtick',1:15)
xTicks = get(gca, 'xtick');

% set some extra markers to indicate trial index
YLIM = ylim;
text(14,YLIM(1)*0.6, {'Trial index'},'horizontalalignment','center')
text(13,YLIM(1)*0.8, {'-1'},'horizontalalignment','center')
text(14,YLIM(1)*0.8, {'0'},'horizontalalignment','center')
text(15,YLIM(1)*0.8, {'1'},'horizontalalignment','center')
A = get(gca);
plot([13 13],[YLIM(1) YLIM(1)+A.TickLength(1)],'k','linewidth',2)
plot([14 14],[YLIM(1) YLIM(1)+A.TickLength(1)],'k','linewidth',2)
plot([15 15],[YLIM(1) YLIM(1)+A.TickLength(1)],'k','linewidth',2)

set(gca,'xtick',[2:3:nbin2use*3],'xticklabel',1:nbin2use)
xlabel('Pupil bin','fontsize',fSet.axLabFontsize);


saveFigName = [bin2use '_' bintype fileExt '_fig4_control'];
for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end


%% Figure 1, copy for CD dataset.
% in same format as fig 4, so can be pasted together

[fH, fSet] = figInit('fig', 3, {'height', 3.5/5; 'width',6/5});

nrow = 3;
ncol = 4;
subplotgap = [0.1 0.065];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% paradigm figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subtightplot(nrow, ncol, 5, subplotgap)

% import paradimg png
fig2import = ['C:\Jochem\Dropbox\Monash\bigDots_st\figures\paradigm_test.png'];
paradigm = imread(fig2import);

image(paradigm)
axis off
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pupil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subtightplot(nrow, ncol, 6, subplotgap)
figInit('ax');

cfg.xlim.Pupil     = [-150 2200];
cfg.ylim.Pupil      = [-0.1 0.08];

hold on

xlim(cfg.xlim.Pupil)
ylim(cfg.ylim.Pupil)
for ibin = 1:nbin2use
    A = boundedline(cfg.t,squeeze(mean(pupil_lp{idx2plot}(:,ibin,:,:),1)),squeeze(std(pupil_lp{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(A,'linewidth',fSet.plLineWidth_in)
end
plot([0 0], ylim ,'k','linewidth',1)


xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
ylabel('Normalized pupil diameter','fontsize',fSet.axLabFontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subtightplot(nrow, ncol, 7, subplotgap)
figInit('ax');

RT = RT;

% idx2plot    = [find(idx_BL_bp) find(idx2plot)];
% idx2plot    = [find(idx2plot)];
% % idx2plot    = [find(idx_BL_lp) find(idx_BL_bp)];
% idx2plot = idx2plot2;

cfg.ylim.RT     = [1000 1200];
cfg.ylim.RT_CV  = [0.175 0.36];
cfg.ylabel.RT   = 'RT (ms)';

% [stats2report, meanFit, CIFit, SEFit, STATS] = fitlme_pupil_singleVar(MAT, 'RT', 0);
% [pstring_chi,starstring] = getSignificanceStrings(stats2report(3), 0, 1, '\it p ');
% betaString = ['\beta_{2} = ' num2str(STATS.Estimate, '%1.3f')];

bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'RT', nSub, nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B * -1, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

ax1 = gca;
hold on
clear hPlotMark

if ~isempty(mfit)
    % plot model fit
    h = plot(mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', fSet.colors(1,:));
    h = patch(([(1:5) flip(1:5)]), [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
% plot se
for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT{idx2plot}(:,ibin))-std(RT{idx2plot}(:,ibin))/sqrt(nSub) ...
        mean(RT{idx2plot}(:,ibin))+std(RT{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.plLineWidth_in);
end
% plot plotMarker
hPlotMark(ibinType) = plot(x, mean(RT{idx2plot}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(1,:),'MarkerFaceColor', fSet.colors(1,:), 'MarkerEdgeColor', 'k');

ylim([cfg.ylim.RT])
h = ylabel(cfg.ylabel.RT,'fontsize',fSet.axLabFontsize);
h.Color = fSet.colors(1,:);
set(gca,'xtick',x)
xlabel('Pupil Bin','fontsize',fSet.axLabFontsize);
% axis square

h = text(0.7, max(mean(RT{idx2plot}(:,1)))*1.03, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize, 'color', fSet.colors(1,:));

% text for RT_CV
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'RT_CV', nSub, nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B * -1, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

text(3.3,max(mean(RT{idx2plot}(:,1)))*1.03, {betaString, pstring_chi }, 'fontsize',fSet.plFontsize, 'color', fSet.colors(2,:))

% plot RT_CV
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k','XTickLabel',[]);
figInit('ax');
ax2.YLim = [cfg.ylim.RT_CV];
ax2.XColor = 'none';
linkaxes([ax1 ax2],'x');
hold on

if ~isempty(mfit)
    h = plot(mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', fSet.colors(2,:));
    h = patch([(1:5) flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT_CV{idx2plot}(:,ibin))-std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ...
        mean(RT_CV{idx2plot}(:,ibin))+std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.plLineWidth_in);
end

hPlotMark(ibinType) = plot(x, mean(RT_CV{idx2plot}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(2,:),'MarkerFaceColor', fSet.colors(2,:), 'MarkerEdgeColor', 'k');
h = ylabel('RT CV','fontsize',fSet.axLabFontsize);
h.Color = fSet.colors(2,:);
% axis square
subplot_ax1 = get(gca);

xlim([0.5 nbin2use+0.5])









saveFigName = [bin2use '_' bintype fileExt '_RT_CD'];

for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            %             print(gcf,['-d' figureFileType{ifiletype} 'c'],[paths.pop 'fig' filesep 'p_level' filesep saveFigName '.' figureFileType{ifiletype}])
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
            
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reduced - CPP (fig5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fH, fSet] = figInit('fig', 5, {'height', 2.5/5; 'width',6/5});

nrow = 2;
ncol = 5;
subplotgap = [0.12 0.058];
subplotgap2 = [0.12 0.06];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CPP onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [1], subplotgap2);
figInit('ax');

hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'CPP_onset',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(nanmean(CPP_onset{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(nanstd(CPP_onset{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.95 1.05])

text(1, min(tmpmean) * 0.965, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('CPP onset latency (ms)', 'fontsize',fSet.axLabFontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CPP slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [2], subplotgap2);
figInit('ax');

hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'CPP_slope2',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(CPPr_csd_slope{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(CPPr_csd_slope{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.8 1.2])

text(1, min(tmpmean) * 0.86, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('CPP build-up rate', 'fontsize',fSet.axLabFontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CPP amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [3], subplotgap2);
figInit('ax');

hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'CPPr_amplitude',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(CPPr_amplitude{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(CPPr_amplitude{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)].* [0.8 1.2] )%

text(1, min(tmpmean) * 0.86, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('CPP amplitude (\muV/m^{2})', 'fontsize',fSet.axLabFontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ITPC band

hsub = subtightplot(nrow, ncol, [4 5], subplotgap);
% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');


x_ITPC = [45:45:225];

ttSPG = SPG_times-abs(t(1)/1000);
ttSPG = ttSPG*1000;
clear LEGEND
hold on

t2test      = [300 550];

tmpCPP = nanmean(CPP_ITPC_bar{idx2plot});
for ibin = 1:nbin2use
    if setAlpha
        plot([-200 700],[tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
    else
        plot([-200 700],[tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:)]);
    end
end
plot([t2test(1) t2test(2)], [0.12 0.12], 'k', 'linewidth', 4)

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'CPP_ITPC',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(x_ITPC, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([x_ITPC flip(x_ITPC)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
clear h h2 LEGEND
for ibin = 1:nbin2use
    if setAlpha
        h(ibin) = boundedline(ttSPG, squeeze(mean(CPP_ITPC_band{idx2plot}(:,ibin, :, :),1)), squeeze(std(CPP_ITPC_band{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    else
        h(ibin) = boundedline(ttSPG, squeeze(mean(CPP_ITPC_band{idx2plot}(:,ibin, :, :),1)), squeeze(std(CPP_ITPC_band{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:));
    end
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % se
    tmpstd = squeeze(std(CPP_ITPC_band{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpCPP(ibin)-tmpstd tmpCPP(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    % plotMarker
    h2(ibin) = plot(x_ITPC(ibin), tmpCPP(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end

yBox = [min(tmpCPP)-min(tmpCPP)*0.07 max(tmpCPP)+min(tmpCPP)*0.07];
xBox = [min(x_ITPC)-min(x_ITPC)*0.6 max(x_ITPC)+min(x_ITPC)*0.6];
text(min(xBox)*0.95, max(yBox)*1.09, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([0 0], [0 0.5] ,'k','linewidth',1)
xlim([-50 700])
ylim([0.1 0.5])
% axis square
xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize)
ylabel({'CPP ITPC'},'fontsize',fSet.axLabFontsize)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveFigName = [bin2use '_' bintype fileExt '_redCPP'];
for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reduced - other var
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fH, fSet] = figInit('fig', 7, {'height', 2.5/5; 'width',6/5});

nrow = 2;
ncol = 9;
subplotgap = [0.12 0.058];
subplotgap2 = [0.12 0.07];

subplotorder = [...
    1 2;...
    3 7;...
    10 12;...
    13 14;...
    15 16;...
    15 16;...
    17 18;...
    ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [subplotorder(1,1):subplotorder(1,2)], subplotgap2);
figInit('ax');

hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'alpha',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
tmpmean = squeeze(mean(alpha{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(alpha{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.85 1.15])

text(1.5, min(tmpmean) * 0.9, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('\alpha power', 'fontsize',fSet.axLabFontsize)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub1 = subtightplot(nrow, ncol, [subplotorder(2,1):subplotorder(2,2)], subplotgap2);
figInit('ax');

cfg.xlim = [-530 80];
cfg.ylim = [0.52 0.76];
cfg.ylim = [-0.2 0.01];


x_beta = [-440:30:-320];

spectral_t = stft_timesr;

hold on

beta2use_mean = beta_base_response_amplitude;
beta2use = beta_base_response;
% beta2use_mean = beta_response_mean;
% beta2use = beta_response;

tmpBeta = mean(beta2use_mean{idx2plot});
tmpBeta_slope = mean(beta_pre_response_slope{idx2plot});
for ibin = 1:nbin2use
    plot([-600 100],[tmpBeta(ibin) tmpBeta(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
%
% plot(x_beta, tmpBeta,'linewidth',4,'color',[0.5 0.5 0.5 0.5])
plot([0 0], [cfg.ylim] ,'k','linewidth',1)

twin_bar    = [-130 -70];
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
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'preRespBeta_base',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    %     plot model fit
    h = plot(x_beta, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([x_beta flip(x_beta)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

clear h h2 LEGEND
for ibin = 1:nbin2use
    % plot line
    h(ibin) = boundedline(spectral_t, squeeze(mean(beta2use{idx2plot}(:,ibin, :, :),1)), squeeze(std(beta2use{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % amplitude se
    tmpstd = squeeze(std(beta2use_mean{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_beta(ibin) x_beta(ibin)],[tmpBeta(ibin)-tmpstd tmpBeta(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot(x_beta(ibin), tmpBeta(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
end


xlim([cfg.xlim])
ylim([cfg.ylim])
set(gca,'fontsize',fSet.plFontsize,'ydir','reverse')
xlabel('Time from response (ms)','fontsize',fSet.axLabFontsize)
ylabel({'Power (dB)'},'fontsize',fSet.axLabFontsize)
% ydir('reverse')

text(xBox(1)*0.99, yBox(2)*0.7, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% beta slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2 = axes('position',[0.37 0.83 0.14 0.09]) ; % inset
figInit('ax',[],{'fontsize',10});
set(gca,'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4])
hold on

hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'preRespBeta_slope',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.1e')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end
if ~isempty(mfit)
    
    % plot model fit
    h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(beta_pre_response_slope{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpBeta_slope(ibin)-tmpstd tmpBeta_slope(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpBeta_slope(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpBeta_slope)*1.15 max(tmpBeta_slope)*0.75])

YTICK = ax2.YTick;
ax2.YTick = linspace(min(YTICK), max(YTICK), 3);
% text(5, max(tmpBeta_slope)*0.85, '*','fontsize',35)

text(5, min(mean(beta_pre_response_slope{idx2plot})) * 1.1, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use,'xticklabel',[],'ydir','reverse')
% xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('LHB Slope', 'fontsize',fSet.axLabFontsize)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2, timecourse and latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [subplotorder(3,1):subplotorder(3,2)], subplotgap2);
figInit('ax');
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
    plot([tmpN2_lat(ibin) tmpN2_lat(ibin)], [-20 1],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
tmpN2_amp = mean(N2_amplitude2plot);
for ibin = 1:nbin2use
    plot([0 400],[tmpN2_amp(ibin) tmpN2_amp(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end


% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'N2c_latency',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

% plot model fit
if ~isempty(mfit)
    clear h
    h = plot(mfit.mean, y_latency, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
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

text(max(xBox)*1.02, max(yBox)*0.85, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

% zero lines
plot([-200 xBox(1)], [0 0] ,' k','linewidth',1)
plot([xBox(2) 800], [0 0] ,' k','linewidth',1)


% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'N2c_amplitude',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

% plot model fit
if ~isempty(mfit)
    clear hM
    hM = plot(x_amplitude, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
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

text(min(xBox)*1.02, max(yBox)*0.7, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

clear h1 h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2c
    h1 = boundedline(cfg.t, squeeze(mean(N2_2plot(:,ibin, :),1)), squeeze(std(N2_2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
    set(h1,'linewidth',fSet.plLineWidth)
    
    % se
    tmpN2_lat_std = std(N2_latency2plot(:,ibin))/sqrt(nSub);
    plot([tmpN2_lat(ibin)-tmpN2_lat_std tmpN2_lat(ibin)+tmpN2_lat_std], [y_latency(ibin) y_latency(ibin)],'linewidth',fSet.plLineWidth_in,'Color','k')
    tmpN2_amp_std = std(N2_amplitude2plot(:,ibin))/sqrt(nSub);
    plot([x_amplitude(ibin) x_amplitude(ibin)], [tmpN2_amp(ibin)-tmpN2_amp_std tmpN2_amp(ibin)+tmpN2_amp_std], 'linewidth',fSet.plLineWidth_in,'Color','k')
    
    
    plot(tmpN2_lat(ibin), y_latency(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    plot(x_amplitude(ibin),tmpN2_amp(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
end

xlim([-20 350])
ylim([-3.2 0.5])

xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
ylabel({'Amplitude (\muV)'},'fontsize',fSet.axLabFontsize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [subplotorder(4,1):subplotorder(4,2)], subplotgap2 + [0 0.02]);
figInit('ax');

hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];

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
    [stats, mfit] = loadStatsR(filename, stats2load, nSub, nbin2use);
    [pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
    switch stats.model
        case 'Linear'
            betaString = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
            model_color = fSet.colors(2,:);
        case 'Quadratic'
            betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
            model_color = fSet.colors(1,:);
        otherwise
            betaString = [];
    end
    if ~isempty(mfit)
        % plot model fit
        h = plot(1:nbin2use, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
        h = patch([1:nbin2use flip(1:nbin2use)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
        h.FaceAlpha = 0.3;
        h.EdgeAlpha = 0.3;
        h.EdgeColor = [h.FaceColor];
    end
    tmpmean = squeeze(mean(data2plot{idx2plot}));
    for ibin = 1:nbin2use
        %     slope
        tmpstd = squeeze(std(data2plot{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
        plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
        
        h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
        LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    end
    
    text(1, min(tmpmean) * LME_text{2}, {[LME_text{1}, betaString, pstring_chi]}, 'fontsize',fSet.plFontsize)
    clear data2plot
end
xlim([0 6])
YLIM = get(gca,'ylim');
ylim([YLIM .* [0.9 1.1]]);

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('N2 ITPC', 'fontsize',fSet.axLabFontsize)


saveFigName = [bin2use '_' bintype fileExt '_redOtherVar'];
for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end


%% reduced - other var. Extra CD plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fH, fSet] = figInit('fig', 6, {'height', 2.5/5; 'width',6/5});

nrow = 2;
ncol = 9;
subplotgap = [0.12 0.058];
subplotgap2 = [0.12 0.07];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pupil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subtightplot(nrow, ncol, [subplotorder(7,1):subplotorder(7,2)], subplotgap2)
figInit('ax');

cfg.xlim.Pupil     = [-150 2200];
cfg.ylim.Pupil      = [-0.07 0.07];

hold on

xlim(cfg.xlim.Pupil)
ylim(cfg.ylim.Pupil)
for ibin = 1:nbin2use
    A = boundedline(cfg.t,squeeze(mean(pupil_lp{idx2plot}(:,ibin,:,:),1)),squeeze(std(pupil_lp{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(A,'linewidth',fSet.plLineWidth_in)
end
plot([0 0], ylim ,'k','linewidth',1)


xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
ylabel('Normalized pupil diameter','fontsize',fSet.axLabFontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subtightplot(nrow, ncol, [subplotorder(8,1):subplotorder(8,2)], subplotgap2)
figInit('ax');
hold on
RT = RT;

cfg.ylim.RT     = [1000 1200];
cfg.ylim.RT_CV  = [0.175 0.32];
cfg.ylabel.RT   = 'RT (ms)';

bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];

clear betaString pstring_chi
for imod = 1:2
    
    if imod == 1
        forceMod = 'Quadratic';
    elseif imod == 2
        forceMod = 'Linear';
    end
    
    [stats, mfit] = loadStatsR(filename, 'RT', nSub, nbin2use, forceMod);
    [pstring_chi{imod},starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
    switch stats.model
        case 'Linear'
            betaString{imod} = ['\beta_{1} = ' num2str(stats.B, '%1.2f')];
            model_color = fSet.colors(2,:);
            %             model_color = [0.3 0.447 1];
        case 'Quadratic'
            betaString{imod} = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
            model_color = fSet.colors(1,:);
            %             model_color = [0.5 0.5 0.5];
        otherwise
            betaString{imod} = [];
    end
    
    
    if ~isempty(mfit)
        % plot model fit
        h = plot(mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
        h = patch(([(1:5) flip(1:5)]), [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
        h.FaceAlpha = 0.3;
        h.EdgeAlpha = 0.3;
        h.EdgeColor = [h.FaceColor];
    end
    
end
ax1 = gca;
hold on
clear hPlotMark
% plot se
for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT{idx2plot}(:,ibin))-std(RT{idx2plot}(:,ibin))/sqrt(nSub) ...
        mean(RT{idx2plot}(:,ibin))+std(RT{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.plLineWidth_in);
end
% plot plotMarker
hPlotMark(ibinType) = plot(1:nbin2use, mean(RT{idx2plot}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(1,:),'MarkerFaceColor', fSet.colors(1,:), 'MarkerEdgeColor', 'k');
ylim([cfg.ylim.RT])
h = ylabel(cfg.ylabel.RT,'fontsize',fSet.axLabFontsize);
h.Color = fSet.colors(1,:);
set(gca,'xtick',x)
xlabel('Pupil Bin','fontsize',fSet.axLabFontsize);
% axis square

h = text(0.7, max(mean(RT{idx2plot}(:,1)))*1.04, {betaString{:}, pstring_chi{:}}, 'fontsize',fSet.plFontsize, 'color', fSet.colors(1,:));
clear betaString pstring_chi

% text for RT_CV
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'RT_CV', nSub, nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
switch stats.model
    case 'Linear'
        betaString = ['\beta_{1} = ' num2str(stats.B * -1, '%1.2f')];
        model_color = fSet.colors(2,:);
    case 'Quadratic'
        betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];
        model_color = fSet.colors(1,:);
    otherwise
        betaString = [];
end

text(3.3,max(mean(RT{idx2plot}(:,1)))*1.04, {betaString, pstring_chi }, 'fontsize',fSet.plFontsize, 'color', fSet.colors(2,:))

% plot RT_CV
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k','XTickLabel',[]);
figInit('ax');
ax2.YLim = [cfg.ylim.RT_CV];
ax2.XColor = 'none';
linkaxes([ax1 ax2],'x');
hold on

if ~isempty(mfit)
    h = plot(mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', fSet.colors(2,:));
    h = patch([(1:5) flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(RT_CV{idx2plot}(:,ibin))-std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ...
        mean(RT_CV{idx2plot}(:,ibin))+std(RT_CV{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'color', 'k','linewidth',fSet.plLineWidth_in);
end

hPlotMark(ibinType) = plot(1:nbin2use, mean(RT_CV{idx2plot}),plotMarker{1},...
    'markersize', fSet.MarkerSize,...
    'color', fSet.colors(2,:),'MarkerFaceColor', fSet.colors(2,:), 'MarkerEdgeColor', 'k');
h = ylabel('RT CV','fontsize',fSet.axLabFontsize);
h.Color = fSet.colors(2,:);
% axis square
subplot_ax1 = get(gca);

xlim([0 nbin2use+1])

saveFigName = [bin2use '_' bintype fileExt '_RT_CD_BL'];

for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            %             print(gcf,['-d' figureFileType{ifiletype} 'c'],[paths.pop 'fig' filesep 'p_level' filesep saveFigName '.' figureFileType{ifiletype}])
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
            
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end

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
saveFigName = [bin2use '_' bintype fileExt '_CPP_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
% set(gca,'color','none')

export_fig([paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.png'],'-transparent')

%N2 - left
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(N2c_topo{idx2plot}(:,:,1,:),1),2)))) 1];
topoplot(squeeze(mean(mean(N2c_topo{idx2plot}(:,:,1,:),1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
saveFigName = [bin2use '_' bintype fileExt '_N2l_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
set(gca,'color','none')

export_fig([paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.png'],'-transparent')


%N2 - right
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(N2c_topo{idx2plot}(:,:,2,:),1),2)))) 1];
topoplot(squeeze(mean(mean(N2c_topo{idx2plot}(:,:,2,:),1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
saveFigName = [bin2use '_' bintype fileExt '_N2r_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
set(gca,'color','none')

export_fig([paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.png'],'-transparent')

%alpha
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(alpha_preTarget_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(alpha_preTarget_topo{idx2plot},1),2))))];
topoplot(squeeze(mean(mean(alpha_preTarget_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
saveFigName = [bin2use '_' bintype fileExt '_Alpha_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
set(gca,'color','none')

export_fig([paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.png'],'-transparent')

%beta
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(beta_base_pre_response_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(beta_base_pre_response_topo{idx2plot},1),2))))];
topoplot(squeeze(mean(mean(beta_base_pre_response_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);

saveFigName = [bin2use '_' bintype fileExt '_Beta_topo'];
A = get(gca);
A.Children(1).LineWidth = 4;
A.Children(2).LineWidth = 4;
A.Children(3).LineWidth = 4;
A.Children(4).LineWidth = 4;
axis xy
% set(gcf,'color','none')
set(gca,'color','none')

export_fig([paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.png'],'-transparent')

