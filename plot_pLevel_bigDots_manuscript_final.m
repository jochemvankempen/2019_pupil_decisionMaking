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

idx_alpha       = strcmpi(allBin2use, 'pretarget_alpha');
idx_alpha_asym  = strcmpi(allBin2use, 'pretarget_alpha_asym');

idx_N2i         = strcmpi(allBin2use, 'N2i_amplitude_regress_iti_side');

idx2plot   = idx_resp;
idx2plot   = idx_BL_lp;


%% plot Settings

chanlocs = readlocs('cap64.loc'); %biosemi
[~,plot_chans, exclude_chans] = getChannelName;

figureFileType = 'png';
% figureFileType = 'jpeg';
figureFileType = 'eps';

figureFileType = {'png','eps'};
clear cfg
cfg.t = t;
cfg.tr = tr;


% plotMarker = {'o','s','d'};
% plotMarker = {'o','d'};
plotMarker = {'o','s','d','^','h'};


%% Figure 1

nrow = 1;
ncol = 3;
subplotgap = [0.06 0.08];

% fig1 = figure(1);clf
[figHandle, fSet] = figInit('fig',1, {'height', 1/3; 'width', 6/5});
% [figHandle, fSet] = figInit('fig',1, {'height', 1/3});
set(gca,'linewidth',2)
% set(gcf,'position',[-1420 154 1092 981])

set(gca,'color','none','XColor','k','YColor','k')
set(gcf,'PaperPositionMode','auto')

set(gca,'fontsize',fSet.plFontsize)


switch allBin2use{idx2plot}
    case 'pupil_lp_baseline_regress_iti_side'
        text_y = 0.855;
    case 'pupil_lp_RT_neg200_200_regress_bl_iti_side'
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

cfg.xlim.Pupil     = [-150 1200];
cfg.ylim.Pupil      = [-0.2 0.45];
cfg.ylim.Pupil      = [-0.05 0.09];

hold on

clear LEGEND
for ibin = 1:nbin2use
    plot(cfg.t,squeeze(mean(squeeze(pupil_lp{idx2plot}(:,ibin,:,:)),1)),'linewidth',fSet.plLineWidth_in,'Color',fSet.colors(ibin,:));
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
    A = boundedline(cfg.t,squeeze(mean(pupil_lp{idx2plot}(:,ibin,:,:),1)),squeeze(std(pupil_lp{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(A,'linewidth',fSet.plLineWidth_in)
end
plot([0 0], ylim ,'k','linewidth',1)


xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
ylabel('Normalized pupil diameter','fontsize',fSet.axLabFontsize);
axis square



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subtightplot(nrow, ncol, 3, subplotgap)
figInit('ax');

RT = RT;

% idx2plot    = [find(idx_BL_bp) find(idx2plot)];
% idx2plot    = [find(idx2plot)];
% % idx2plot    = [find(idx_BL_lp) find(idx_BL_bp)];
% idx2plot = idx2plot2;

cfg.ylim.RT     = [470 620];
cfg.ylim.RT_CV  = [0.195 0.28];
cfg.ylabel.RT   = 'RT (ms)';

% [stats2report, meanFit, CIFit, SEFit, STATS] = fitlme_pupil_singleVar(MAT, 'RT', 0);
% [pstring_chi,starstring] = getSignificanceStrings(stats2report(3), 0, 1, '\it p ');
% betaString = ['\beta_{2} = ' num2str(STATS.Estimate, '%1.3f')];

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
subplotgap = [0.07 0.11];
subplotgap = [0.06 0.1];

% idx2plot    = [find(idx2plot)];
% idx2plot = idx2plot2;
figsize = 4/5;

[figHandle, fSet] = figInit('fig', 2, {'width',6/5 ;'height', figsize});

nsubplot = 3; % onset, response, ITPC
nrow = 3;
ncol = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% stim locked/ onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub1 = subtightplot(nrow, ncol, [1:5], subplotgap);
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
    tmpCPPstd = nanstd(CPP_onset2plot(:,ibin))/sqrt(nSub);
    plot([tmpCPP(ibin)-tmpCPPstd tmpCPP(ibin)+tmpCPPstd], [y_onset(ibin) y_onset(ibin)],'linewidth',fSet.plLineWidth_in,'color','k');
    
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

plot([0 0], [-2 30] ,'-k','linewidth',1)
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
if 1
    hsub2 = subtightplot(nrow, ncol, [6:10], subplotgap);
    figInit('ax');
    hold on
    
    CPP2plot        = squeeze(CPPr_csd{idx2plot});
    CPP_amp2plot    = squeeze(CPPr_amplitude{idx2plot});
    CPP_slope2plot  = squeeze(CPPr_csd_slope{idx2plot});
    
    % location to plot plotMarker CPP threshold
    x_amp = [-320:15:-250]+140;
    
    t2test_amp      = [-100 0]; % time threshold is computed on
    t2test_slope    = [-100 0]; % time slope is computed on
    
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
    plot([t2test_amp(1) t2test_amp(2)], [2 2], 'k', 'linewidth', 4)
    % plot([t2test_slope(1) t2test_slope(2)], [1.5 1.5], 'color',[0.5 0.5 0.5], 'linewidth', 4)
    
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
        tmpCPPstd = std(CPP_amp2plot(:,ibin))/sqrt(nSub);
        plot([x_amp(ibin) x_amp(ibin)], [tmpCPP(ibin)-tmpCPPstd tmpCPP(ibin)+tmpCPPstd],'linewidth',fSet.plLineWidth_in,'color','k');
        h2(ibin) = plot(x_amp(ibin), tmpCPP(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
        LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    end
    
    plot([0 0], [-0 40] ,'-k','linewidth',1)
    plot([-100 600], [0 0] ,'-k','linewidth',1)
    xlim([-350 50])
    ylim([-0 30])
    
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
    ylabel({'CPP \muV/m^{2}'},'fontsize',fSet.axLabFontsize)
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub3 = subtightplot(nrow, ncol, [11 12], subplotgap);
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

hsub4 = subtightplot(nrow, ncol, [13:15], subplotgap);
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
    tmpCPPstd = squeeze(std(CPP_ITPC_band{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpCPP(ibin)-tmpCPPstd tmpCPP(ibin)+tmpCPPstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
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
ylabel({'ITPC'},'fontsize',fSet.axLabFontsize)


saveFigName = [bin2use '_' bintype fileExt '_CPP'];
for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end


%% CPP resp ALL (for some reason, the figure is not saved as vector graphics when the response subplot is in it, so here it is plotted separately in the same proportions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setAlpha = 1;
subplotgap = [0.07 0.11];

% idx2plot    = [find(idx2plot)];
% idx2plot = idx2plot2;

% [figHandle, fSet] = figInit('fig', 2, {'width',1 ;'height', 1});
[figHandle, fSet] = figInit('fig', 2, {'width',6/5 ;'height', 4/5});

nsubplot = 3; % onset, response, ITPC
nrow = 3;
ncol = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% stim locked/ onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% resp locked/ threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hsub2 = subtightplot(nrow, ncol, [6:10], subplotgap);
    figInit('ax');
    hold on
    
    CPP2plot        = squeeze(CPPr_csd{idx2plot});
    CPP_amp2plot    = squeeze(CPPr_amplitude{idx2plot});
    CPP_slope2plot  = squeeze(CPPr_csd_slope{idx2plot});
    
    % location to plot plotMarker CPP threshold
    x_amp = [-320:10:-280]+140;
    
    t2test_amp      = [-100 0]; % time threshold is computed on
    t2test_slope    = [-100 0]; % time slope is computed on
    
    % plot lines indicating threshold
    tmpCPP = mean(CPP_amp2plot);
    for ibin = 1:nbin2use
        if setAlpha
            plot([max(x_amp)+min(x_amp)*0.28 100], [tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
        else
            plot([max(x_amp)+min(x_amp)*0.28 100], [tmpCPP(ibin) tmpCPP(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:)]);
        end
    end
    
    % plot times amplitude and slope are computed over
    plot([t2test_amp(1) t2test_amp(2)], [2 2], 'k', 'linewidth', 4)
    % plot([t2test_slope(1) t2test_slope(2)], [1.5 1.5], 'color',[0.5 0.5 0.5], 'linewidth', 4)
    
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
        tmpCPPstd = std(CPP_amp2plot(:,ibin))/sqrt(nSub);
        plot([x_amp(ibin) x_amp(ibin)], [tmpCPP(ibin)-tmpCPPstd tmpCPP(ibin)+tmpCPPstd],'linewidth',fSet.plLineWidth_in,'color','k');
        h2(ibin) = plot(x_amp(ibin), tmpCPP(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
        LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    end
    
    plot([0 0], [-0 40] ,'-k','linewidth',1)
    plot([-100 600], [0 0] ,'-k','linewidth',1)
    xlim([-300 20])
    ylim([-0 30])
    
    yBox = [min(tmpCPP)-min(tmpCPP)*0.1 max(tmpCPP)+min(tmpCPP)*0.1];
    xBox = [min(x_amp)-min(x_amp)*0.28 max(x_amp)+min(x_amp)*0.28];
    
    h = rectangle('Position',[xBox(2) yBox(1) abs(diff(xBox)) diff(yBox)]);
    h.LineWidth = 1;
    h.LineStyle = '--';
    
    % get stats
    bin2use = allBin2use{idx2plot};
    filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
    [stats, mfit] = loadStatsR(filename, 'CPPr_amplitude',nSub,nbin2use);
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
    text(min(xBox)*0.99, max(yBox)*1.07, {pstring_chi}, 'fontsize',fSet.plFontsize)
    
    xlabel('Time from response (ms)','fontsize',fSet.axLabFontsize)
    ylabel({'CPP \muV/m^{2}'},'fontsize',fSet.plFontsize)
%     ylabel({'CSD ()'},'fontsize',fSet.axLabFontsize)
    %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% resp locked/ slope
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    axes('position',[0.14 0.48 0.2 0.09]) ; % inset
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



saveFigName = [bin2use '_' bintype fileExt '_CPPresp'];
for ifiletype = 1:length(figureFileType)
    switch figureFileType{ifiletype}
        case 'eps'
            saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
        otherwise
            print(gcf,['-d' figureFileType{ifiletype} ],[paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.' figureFileType{ifiletype}])
    end
end


%%%%%% 
%% fig 3 without control
%%%%%%%%%%%%%%%%%%%%%


[fH, fSet] = figInit('fig', 3, {'height', 2.5/5; 'width',6/5});

nrow = 2;
ncol = 4;
subplotgap = [0.12 0.058];
subplotgap2 = [0.12 0.09];
subplotgap3 = [0.12 0.03];
% subplotgap3 = [0.12 0.12];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% alpha 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hsub = subtightplot(nrow, ncol, [1], subplotgap2);
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
hsub = subtightplot(nrow, ncol, [2], subplotgap2);

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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2, timecourse and latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N2c2plot = squeeze(N2c{idx2plot});
N2c_latency2plot = squeeze(N2c_latency{idx2plot});
N2c_amplitude2plot = squeeze(N2c_amplitude{idx2plot});

alignN2c = 0;
if alignN2c % align N2c to a N2c latency value
    bin2align = 3;
    tWin = [-0.500 0.400];% chop out this time window
    tWin = tWin / (1/fs);
    alignN2c = zeros(size(N2c2plot,1),size(N2c2plot,2),length(tWin(1):tWin(2)));
    
    for isub = 1:nSub
        
%         t2align = N2c_latency2plot(isub, bin2align);
        t2align = mean(N2c_latency2plot(isub, :));
        [~,tIdx] = min(abs(cfg.t - t2align));
        
        alignN2c(isub, :, :) = N2c2plot(isub,:,(tIdx + tWin(1)):(tIdx + tWin(2)));
        
%         N2c_latency2plot(isub,:) = N2c_latency2plot(isub,:) - N2c_latency2plot(isub, bin2align);
        N2c_latency2plot(isub,:) = N2c_latency2plot(isub,:) - mean(N2c_latency2plot(isub, :));
    end
end



% x and y values where to plot N2 latency and amplitude
y_latency = [0.3:-0.3:-0.9];
x_amplitude = [40:20:120];

% times N2 values have been computed on
t2test_amp      = 266 + [-50 50];
t2test_lat      = [150 400];

hsub1 = subtightplot(nrow, ncol, [3 4], subplotgap);
figInit('ax');

hold on

plot([0 0], [-4 1] ,' k','linewidth',1)

plot(t2test_amp, [-2.85 -2.85] ,' k','linewidth',4)

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
h = plot(mfit.mean, y_latency, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
h = patch([mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], [y_latency flip(y_latency)] , h.Color);
h.FaceAlpha = 0.3;
h.EdgeAlpha = 0.3;
h.EdgeColor = [h.FaceColor];
end

% plot lines to indicate N2c lat and amp
tmpN2c_lat = mean(N2c_latency2plot);
for ibin = 1:nbin2use
    plot([tmpN2c_lat(ibin) tmpN2c_lat(ibin)], [-20 1],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
tmpN2c_amp = mean(N2c_amplitude2plot);
for ibin = 1:nbin2use
    plot([0 400],[tmpN2c_amp(ibin) tmpN2c_amp(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end

clear h h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2c
    if ~alignN2c
        h(ibin) = boundedline(cfg.t, squeeze(mean(N2c2plot(:,ibin, :),1)), squeeze(std(N2c2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
    else
%         h(ibin) = boundedline(tWin(1):tWin(2), squeeze(mean(alignN2c(:,ibin, :),1)), squeeze(std(alignN2c(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
        h(ibin) = plot(tWin(1):tWin(2), squeeze(mean(alignN2c(:,ibin, :),1)),'color',fSet.colors(ibin,:));
    end
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % se
    tmpN2c_lat_std = std(N2c_latency2plot(:,ibin))/sqrt(nSub);
    plot([tmpN2c_lat(ibin)-tmpN2c_lat_std tmpN2c_lat(ibin)+tmpN2c_lat_std], [y_latency(ibin) y_latency(ibin)],'linewidth',fSet.plLineWidth_in,'Color','k')
    h2(ibin) = plot(tmpN2c_lat(ibin), y_latency(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    %     plot([tmpCPP(ibin)-tmpCPPstd tmpCPP(ibin)+tmpCPPstd], [y_onset(ibin) y_onset(ibin)],'linewidth',4,'color',[fSet.colors(ibin,:)]);
    
    tmpN2c_amp_std = std(N2c_amplitude2plot(:,ibin))/sqrt(nSub);
    plot([x_amplitude(ibin) x_amplitude(ibin)], [tmpN2c_amp(ibin)-tmpN2c_amp_std tmpN2c_amp(ibin)+tmpN2c_amp_std], 'linewidth',fSet.plLineWidth_in,'Color','k')
    h2(ibin) = plot(x_amplitude(ibin),tmpN2c_amp(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
    
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    
end
% [hLeg,icons,plots] = legend([h2], LEGEND,'location','southwest','fontsize',fSet.axLabFontsize , 'box','off');
% 
% Pos = hLeg.Position;
% hLeg.Position(1) = Pos(1)*0.78;
% % hLeg.Position(2) = Pos(2)+Pos(2)*0.1;
% 
% icons(1).FontSize = fSet.axLabFontsize;
% icons(2).FontSize = fSet.axLabFontsize;
% icons(3).FontSize = fSet.axLabFontsize;
% icons(4).FontSize = fSet.axLabFontsize;
% icons(5).FontSize = fSet.axLabFontsize;
% % POS = icons(1).Position;
% % icons(1).Position = POS + [0.1 0 0];
% % POS = icons(2).Position;
% % icons(2).Position = POS + [0.1 0 0];
% % POS = icons(3).Position;
% % icons(3).Position = POS + [0.1 0 0];
% % POS = icons(4).Position;
% % icons(4).Position = POS + [0.1 0 0];
% % POS = icons(5).Position;
% % icons(5).Position = POS + [0.1 0 0];
% 
% XDAT = icons(6).XData;
% icons(6).LineStyle = '-';
% icons(6).LineWidth = fSet.plLineWidth;
% % icons(6).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];
% icons(8).LineStyle = '-';
% icons(8).LineWidth = fSet.plLineWidth;
% % icons(8).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];
% icons(10).LineStyle = '-';
% icons(10).LineWidth = fSet.plLineWidth;
% % icons(10).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];
% icons(12).LineStyle = '-';
% icons(12).LineWidth = fSet.plLineWidth;
% % icons(12).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];
% icons(14).LineStyle = '-';
% icons(14).LineWidth = fSet.plLineWidth;
% icons(14).XData = [-XDAT(1) XDAT(2)+XDAT(1)*2];

if ~alignN2c
    xlim([-20 350])
else
    xlim([-50 50])
end
ylim([-3 0.5])
% ylim([-20 3])%CSD

xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
% ylabel({'N2c \muV'},'fontsize',fSet.plFontsize)
ylabel({'Amplitude (\muV)'},'fontsize',fSet.axLabFontsize);

xBox = [min(tmpN2c_lat)-min(tmpN2c_lat)*0.05 max(tmpN2c_lat)+min(tmpN2c_lat)*0.05];
yBox = [min(y_latency)+min(y_latency)*0.22 max(y_latency)-min(y_latency)*0.22];

text(max(xBox)*1.02, max(yBox)*0.85, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([-200 xBox(1)], [0 0] ,' k','linewidth',1)
plot([xBox(2) 800], [0 0] ,' k','linewidth',1)

yBox = [min(tmpN2c_amp)+min(tmpN2c_amp)*0.12 max(tmpN2c_amp)-min(tmpN2c_amp)*0.12];
xBox = [min(x_amplitude)-min(x_amplitude)*0.5 max(x_amplitude)+min(x_amplitude)*0.5];

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

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
if ~isempty(mfit)
    % plot model fit
    h = plot(mfit.mean, x_amplitude, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], [y_latency flip(y_latency)] , h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
text(min(xBox)*1.02, max(yBox)*0.85, {pstring_chi}, 'fontsize',fSet.plFontsize)




% text(xBox(2)*0.6, yBox(2)*0.94, 'n.s.','fontsize',25)

% inset with topoplot
% axes('position',[0.15 0.5 0.4 0.35]) ; % inset
% maplimits = [min(min(squeeze(mean(mean(N2c_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(N2c_topo{idx2plot},1),2))))];
% topoplot(squeeze(mean(mean(N2c_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub3 = subtightplot(nrow, ncol, [5], subplotgap2);
figInit('ax');

% hsub3 = subplot(nrow, ncol, [5]);
colormap(hsub3,'Default')

t2test      = [200 400];
f2test      = [0.1 4]; % in reality 2 4, but this translates to this in the plot

ttSPG = SPG_times-abs(t(1)/1000);
ttSPG = ttSPG * 1000;
maplimits = [min(min(min(squeeze(mean(N2c_ITPC{idx2plot},1))))) max(max(max(squeeze(mean(N2c_ITPC{idx2plot},1)))))];
maplimits = [0 max(max(max(squeeze(mean(N2c_ITPC{idx2plot},1)))))];

contourf(ttSPG, SPG_freq, squeeze(mean(mean(N2c_ITPC{idx2plot},1),2))',20,'LineColor','none');
set(gca,'clim',maplimits)
ylim([0 18])
axis xy
xlim([-50 450])
% ylim([-0.1 34])
hold on

h = rectangle('Position',[t2test(1) f2test(1) diff(t2test) diff(f2test)]);%, '-w','linewidth',3)
h.LineWidth = fSet.plLineWidth;
h.EdgeColor = [1 1 1];

ylabel('Frequency (Hz)','fontsize',fSet.axLabFontsize)
xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize)
figInit('ax');
set(gca,'color','white')
B=colorbar;
% set(B, 'Position', [.345 .375 .01981 .1], 'Limits', maplimits)
xpos = hsub3.Position(3) + hsub3.Position(1) - .01981 - 0.01;
ypos = hsub3.Position(4) + hsub3.Position(2) - .1 - 0.01;
set(B, 'Position', [xpos ypos .01981 .1], 'Limits', maplimits)
B.AxisLocation = 'in';
B.FontWeight = 'bold';
% title(B,'ITPC','fontsize',fSet.axFontsize)
% axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub4 = subtightplot(nrow, ncol, [6], subplotgap);
% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');
set(gca,'color','white')

x_ITPC = [-80:50:120]+5;

ttSPG = SPG_times-abs(t(1)/1000);
ttSPG = ttSPG*1000;
clear LEGEND
hold on

tmpN2c = nanmean(N2c_ITPC_bar{idx2plot});
for ibin = 1:nbin2use
        plot([-200 700],[tmpN2c(ibin) tmpN2c(ibin)],'linewidth',2,'color',[fSet.colors(ibin,:) 0.4]);
end
plot([t2test(1) t2test(2)], [0.11 0.11], 'k', 'linewidth', 4)

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
%     plot model fit
    h = plot(x_ITPC, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([x_ITPC flip(x_ITPC)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end
yBox = [min(tmpN2c)-min(tmpN2c)*0.07 max(tmpN2c)+min(tmpN2c)*0.07];
xBox = [min(x_ITPC)-35 max(x_ITPC)+35];
text(10, max(yBox)*1.04, {pstring_chi}, 'fontsize',fSet.plFontsize)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';


clear h h2 LEGEND
for ibin = 1:nbin2use
    h(ibin) = boundedline(ttSPG, squeeze(mean(N2c_ITPC_band{idx2plot}(:,ibin, :, :),1)), squeeze(std(N2c_ITPC_band{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % se
    tmpCPPstd = squeeze(std(N2c_ITPC_band{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpN2c(ibin)-tmpCPPstd tmpN2c(ibin)+tmpCPPstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    % plotMarker
    h2(ibin) = plot(x_ITPC(ibin), tmpN2c(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end


% plot([0 0], [0 0.5] ,'k','linewidth',1)

plot([0 0], [0 yBox(1)] ,' k','linewidth',1)
plot([0 0], [yBox(2) 0.5] ,' k','linewidth',1)

xlim([-130 400])
ylim([0.1 0.33])
% axis square
xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize)
ylabel('ITPC','fontsize',fSet.axLabFontsize)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub4 = subtightplot(nrow, ncol, [7 ], subplotgap3 .* [1 0.2]);
% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');

cfg.xlim.beta = [-500 80];
cfg.ylim.beta = [0.5 0.76];
cfg.ylim.beta = [0.57 0.76];


x_beta = [-440:50:-240];

spectral_t = stft_timesr;

hold on

tmpBeta = mean(beta_response_mean{idx2plot});
tmpBeta_slope = mean(beta_pre_response_slope{idx2plot});
for ibin = 1:nbin2use
    plot([-600 100],[tmpBeta(ibin) tmpBeta(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
%
% plot(x_beta, tmpBeta,'linewidth',4,'color',[0.5 0.5 0.5 0.5])
plot([0 0], [0 1] ,'k','linewidth',1)

twin_bar    = [-100 0];
t2test_slope= [-300 0];

% plot([twin_bar(1) twin_bar(2)], [-0.13 -0.13], 'k', 'linewidth', 6)
% plot([twin_bar(1) twin_bar(2)], [0.515 0.515], 'k', 'linewidth', 4)
% plot([t2test_slope(1) t2test_slope(2)], [0.505 0.505], 'color',[0.5 0.5 0.5], 'linewidth', 4)
plot([twin_bar(1) twin_bar(2)], cfg.ylim.beta(1) + [0.012 0.012], 'k', 'linewidth', 4)
plot([t2test_slope(1) t2test_slope(2)], cfg.ylim.beta(1) + [0.005 0.005], 'color',[0.5 0.5 0.5], 'linewidth', 4)

yBox = [min(tmpBeta)-min(tmpBeta)*0.02 max(tmpBeta)+min(tmpBeta)*0.02];
xBox = [min(x_beta)+min(x_beta)*0.075 max(x_beta)-min(x_beta)*0.075];

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'preRespBeta',nSub,nbin2use);
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
    h(ibin) = boundedline(spectral_t, squeeze(mean(beta_response{idx2plot}(:,ibin, :, :),1)), squeeze(std(beta_response{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % amplitude se
    tmpstd = squeeze(std(beta_response_mean{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_beta(ibin) x_beta(ibin)],[tmpBeta(ibin)-tmpstd tmpBeta(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot(x_beta(ibin), tmpBeta(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');    
    
end


xlim([cfg.xlim.beta])
ylim([cfg.ylim.beta])
set(gca,'fontsize',fSet.plFontsize)
xlabel('Time from response (ms)','fontsize',fSet.axLabFontsize)
ylabel({'Power (dB)'},'fontsize',fSet.axLabFontsize)


text(xBox(1)*0.99, yBox(1)*0.985, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% beta slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub4 = subtightplot(nrow, ncol, [8 ], subplotgap2);
figInit('ax');

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
ylim([min(tmpBeta_slope)*1.15 max(tmpBeta_slope)*0.85])

YTICK = ax2.YTick;
ax2.YTick = linspace(min(YTICK), max(YTICK), 3);
% text(5, max(tmpBeta_slope)*0.85, '*','fontsize',35)

text(3, min(mean(beta_pre_response_slope{idx2plot})) * 1.1, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('LHB Slope', 'fontsize',fSet.axLabFontsize)




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
%% control, fig 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fH, fSet] = figInit('fig', 3, {'height', 3.5/5; 'width',6/5});

nrow = 3;
ncol = 5;
subplotgap = [0.1 0.065];
subplotgap2 = [0.1 0.065];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% alpha_asym
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hsub = subtightplot(nrow, ncol, [1], subplotgap);

% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');

% cfg.ylim.alpha = [2.1 3];
hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'alpha_asym',nSub,nbin2use);
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
    h = plot(1:5, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([1:5 flip(1:5)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
end

for ibin = 1:nbin2use
    plot([x(ibin) x(ibin)],...
        [mean(alpha_asym{idx2plot}(:,ibin))-std(alpha_asym{idx2plot}(:,ibin))/sqrt(nSub) mean(alpha_asym{idx2plot}(:,ibin))+std(alpha_asym{idx2plot}(:,ibin))/sqrt(nSub) ],...
        'linewidth',fSet.plLineWidth_in,'color','k');
    
    plot(x(ibin), mean(alpha_asym{idx2plot}(:,ibin)),plotMarker{ibin}, ...
        'color', fSet.colors(ibin,:),'markersize', fSet.MarkerSize,...
        'MarkerFaceColor', fSet.colors(ibin,:),'MarkerEdgeColor','k')
end

text(3.3, min(mean(alpha_asym{idx2plot}))*0.65, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)


% ylim([cfg.ylim.alpha])
h = ylabel('\alpha Power Asymmetry','fontsize',fSet.axLabFontsize);
% xlim([0.5 5.5])
xlim([0 nbin2use+1])

set(gca,'xtick',x)
xlabel('Pupil Bin','fontsize',fSet.axLabFontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% alpha asym, pupil response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [2], subplotgap);
figInit('ax');
hold on

tmpPupil = mean(pupil_lp_RT_neg200_200{idx_alpha_asym});

% get stats
bin2use = allBin2use{idx_alpha_asym};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'plevelPupil_RT200',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
text(3, min(tmpPupil)*0.8, {pstring_chi}, 'fontsize',fSet.plFontsize)

for ibin = 1:nbin2use
    % se
    tmpPupil_std = squeeze(std(pupil_lp_RT_neg200_200{idx_alpha_asym}(:,ibin)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpPupil(ibin)-tmpPupil_std tmpPupil(ibin)+tmpPupil_std],'linewidth',fSet.plLineWidth_in,'color','k');
    % plotMarker
    h2(ibin) = plot((ibin), tmpPupil(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize, 'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k','LineWidth',1);
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpPupil)*0.75 max(tmpPupil)*1.15])

set(gca,'xtick',1:nbin2use)
xlabel('\alpha asymmetry bin','fontsize',fSet.axLabFontsize)
ylabel('Pupil response','fontsize',fSet.axLabFontsize)

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2i, timecourse and latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if 1
N2i2plot = squeeze(N2i{idx2plot});
N2i_latency2plot = squeeze(N2i_latency{idx2plot});
N2i_amplitude2plot = squeeze(N2i_amplitude{idx2plot});


% x and y values where to plot N2 latency and amplitude
% y_latency = [0.5:-0.17:-0.18];
y_latency = [0.7:-0.2:-0.1]-0.1;
x_amplitude = [50:27:160];

% times N2 values have been computed on
t2test_amp      = 340 + [-50 50];        
t2test_lat      = [150 400];

hsub1 = subtightplot(nrow, ncol, [3 4], subplotgap);
figInit('ax');

hold on

plot([0 0], [-4 1] ,' k','linewidth',1)

plot(t2test_amp, [-1.55 -1.55] ,' k','linewidth',4)

% plot lines to indicate N2c lat and amp
tmpN2i_lat = mean(N2i_latency2plot);
for ibin = 1:nbin2use
    plot([tmpN2i_lat(ibin) tmpN2i_lat(ibin)], [-20 1],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
tmpN2i_amp = mean(N2i_amplitude2plot);
for ibin = 1:nbin2use
    plot([0 400],[tmpN2i_amp(ibin) tmpN2i_amp(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
% get stats
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

xBox = [min(tmpN2i_lat)-min(tmpN2i_lat)*0.05 max(tmpN2i_lat)+min(tmpN2i_lat)*0.05];
yBox = [min(y_latency)-max(y_latency)*0.22 max(y_latency)+max(y_latency)*0.22];
text(max(xBox)*1.02, max(yBox)*0.85, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

plot([-200 xBox(1)], [0 0] ,' k','linewidth',1)
plot([xBox(2) 800], [0 0] ,' k','linewidth',1)

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'N2i_amplitude',nSub,nbin2use);
[pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');

yBox = [min(tmpN2i_amp)+min(tmpN2i_amp)*0.22 max(tmpN2i_amp)-min(tmpN2i_amp)*0.22];
xBox = [min(x_amplitude)-min(x_amplitude)*0.6 max(x_amplitude)+min(x_amplitude)*0.6];

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';


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
text(min(xBox)*1.02, min(yBox)*1.15, [betaString ' ' pstring_chi], 'fontsize',fSet.plFontsize)

if ~isempty(mfit)    
    % plot model fit
    h = plot(x_amplitude, mfit.mean, 'linewidth', fSet.plLineWidth_in, 'color', model_color);
    h = patch([x_amplitude flip(x_amplitude)], [mfit.mean flip(mfit.mean)] + [-mfit.se flip(mfit.se)], h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];

end

clear h h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2c
    h(ibin) = boundedline(cfg.t, squeeze(mean(N2i2plot(:,ibin, :),1)), squeeze(std(N2i2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');

    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % se
    tmpN2i_lat_std = std(N2i_latency2plot(:,ibin))/sqrt(nSub);
    plot([tmpN2i_lat(ibin)-tmpN2i_lat_std tmpN2i_lat(ibin)+tmpN2i_lat_std], [y_latency(ibin) y_latency(ibin)],'linewidth',fSet.plLineWidth_in,'Color','k')
    h2(ibin) = plot(tmpN2i_lat(ibin), y_latency(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    %     plot([tmpCPP(ibin)-tmpCPPstd tmpCPP(ibin)+tmpCPPstd], [y_onset(ibin) y_onset(ibin)],'linewidth',4,'color',[fSet.colors(ibin,:)]);
    
    tmpN2i_amp_std = std(N2i_amplitude2plot(:,ibin))/sqrt(nSub);
    plot([x_amplitude(ibin) x_amplitude(ibin)], [tmpN2i_amp(ibin)-tmpN2i_amp_std tmpN2i_amp(ibin)+tmpN2i_amp_std], 'linewidth',fSet.plLineWidth_in,'Color','k')
    h2(ibin) = plot(x_amplitude(ibin),tmpN2i_amp(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    
    
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
    
end

xlim([-20 500])
ylim([-1.7 0.75])

xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
% ylabel({'N2c \muV'},'fontsize',fSet.plFontsize)
ylabel({'Amplitude (\muV)'},'fontsize',fSet.axLabFontsize);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2i, pupil response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [5], subplotgap);
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







saveFigName = [bin2use '_' bintype fileExt '_control'];
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

%

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
% axis square



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

%% fig 5
%%%%% development of RT across trials

[fH, fSet] = figInit('fig', 3, {'height', 3.5/5; 'width',6/5});

nrow = 3;
ncol = 4;
subplotgap = [0.1 0.058];


allRTnames = {'RT_min5','RT_min4','RT_min3','RT_min2','RT_min1','RT_0','RT_plus1','RT_plus2','RT_plus3','RT_plus4','RT_plus5'};
x2plot = -5:5;

hsub = subtightplot(nrow, ncol, [9 10 ], subplotgap);
figInit('ax');


cfg.xlim     = [-5 5];
cfg.ylim     = [-0.2 0.2];

hold on

for ibin = 1:nbin2use
    A = boundedline(x2plot,squeeze(mean(RT_window{idx2plot}(:,ibin,:,:),1)),squeeze(std(RT_window{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(A,'linewidth',fSet.plLineWidth_in)
end
% plot([0 0], cfg.ylim ,'k','linewidth',1)
% plot(x2plot(signRTbin), cfg.ylim(1)*0.9, 'k*','markersize',10)
% sigstar(num2cell(find(signRTbin)), pRTbin(signRTbin))
xlim(cfg.xlim)
ylim(cfg.ylim)
set(gca,'xtick',-5:5)
xlabel('Trial relative to pupil response','fontsize',fSet.axLabFontsize);
ylabel('RT (Z-score)','fontsize',fSet.axLabFontsize);



RT2compare = [];
clear P
for itrial = 1:11
    RT2compare = [mean(RT_window{idx2plot}(:,[1 2 3 4],:,itrial),2), RT_window{idx2plot}(:,[5],:,itrial)];
    
    [H,P(itrial,1),CI,STATS] = ttest(RT2compare(:,1), RT2compare(:,2));
end

signRTbin = false(11,1);
pt = FDR(P,0.05);
% pt = 0.05/numel(P) ;
signRTbin(P<pt) = true;
P(~signRTbin) = 1;
[pstring_chi,starstring] = getSignificanceStrings(P, 0, 1);

YLIM = get(gca,'ylim');
for iP = 1:length(P)
    text(x2plot(iP), YLIM(2)*0.95, starstring{iP},'fontsize',fSet.axLabFontsize,'HorizontalAlignment','center')
end 

% axes('position',[0.55 0.6 0.4 0.4]) ; % inset
hsub = subtightplot(nrow, ncol, [11 12 ], subplotgap);

figInit('ax');
% set(gca,'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4])
hold on

% bin2use = allBin2use{idx2plot};
% filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
% [stats, mfit] = loadStatsR(filename, 'RT_0',nSub,nbin2use);
% [pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, '\it p ');
% betaString = ['\beta_{2} = ' num2str(stats.B, '%1.2f')];

RT2compare = [];
x2plot = 1:3:nbin2use*3;
clear P
for ibin = 1:nbin2use
    RT2compare = [RTmin1{idx2plot}(:,ibin), RT0{idx2plot}(:,ibin), RTplus1{idx2plot}(:,ibin)];
    
    [H,P(ibin,1),CI,STATS] = ttest(RT2compare(:,1), RT2compare(:,2));
    [H,P(ibin,2),CI,STATS] = ttest(RT2compare(:,2), RT2compare(:,3));
%     [H,P(ibin,3),CI,STATS] = ttest(RT2compare(:,3), RT2compare(:,4));
       
    
    h = bar([x2plot(ibin):x2plot(ibin)+2], mean(RT2compare));
    h.FaceColor = fSet.colors(ibin,:);
    
    plot([x2plot(ibin) x2plot(ibin)], mean(RT2compare(:,1)) + [-std(RT2compare(:,1))/sqrt(nSub) std(RT2compare(:,1))/sqrt(nSub)],'linewidth',1,'color','k')
    plot([x2plot(ibin)+1 x2plot(ibin)+1], mean(RT2compare(:,2)) + [-std(RT2compare(:,2))/sqrt(nSub) std(RT2compare(:,2))/sqrt(nSub)],'linewidth',1,'color','k')
    plot([x2plot(ibin)+2 x2plot(ibin)+2], mean(RT2compare(:,3)) + [-std(RT2compare(:,3))/sqrt(nSub) std(RT2compare(:,3))/sqrt(nSub)],'linewidth',1,'color','k')
%     plot([x2plot(ibin)+3 x2plot(ibin)+3], mean(RT2compare(:,4)) + [-std(RT2compare(:,4))/sqrt(nSub) std(RT2compare(:,4))/sqrt(nSub)],'linewidth',1,'color','k')
     
end

signRTbin = false(nbin2use,2);
% pt = FDR(P,0.05);
pt = 0.05/numel(P) ;
signRTbin(P<pt) = true;
% signRTbin = find(signRTbin);
clear groups
for ibin = 1:nbin2use
    
    for icomp = 1:2
        if signRTbin(ibin,icomp)
            groups{ibin,icomp} = [x2plot(ibin) x2plot(ibin)+1] + (icomp-1);
        end
    end

end
groups = groups(signRTbin);
sigstar(groups(:), P(signRTbin))

% P = P(:);


set(gca,'xtick',1:15)
xTicks = get(gca, 'xtick');



YLIM = ylim;
text(14,YLIM(1)*0.6, {'Trial index'},'horizontalalignment','center')
text(13,YLIM(1)*0.8, {'-1'},'horizontalalignment','center')
text(14,YLIM(1)*0.8, {'0'},'horizontalalignment','center')
text(15,YLIM(1)*0.8, {'1'},'horizontalalignment','center')
A = get(gca);
plot([13 13],[YLIM(1) YLIM(1)+A.TickLength(1)],'k','linewidth',2)
plot([14 14],[YLIM(1) YLIM(1)+A.TickLength(1)],'k','linewidth',2)
plot([15 15],[YLIM(1) YLIM(1)+A.TickLength(1)],'k','linewidth',2)
% get(gca,'
% xtickloc = 

set(gca,'xtick',[2:3:nbin2use*3],'xticklabel',1:nbin2use)
xlabel('Pupil bin','fontsize',fSet.axLabFontsize);
% xlim([0 16])

% axis square


saveFigName = [bin2use '_' bintype fileExt '_RT_acrossTrials'];

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
%% reduced - CPP
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
    tmpCPPstd = squeeze(std(CPP_ITPC_band{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_ITPC(ibin) x_ITPC(ibin)],[tmpCPP(ibin)-tmpCPPstd tmpCPP(ibin)+tmpCPPstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
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

[fH, fSet] = figInit('fig', 7, {'height', 2/5; 'width',6/5});

nrow = 2;
ncol = 9;
subplotgap = [0.12 0.058];
subplotgap2 = [0.12 0.07];

subplotorder = [...
    1 2;...
    3 4;...
    5 7;...
    8 9;...
    10 12;...
    13 14;...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% alpha - asym
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [subplotorder(2,1):subplotorder(2,2)], subplotgap2);
figInit('ax');

hold on

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'alpha_asym',nSub,nbin2use);
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
tmpmean = squeeze(mean(alpha_asym{idx2plot}));
for ibin = 1:nbin2use
    %     slope
    tmpstd = squeeze(std(alpha_asym{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [(ibin) (ibin)],[tmpmean(ibin)-tmpstd tmpmean(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot((ibin), tmpmean(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    LEGEND{ibin} = ['Bin ' num2str(ibin) ];
end
xlim([0 6])
ylim([min(tmpmean) max(tmpmean)] .* [0.1 1.4])

text(1.5, min(tmpmean) * 0.4, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('\alpha power asymmetry', 'fontsize',fSet.axLabFontsize)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N2, timecourse and latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [subplotorder(3,1):subplotorder(3,2)], subplotgap2);
figInit('ax');
hold on

N2c2plot = squeeze(N2c{idx2plot});
N2c_latency2plot = squeeze(N2c_latency{idx2plot});
N2c_amplitude2plot = squeeze(N2c_amplitude{idx2plot});

% x and y values where to plot N2 latency and amplitude
y_latency = [0.3:-0.3:-0.9];
x_amplitude = [40:30:160];

% times N2 values have been computed on
t2test_amp      = 266 + [-50 50];
t2test_lat      = [150 400];

plot([0 0], [-4 1] ,' k','linewidth',1)
plot(t2test_amp, [-3 -3] ,' k','linewidth',4)

% plot lines to indicate N2c lat and amp
tmpN2c_lat = mean(N2c_latency2plot);
for ibin = 1:nbin2use
    plot([tmpN2c_lat(ibin) tmpN2c_lat(ibin)], [-20 1],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
tmpN2c_amp = mean(N2c_amplitude2plot);
for ibin = 1:nbin2use
    plot([0 400],[tmpN2c_amp(ibin) tmpN2c_amp(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
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
xBox = [min(tmpN2c_lat)-min(tmpN2c_lat)*0.05 max(tmpN2c_lat)+min(tmpN2c_lat)*0.05];
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
yBox = [min(tmpN2c_amp)+min(tmpN2c_amp)*0.12 max(tmpN2c_amp)-min(tmpN2c_amp)*0.12];
xBox = [min(x_amplitude)-min(x_amplitude)*0.5 max(x_amplitude)+min(x_amplitude)*0.5];

clear h
h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

text(min(xBox)*1.02, max(yBox)*0.7, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

clear h1 h2 LEGEND
for ibin = 1:nbin2use
    %     plot N2c
    h1 = boundedline(cfg.t, squeeze(mean(N2c2plot(:,ibin, :),1)), squeeze(std(N2c2plot(:,ibin, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'transparency', 0.2,'alpha');
    set(h1,'linewidth',fSet.plLineWidth)
    
    % se
    tmpN2c_lat_std = std(N2c_latency2plot(:,ibin))/sqrt(nSub);
    plot([tmpN2c_lat(ibin)-tmpN2c_lat_std tmpN2c_lat(ibin)+tmpN2c_lat_std], [y_latency(ibin) y_latency(ibin)],'linewidth',fSet.plLineWidth_in,'Color','k')
    tmpN2c_amp_std = std(N2c_amplitude2plot(:,ibin))/sqrt(nSub);
    plot([x_amplitude(ibin) x_amplitude(ibin)], [tmpN2c_amp(ibin)-tmpN2c_amp_std tmpN2c_amp(ibin)+tmpN2c_amp_std], 'linewidth',fSet.plLineWidth_in,'Color','k')
    
    
    plot(tmpN2c_lat(ibin), y_latency(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
    plot(x_amplitude(ibin),tmpN2c_amp(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');
      
end

xlim([-20 350])
ylim([-3.2 0.5])
% ylim([-20 3])%CSD

xlabel('Time from motion onset (ms)','fontsize',fSet.axLabFontsize);
ylabel({'Amplitude (\muV)'},'fontsize',fSet.axLabFontsize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ITPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub = subtightplot(nrow, ncol, [subplotorder(4,1):subplotorder(4,2)], subplotgap2);
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

text(1.5, min(tmpmean) * 0.85, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('N2c ITPC', 'fontsize',fSet.axLabFontsize)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub4 = subtightplot(nrow, ncol, [subplotorder(5,1):subplotorder(5,2)], subplotgap2);
% hsub4 = subplot(nrow, ncol, [6]);
figInit('ax');

cfg.xlim = [-500 80];
cfg.ylim = [0.57 0.76];


x_beta = [-440:40:-280];

spectral_t = stft_timesr;

hold on

tmpBeta = mean(beta_response_mean{idx2plot});
tmpBeta_slope = mean(beta_pre_response_slope{idx2plot});
for ibin = 1:nbin2use
    plot([-600 100],[tmpBeta(ibin) tmpBeta(ibin)],'linewidth',fSet.plLineWidth_in,'color',[fSet.colors(ibin,:) 0.4]);
end
%
% plot(x_beta, tmpBeta,'linewidth',4,'color',[0.5 0.5 0.5 0.5])
plot([0 0], [0 1] ,'k','linewidth',1)

twin_bar    = [-100 0];
t2test_slope= [-300 0];

% plot([twin_bar(1) twin_bar(2)], [-0.13 -0.13], 'k', 'linewidth', 6)
% plot([twin_bar(1) twin_bar(2)], [0.515 0.515], 'k', 'linewidth', 4)
% plot([t2test_slope(1) t2test_slope(2)], [0.505 0.505], 'color',[0.5 0.5 0.5], 'linewidth', 4)
plot([twin_bar(1) twin_bar(2)], cfg.ylim(1) + [0.012 0.012], 'k', 'linewidth', 4)
plot([t2test_slope(1) t2test_slope(2)], cfg.ylim(1) + [0.005 0.005], 'color',[0.5 0.5 0.5], 'linewidth', 4)

yBox = [min(tmpBeta)-min(tmpBeta)*0.02 max(tmpBeta)+min(tmpBeta)*0.02];
xBox = [min(x_beta)+min(x_beta)*0.075 max(x_beta)-min(x_beta)*0.075];

h = rectangle('Position',[xBox(1) yBox(1) abs(diff(xBox)) diff(yBox)]);
h.LineWidth = 1;
h.LineStyle = '--';

% get stats
bin2use = allBin2use{idx2plot};
filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
[stats, mfit] = loadStatsR(filename, 'preRespBeta',nSub,nbin2use);
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
    h(ibin) = boundedline(spectral_t, squeeze(mean(beta_response{idx2plot}(:,ibin, :, :),1)), squeeze(std(beta_response{idx2plot}(:,ibin, :, :),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
    set(h(ibin),'linewidth',fSet.plLineWidth)
    
    % amplitude se
    tmpstd = squeeze(std(beta_response_mean{idx2plot}(:,ibin, :, :)))/sqrt(nSub);
    plot( [x_beta(ibin) x_beta(ibin)],[tmpBeta(ibin)-tmpstd tmpBeta(ibin)+tmpstd],'linewidth',fSet.plLineWidth_in,'color','k');
    
    h2(ibin) = plot(x_beta(ibin), tmpBeta(ibin), plotMarker{ibin},'markersize',fSet.MarkerSize,'color',fSet.colors(ibin,:),'MarkerFaceColor',fSet.colors(ibin,:),'MarkerEdgeColor','k');    
    
end


xlim([cfg.xlim])
ylim([cfg.ylim])
set(gca,'fontsize',fSet.plFontsize)
xlabel('Time from response (ms)','fontsize',fSet.axLabFontsize)
ylabel({'Power (dB)'},'fontsize',fSet.axLabFontsize)


text(xBox(1)*0.99, yBox(1)*0.985, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% beta slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsub4 = subtightplot(nrow, ncol, [subplotorder(6,1):subplotorder(6,2)], subplotgap2);
figInit('ax');

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
ylim([min(tmpBeta_slope)*1.2 max(tmpBeta_slope)*0.85])

text(3, min(mean(beta_pre_response_slope{idx2plot})) * 1.15, {betaString, pstring_chi}, 'fontsize',fSet.plFontsize)

set(gca,'xtick',1:nbin2use)
xlabel('Pupil bin', 'fontsize',fSet.axLabFontsize)
ylabel('LHB Slope', 'fontsize',fSet.axLabFontsize)




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

[fH, fSet] = figInit('fig', 6, {'height', 2/5; 'width',6/5});

nrow = 2;
ncol = 9;
subplotgap = [0.12 0.058];
subplotgap2 = [0.12 0.07];

% hsub4 = subtightplot(nrow, ncol, [9 10 ], subplotgap2);

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
% axis square



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subtightplot(nrow, ncol, [subplotorder(8,1):subplotorder(8,2)], subplotgap2)
figInit('ax');
hold on
RT = RT;

% idx2plot    = [find(idx_BL_bp) find(idx2plot)];
% idx2plot    = [find(idx2plot)];
% % idx2plot    = [find(idx_BL_lp) find(idx_BL_bp)];
% idx2plot = idx2plot2;

cfg.ylim.RT     = [1000 1200];
cfg.ylim.RT_CV  = [0.175 0.32];
cfg.ylabel.RT   = 'RT (ms)';

% [stats2report, meanFit, CIFit, SEFit, STATS] = fitlme_pupil_singleVar(MAT, 'RT', 0);
% [pstring_chi,starstring] = getSignificanceStrings(stats2report(3), 0, 1, '\it p ');
% betaString = ['\beta_{2} = ' num2str(STATS.Estimate, '%1.3f')];

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


%alpha asym
figure(6), clf
set(gcf,'color','none')
set(gca,'color','none')
set(gca,'linewidth',2)
set(gcf,'position',[-1420 154 1092 981])
set(gcf,'PaperPositionMode','auto')

maplimits = [min(min(squeeze(mean(mean(alpha_asym_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(alpha_asym_topo{idx2plot},1),2))))];
topoplot(squeeze(mean(mean(alpha_asym_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);
saveFigName = [bin2use '_' bintype fileExt '_Alpha_asym_topo'];
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

maplimits = [min(min(squeeze(mean(mean(beta_pre_response_topo{idx2plot},1),2)))) max(max(squeeze(mean(mean(beta_pre_response_topo{idx2plot},1),2))))];
topoplot(squeeze(mean(mean(beta_pre_response_topo{idx2plot},1),2)),chanlocs,'maplimits',maplimits,'electrodes','off','plotchans',plot_chans,'numcontour',0);

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

