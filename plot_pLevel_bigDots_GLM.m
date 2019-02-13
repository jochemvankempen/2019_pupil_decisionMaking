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
% plot_pLevel_bigDots_GLM

paths.fig = [paths.pop 'fig' filesep 'manuscript' filesep 'GLM' filesep];
if ~exist(paths.fig,'dir')
    mkdir(paths.fig)
end

figureFileType = {'png','svg'};
clear cfg
cfg.t = t;
cfg.tr = tr;
cfg.t_pupil = t_pupil;
cfg.tr_pupil = tr_pupil;
plotMarker = {'o','s','d','^','h'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary figure 1. plot the average BIC scores, determine which model fits best. Also plot effect size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, selectedModel] = min(mean(allGLM_conc_pupil_BIC));

[figHandle, fSet] = figInit('fig', 1, {...
    'figureType', 'Manuscript'; ...
    'width', 22; ...
    'height', 10});

subplotgap = [0.08 0.08];
subplotmargin = [0.1 0.1];

axH = subtightplot(1,3,[1 2],subplotgap, subplotmargin, subplotmargin);
figInit('ax');
hold on
plot_subplot_label(axH, fSet, 'A', [0.08 0.05])

hBar = bar(mean(allGLM_conc_pupil_BIC) - mean(allGLM_conc_pupil_BIC(:,selectedModel)));
% plot( 1:size(allGLM_conc_pupil_BIC,2), std(allGLM_conc_pupil_BIC(:,imodel) - allGLM_conc_pupil_BIC(:,selectedModel)))

% plot se
for imodel = 1:size(allGLM_conc_pupil_BIC,2)
    plot([imodel imodel],...
        [(mean(allGLM_conc_pupil_BIC(:,imodel)) - mean(allGLM_conc_pupil_BIC(:,selectedModel))) - (std(allGLM_conc_pupil_BIC(:,imodel) - allGLM_conc_pupil_BIC(:,selectedModel))/sqrt(nSub)) ...
        (mean(allGLM_conc_pupil_BIC(:,imodel)) - mean(allGLM_conc_pupil_BIC(:,selectedModel))) + (std(allGLM_conc_pupil_BIC(:,imodel) - allGLM_conc_pupil_BIC(:,selectedModel))/sqrt(nSub)) ],...
        'color', 'k','linewidth',fSet.LineWidth);
    
    mod2compare = allGLM_conc_pupil_BIC(:, imodel) - allGLM_conc_pupil_BIC(:,selectedModel);
    [P, H, STATS] = signrank(mod2compare, 0, 'method', 'exact');
    [pstring_chi,starstring] = getSignificanceStrings(P, 0, 1, 'p ');
    text(imodel, max(hBar.YData) * 0.01, [starstring], 'fontsize', fSet.Fontsize_title + 3, 'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', 'Ariel')
end

hBar.FaceColor = [0.5 0.5 0.5];
set(gca,'xtick', 1:size(allGLM_conc_pupil_BIC,2))
xlabel('Model number', 'fontsize', fSet.Fontsize_text)
ylabel('\Delta BIC', 'fontsize', fSet.Fontsize_text)

axH = subtightplot(1,3,[3],[0.08 0.08], [0.1 0.1], [0.1 0.1]);
figInit('ax');
hold on
plot_subplot_label(axH, fSet, 'B')

switch selectedModel
    case 1
        tstat2plot = allGLM_conc_pupil_StimResp_tstat;        
        modelName = ['Model ' num2str(selectedModel) ': Onset-Response'];        
    case 2
        tstat2plot = allGLM_conc_pupil_Boxc_tstat;
        modelName = ['Model ' num2str(selectedModel) ': Boxcar'];        
    case 3
        tstat2plot = allGLM_conc_pupil_Ramp_tstat;
        modelName = ['Model ' num2str(selectedModel) ': Ramp'];        
    case 4
        tstat2plot = allGLM_conc_pupil_Ramp2thresh_tstat;
        modelName = ['Model ' num2str(selectedModel) ': Ramp to threshold'];        
end

hold on
hBar = bar(mean(tstat2plot) );
% plot se
for icomponent = 1:size(tstat2plot,2)
    plot([icomponent icomponent],...
        [mean(tstat2plot(:,icomponent)) - (std(tstat2plot(:,icomponent))/sqrt(nSub)) ...
        mean(tstat2plot(:,icomponent)) + (std(tstat2plot(:,icomponent))/sqrt(nSub))],...
        'color', 'k','linewidth',fSet.LineWidth);
    
    [P, H, STATS] = signrank(tstat2plot(:, icomponent), 0, 'method', 'exact');
    [pstring_chi,starstring] = getSignificanceStrings(P, 0, 1, '\it p ');
    text(icomponent, max(hBar.YData) * 0.01, [starstring], 'fontsize', fSet.Fontsize_title + 3, 'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', 'Ariel')

end

hBar.FaceColor = [0.5 0.5 0.5];
set(gca,'xtick', 1:size(tstat2plot,2))
set(gca,'xticklabel',{'Onset', 'Sustained', 'Response'}, 'xticklabelrotation', 25, 'fontsize', fSet.Fontsize_text)
ylabel('Effect size (t)', 'fontsize', fSet.Fontsize_text)
title(modelName, 'fontsize', fSet.Fontsize_title)

saveFigName = ['GLM_conc_modelComparison' fileExt_GLM];

figSave(saveFigName, paths.fig, {'png', 'svg'})




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the average BIC scores per bin, determine which model fits best. Also plot effect size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelIdx = 1:3;

[figHandle, fSet] = figInit('fig', 1, {...
    'figureType', 'Manuscript'; ...
    'width', 22; ...
    'height', 22});

subplotgap = [0.08 0.08];
subplotmargin = [0.1 0.1];

subplotIdx = [1:2:10];
nbin = 5;

for ibin = 1:nbin
    
    %%% BIC scores
    axH = subtightplot(nbin,2,subplotIdx(ibin),subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
    plot_subplot_label(axH, fSet, 'A')
    BIC2plot = squeeze(allGLM_bin_pupil_BIC{3}(:,ibin,:));
    
    [~, selectedModel] = min(mean(BIC2plot(:,modelIdx)));
    
    hold on
    hBar = bar(mean(BIC2plot(:,modelIdx)) - mean(BIC2plot(:,selectedModel)));
    
    % plot se
    for imodel = 1:size(BIC2plot(:,modelIdx),2)
        plot([imodel imodel],...
            [(mean(BIC2plot(:,imodel)) - mean(BIC2plot(:,selectedModel))) - (std(BIC2plot(:,imodel) - BIC2plot(:,selectedModel))/sqrt(nSub)) ...
            (mean(BIC2plot(:,imodel)) - mean(BIC2plot(:,selectedModel))) + (std(BIC2plot(:,imodel) - BIC2plot(:,selectedModel))/sqrt(nSub)) ],...
            'color', 'k','linewidth',fSet.LineWidth);
        
        mod2compare = BIC2plot(:, imodel) - BIC2plot(:,selectedModel);
        [P, H, STATS] = signrank(mod2compare, 0, 'method', 'exact');
        [pstring_chi,starstring] = getSignificanceStrings(P, 0, 1, '\it p ');
        text(imodel, max(hBar.YData) * 0.01, starstring, 'fontsize', fSet.Fontsize_title, 'Color', 'w', 'HorizontalAlignment', 'center')
    end
    
    hBar.FaceColor = [0.5 0.5 0.5];
    set(gca,'xtick', 1:size(BIC2plot(:,modelIdx),2))
    xlabel('Model number', 'fontsize', fSet.Fontsize_text)
    ylabel('\Delta BIC', 'fontsize', fSet.Fontsize_text)
    
    %%% effect size
    subtightplot(nbin,2,subplotIdx(ibin)+1,[0.08 0.08], [0.1 0.1], [0.1 0.1])
    figInit('ax');
    hold on
    
    tstat2plot = squeeze(allGLM_bin_pupil_Ramp_tstat(:,ibin,:));
    
    hBar = bar(mean(tstat2plot) );
    % plot se
    for icomponent = 1:size(tstat2plot,2)
        plot([icomponent icomponent],...
            [mean(tstat2plot(:,icomponent)) - (std(tstat2plot(:,icomponent))/sqrt(nSub)) ...
            mean(tstat2plot(:,icomponent)) + (std(tstat2plot(:,icomponent))/sqrt(nSub))],...
            'color', 'k','linewidth',fSet.LineWidth);
        
        [P, H, STATS] = signrank(tstat2plot(:, icomponent), 0, 'method', 'exact');
        [pstring_chi,starstring] = getSignificanceStrings(P, 0, 1, '\it p ');
        text(icomponent, max(hBar.YData) * 0.01, starstring, 'fontsize', fSet.Fontsize_title, 'Color', 'w', 'HorizontalAlignment', 'center')
        
    end
    
    hBar.FaceColor = [0.5 0.5 0.5];
    set(gca,'xtick', 1:size(tstat2plot,2))
    if ibin == nbin
        set(gca,'xticklabel',{'Onset', 'Sustained', 'Response'}, 'xticklabelrotation', 25)
    end
    
    ylabel('Effect size (t)', 'fontsize', fSet.Fontsize_text)
    
%     title(modelName, 'fontsize', fSet.Fontsize_title)
    
    
end

saveFigName = ['GLM_bin_modelComparison' fileExt_GLM];

figSave(saveFigName, paths.fig, {'png', 'svg'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary figure 3. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

panelLabeldisplacement = [0.06 0.06];
subplotgap = [0.12 0.07];
subplotmargin = [0.1 0.1];


[figHandle, fSet] = figInit('fig', 3, {...
    'figureType', 'Manuscript'; ...
    'width', 28; ...
    'height', 18});

maxYLIM = [0 30];

ncol = 5;
nrow = 3;

for iplot = 1:2
    
    switch iplot
        case 1
            tmpSelectedModel = selectedModel;
            panelLabel = {'A', 'B'};
        case 2
            tmpSelectedModel = 1;
            panelLabel = {'C', 'D'};
    end
    
    switch tmpSelectedModel
        case 1
            VIF2plot_conc = allGLM_conc_pupil_StimResp_VIF;
            VIF2plot_st = allGLM_pupil_StimResp_VIF;
            
            labels = {'Target Onset', 'Response'};
            title2plot = 'Target-Response';
            
        case 3
            VIF2plot_conc = allGLM_conc_pupil_Ramp_VIF;
            VIF2plot_st = allGLM_pupil_Ramp_VIF;
            
            labels = {'Target Onset', 'Sustained', 'Response'};
            title2plot = 'Linear up-ramp';
                        
        otherwise
            error
            
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% VIF for across trial GLM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axH = subtightplot(nrow,ncol,[1 + 2 * (iplot-1)], subplotgap, subplotmargin, subplotmargin);
    plot_subplot_label(axH, fSet, panelLabel{1}, panelLabeldisplacement)
    figInit('ax');
    hold on
    
    hBar = bar(mean(VIF2plot_conc, 1));
    
    % plot se
    clear tmpsub tmpbin tmpcomp
    for icomp = 1:size(VIF2plot_conc,2)
        plot([icomp icomp],...
            [(squeeze(mean(VIF2plot_conc(:, icomp), 1)) - squeeze(std(VIF2plot_conc(:, icomp), [], 1))/sqrt(nSub)) ...
            (squeeze(mean(VIF2plot_conc(:, icomp), 1)) + squeeze(std(VIF2plot_conc(:, icomp), [], 1))/sqrt(nSub)) ] , ...
            'color', 'k','linewidth',fSet.LineWidth);
    end
    hBar.FaceColor = [0.5 0.5 0.5];
    
    ylim(maxYLIM)
    XLIM = get(gca,'xlim');
    xlim(XLIM + [0.6 -0.6])
    set(gca,'xtick', 1:size(VIF2plot_conc,2), 'xticklabel', labels, 'xticklabelrotation', 25, 'fontsize', fSet.Fontsize_text)
    ylabel('VIF', 'fontsize', fSet.Fontsize_text)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% single trial VIF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axH = subtightplot(nrow,ncol,[2 + 2 * (iplot-1)], subplotgap, subplotmargin, subplotmargin);
    plot_subplot_label(axH, fSet, panelLabel{2}, panelLabeldisplacement)
    figInit('ax');
    hold on
    
    subVIF = NaN(nSub, size(VIF2plot_st,2));
    for isub = 1:nSub
        
        trIdx = allSubject == isub;
        subVIF(isub,:) = nanmean(VIF2plot_st(trIdx, :));
        
    end
    hBar = bar([nanmean(subVIF, 1)]);
    
    % plot se
    clear tmpsub tmpbin tmpcomp
    for icomp = 1:size(VIF2plot_conc,2)
        plot([icomp icomp],...
            [(squeeze(nanmean(subVIF(:, icomp), 1)) - squeeze(std(subVIF(:, icomp), [], 1))/sqrt(nSub)) ...
            (squeeze(nanmean(subVIF(:, icomp), 1)) + squeeze(std(subVIF(:, icomp), [], 1))/sqrt(nSub)) ] , ...
            'color', 'k','linewidth',fSet.LineWidth);
    end
    hBar.FaceColor = [0.5 0.5 0.5];
    ylim(maxYLIM)
    xlim(XLIM + [0.6 -0.6])
    set(gca,'xtick', 1:size(VIF2plot_conc,2), 'xticklabel', labels, 'xticklabelrotation', 25, 'fontsize', fSet.Fontsize_text)
    ylabel('VIF', 'fontsize', fSet.Fontsize_text)
    % title({'Single-trial model:', title2plot}, 'fontsize', fSet.Fontsize_title)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot relationship between VIF and RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axH = subtightplot(nrow,ncol,5, subplotgap, subplotmargin, subplotmargin);
plot_subplot_label(axH, fSet, 'E', panelLabeldisplacement)
figInit('ax');

VIF2plot_st = allGLM_pupil_Ramp_VIF;

plot(allRT_zscore, VIF2plot_st(:,1), '.', 'MarkerSize', fSet.MarkerSize)

ylim([0 150])
xlabel('RT (z-score)', 'fontsize', fSet.Fontsize_text)
ylabel({'Target onset','component VIF'}, 'fontsize', fSet.Fontsize_text)
figInit('ax');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot pupil timecourse and behavioural performance for excludsion of VIF>5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Pupil
subplotIdx = [6 8 11];
for iplot = 1:3
    
    switch iplot
        case 1
            idx2plot        = strcmpi(allBin2use, 'GLM_pupil_Ramp_stim_VIF5_regress_bl_iti_side');
            cfg.ylim_RT     = [550 680];
            cfg.ylim_RT_CV  = [0.12 0.28];
        case 2
            idx2plot        = strcmpi(allBin2use, 'GLM_pupil_Ramp_sust_regress_bl_iti_side');
            cfg.ylim_RT     = [350 690];
            cfg.ylim_RT_CV  = [0.10 0.33];
        case 3
            idx2plot        = strcmpi(allBin2use, 'GLM_pupil_StimResp_stim_regress_bl_iti_side');
            cfg.ylim_RT     = [500 610];
            cfg.ylim_RT_CV  = [0.20 0.31];
    end
    filename_R = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{idx2plot} '_' bintype fileExt_preprocess  fileExt_CDT fileExt_GLM  ];

    cfg.ylim = [-0.2 0.3];
    cfg.xlim = [-250 2500];

    pupil2plot = pupil_bp{idx2plot};
    
    axH = subtightplot(nrow, ncol, subplotIdx(iplot), subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
    plot_subplot_label(axH, fSet, subplotIdx(iplot), panelLabeldisplacement)
    
    hold on
        
    for ibin = 1:nbin2use
        A = boundedline(cfg.t_pupil(1:size(pupil2plot,4)),squeeze(mean(pupil2plot(:,ibin,:,:),1)),squeeze(std(pupil2plot(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
        set(A,'linewidth',fSet.LineWidth_in)
    end
    plot([0 0], cfg.ylim ,'k','linewidth',1)
    
    xlim(cfg.xlim)
    ylim(cfg.ylim)
    
    xlabel('Time (ms)','fontsize',fSet.Fontsize_text);
    ylabel('Pupil diameter (a.u.)','fontsize',fSet.Fontsize_text);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% RT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    axH = subtightplot(nrow, ncol, subplotIdx(iplot)+1, subplotgap, subplotmargin, subplotmargin);
    pos = axH.Position;
    pos(1) = pos(1) - 0.015;
    axH.Position = pos;
    figInit('ax');
    plot_subplot_label(axH, fSet, subplotIdx(iplot)+1, panelLabeldisplacement + [-0.01 0])
    hold on
        
    [stats, mfit] = loadStatsR(filename_R, 'RT',nSub,nbin2use);
       
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
    plot(x, mean(RT{idx2plot}),plotMarker{1},...
        'markersize', fSet.MarkerSize,...
        'color', fSet.colors(1,:),'MarkerFaceColor', fSet.colors(1,:), 'MarkerEdgeColor', 'k');
    
    ylim([cfg.ylim_RT])
    h = ylabel('RT (ms)','fontsize',fSet.Fontsize_text);
    h.Color = fSet.colors(1,:);
    set(gca,'xtick',x)
    xlabel('Pupil Bin','fontsize',fSet.Fontsize_text);
%     axis square
        
    % text for RT_CV
    [stats, mfit] = loadStatsR(filename_R, 'RT_CV',nSub,nbin2use);
    
    ax1 = gca;
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
%     axis square
    subplot_ax1 = get(gca);
    
    xlim([0.5 nbin2use+0.5])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot rsquare for different single trial models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

example_sub = 22;
rsquare2plot1 = allGLM_pupil_StimResp_rsquare; BIC2plot1 = 1;
rsquare2plot2 = allGLM_pupil_Ramp_rsquare; BIC2plot2 = 3;
xticklabels = {'Target-Response','Ramp'};
        
axH = subtightplot(nrow,ncol,10, subplotgap, subplotmargin, subplotmargin);
plot_subplot_label(axH, fSet, 'J', panelLabeldisplacement + [-0.01 0])
figInit('ax');
hold on

xlabel('\Delta R^2', 'fontsize', fSet.Fontsize_text)

subRsquare = NaN(nSub, 2);
for isub = 1:nSub
    
    trIdx = allSubject == isub;
    subRsquare(isub,:) = nanmean([rsquare2plot1(trIdx, :) rsquare2plot2(trIdx, :)]);
    
    if isub == example_sub
        trIdx2use = trIdx;
    end
    
    % check whether signrank is different from 0
    [P, H, STATS] = signrank(rsquare2plot2(trIdx, :) - rsquare2plot1(trIdx, :), 0, 'method', 'exact','tail','right');
    
    if H == 0
        error
    end
end

% if bar
if 0
    hplot = bar([nanmean(subRsquare, 1)]);
    
    % plot se
    clear tmpsub tmpbin tmpcomp
    for icomp = 1:size(subRsquare,2)
        plot([icomp icomp],...
            [(squeeze(nanmean(subRsquare(:, icomp), 1)) - squeeze(std(subRsquare(:, icomp), [], 1))/sqrt(nSub)) ...
            (squeeze(nanmean(subRsquare(:, icomp), 1)) + squeeze(std(subRsquare(:, icomp), [], 1))/sqrt(nSub)) ] , ...
            'color', 'k','linewidth',fSet.LineWidth);
    end
    hplot.FaceColor = [0.5 0.5 0.5];
    
else
    
    hist(subRsquare(:,2) - subRsquare(:,1), 10);
    hplot = get(gca);
    hplot.Children.FaceColor = [0.5 0.5 0.5];
    xlim([0 0.1])
    
    ylabel('#Subjects', 'fontsize', fSet.Fontsize_text)
end

axes('position',[0.89 0.52 0.06 0.07]) ; % inset
figInit('ax',[],{'fontsize',10});
set(gca,'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4])
hold on

hold on
hist(rsquare2plot2(trIdx2use, :) - rsquare2plot1(trIdx2use, :), 10);
hplot = get(gca);
hplot.Children.FaceColor = [0.5 0.5 0.5];
meanDR2 = nanmean(rsquare2plot2(trIdx2use, :) - rsquare2plot1(trIdx2use, :));
YLIM = get(gca,'ylim');
plot(meanDR2, YLIM(2) * 0.95, 'v', 'markerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k')
ylabel('#Trials', 'fontsize', fSet.Fontsize_text_in)
xlabel('\Delta R^2', 'fontsize', fSet.Fontsize_text_in)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot fitted pupil diameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx2plot        = strcmpi(allBin2use, 'GLM_pupil_Ramp_stim_regress_bl_blPhase_iti_side');

cfg.t_pupilFit = cfg.t_pupil(cfg.t_pupil>=pupilSet.prestim_range(1));

cfg.ylim = [-0.2 0.4];
cfg.xlim = [-250 2500];

for iplot = 13:14
    if iplot == 13
        pupil2plot = pupil_bp;
        tt2use = cfg.t_pupil;
    elseif iplot == 14
        pupil2plot = pupil_GLM_Ramp_yhat;
        tt2use = cfg.t_pupilFit;
    end
    
    axH = subtightplot(nrow,ncol,iplot,subplotgap, subplotmargin, subplotmargin);
    plot_subplot_label(axH, fSet, iplot, panelLabeldisplacement)
    figInit('ax');
    hold on

    for ibin = 1:nbin2use
        A = boundedline(tt2use(1:size(pupil2plot{1},4)),squeeze(nanmean(pupil2plot{idx2plot}(:,ibin,:,:),1)),squeeze(nanstd(pupil2plot{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
        set(A,'linewidth',fSet.LineWidth_in)
    end
    plot([0 0], cfg.ylim ,'k','linewidth',1)
    
    xlabel('Time (ms)','fontsize',fSet.Fontsize_text);
    ylabel('Pupil diameter (a.u.)','fontsize',fSet.Fontsize_text);
    xlim(cfg.xlim)
    ylim(cfg.ylim)
end

saveFigName = ['GLM_VIF' fileExt_GLM];
figSave(saveFigName, paths.fig, {'png', 'svg'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Figure 4. Plot the average weight for target onset component for binned GLM scores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[figHandle, fSet] = figInit('fig', 2, {...
    'figureType', 'Manuscript'; ...
    'width', 24; ...
    'height', 16});

subplotgap = [0.08 0.06];
subplotmargin = [0.1 0.1];

componentIdx = 1;
ncol = 4;
nrow = 3;
subplotIdx = [1 5 9 3];


if pupilSet.orthogonalise_predictors==0
    plots = 1:3;
else
    plots = 4;
end
for iplot = plots

    if iplot==1
        dat2use = allGLM_bin_pupil_Ramp_beta{3};
        vif2use = allGLM_bin_pupil_Ramp_VIF{3};
    elseif iplot==2
        dat2use = allGLM_bin_pupil_Ramp_beta{2};
        vif2use = allGLM_bin_pupil_Ramp_VIF{2};
    elseif iplot==3
        dat2use = allGLM_bin_pupil_Ramp_beta{1};
        vif2use = allGLM_bin_pupil_Ramp_VIF{1};
    elseif iplot==4
        dat2use = allGLM_bin_pupil_Ramp_beta{3};
        vif2use = allGLM_bin_pupil_Ramp_VIF{3};
    end

    axH = subtightplot(nrow,ncol,[subplotIdx(iplot) ], subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
    plot_subplot_label(axH, fSet, iplot)

    hold on
    hBar = bar(squeeze(mean(dat2use(:, :, componentIdx), 1)));
    % plot( 1:size(allGLM_conc_pupil_BIC,2), std(allGLM_conc_pupil_BIC(:,imodel) - allGLM_conc_pupil_BIC(:,selectedModel)))
    
    % plot se
    clear tmpsub tmpbin tmpcomp
    for ibin = 1:size(dat2use,2)
        plot([ibin ibin],...
            [(squeeze(mean(dat2use(:, ibin, componentIdx), 1)) - squeeze(std(dat2use(:, ibin, componentIdx), [], 1))/sqrt(nSub)) ...
            (squeeze(mean(dat2use(:, ibin, componentIdx), 1)) + squeeze(std(dat2use(:, ibin, componentIdx), [], 1))/sqrt(nSub)) ] , ...
            'color', 'k','linewidth',fSet.LineWidth);
        
%         [P, H, STATS] = signrank(dat2use(:, ibin, componentIdx), 0, 'method', 'exact');
%         [pstring_chi,starstring] = getSignificanceStrings(P, 0, 1, '\it p ');
%         text(ibin, max(hBar.YData) * 0.02, starstring, 'fontsize', fSet.Fontsize_title + 3, 'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', 'Ariel')

        tmpsub(:,ibin) = 1:nSub;
        tmpbin(:,ibin) = ones(nSub, 1) * ibin;
        tmpcomp(:,ibin) = squeeze(dat2use(:, ibin, componentIdx));
    end
    
    % test trend with lme
    lme_table = table(tmpsub(:), tmpbin(:), tmpcomp(:), 'VariableNames',{'Subject','Bin','Y'});
    
    [stats2report, meanFit, CIFit, SEFit, STATS, model] = fitlme_singleVar_sequential(lme_table, 0);
    
    % plot model fit
    plot(1:size(dat2use,2), meanFit, 'Color', [0.3 0.3 0.3], 'linewidth', fSet.LineWidth)
    %     h = patch(([(1:5) flip(1:5)]), [meanFit flip(meanFit)] + [-std(fit_Y)/sqrt(nSub) flip(std(fit_Y)/sqrt(nSub))], 'k');
    h = patch(SEFit(1,:), SEFit(2,:), 'k');
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
    
    switch model
        case 'linear'
            betaString = ['\beta_{1} = ' num2str(STATS.Estimate, '%1.1e')];
        case 'quadratic'
            betaString = ['\beta_{2} = ' num2str(STATS.Estimate, '%1.2f')];
        otherwise
            betaString = [];
    end
    
%     betaString = ['\beta_{1} = ' num2str(STATS.Estimate, '%1.2f')];
    [pstring_chi,starstring] = getSignificanceStrings(STATS.pValue, 0, 1, '\it p ');
    
    YLIM = get(gca,'ylim');
    ylim(YLIM.*[1 1.3])
    XLIM = get(gca,'xlim');
    h = text(mean(XLIM)*0.1, max(mean(dat2use(:, :, componentIdx), 1)) * 1.4, {betaString, pstring_chi}, 'fontsize',fSet.Fontsize_text_in, 'color', 'k');
    hBar.FaceColor = [0.5 0.5 0.5];
    
    set(gca,'xtick', 1:size(dat2use,2))
    xlabel('RT bin', 'fontsize', fSet.Fontsize_text)
    ylabel('Target onset \beta weight', 'fontsize', fSet.Fontsize_text)
    
    
    %%% VIF values for binned GLM
    axH = subtightplot(nrow,ncol,[subplotIdx(iplot)+1 ], subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
    hold on

    allVIF = squeeze(mean(vif2use,2));
    
    hBar = bar(squeeze(mean(allVIF, 1)));
    hBar.FaceColor = [0.5 0.5 0.5];
    % plot se
    clear tmpsub tmpbin tmpcomp
    for icomp = 1:size(allVIF,2)
        plot([icomp icomp],...
            [(squeeze(mean(allVIF(:,icomp), 1)) - squeeze(std(allVIF(:,icomp), [], 1))/sqrt(nSub)) ...
            (squeeze(mean(allVIF(:,icomp), 1)) + squeeze(std(allVIF(:,icomp), [], 1))/sqrt(nSub)) ] , ...
            'color', 'k','linewidth',fSet.LineWidth);
    end
    
    set(gca,'xtick', 1:size(allVIF,2))
    set(gca,'xticklabel',{'Onset', 'Sustained', 'Response'}, 'xticklabelrotation', 25, 'fontsize', fSet.Fontsize_text)
    ylabel('VIF', 'fontsize', fSet.Fontsize_text)

    ylim([0 20])
    
end

if pupilSet.orthogonalise_predictors
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Pupil
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    axH = subtightplot(nrow,ncol,[subplotIdx(iplot)+4], subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
    plot_subplot_label(axH, fSet, iplot+1)
    hold on
    
    idx2plot        = strcmpi(allBin2use, 'GLM_pupil_Ramp_stim_regress_bl_blPhase_iti_side');
    
    cfg.xlim = [-250 3000];
    cfg.ylim = [-0.2 0.4];
     
    for ibin = 1:nbin2use
        A = boundedline(cfg.t_pupil(1:size(pupil_bp{idx2plot},4)),squeeze(mean(pupil_bp{idx2plot}(:,ibin,:,:),1)),squeeze(std(pupil_bp{idx2plot}(:,ibin,:,:),[],1))/sqrt(nSub),'cmap',fSet.colors(ibin,:),'alpha');
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
    
    axH = subtightplot(nrow,ncol,[subplotIdx(iplot)+5], subplotgap, subplotmargin, subplotmargin);
    figInit('ax');
     
    cfg.ylim_RT     = [470 610];
    cfg.ylim_RT_CV  = [0.18 0.29];
    
    filename_R = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' allBin2use{idx2plot} '_' bintype fileExt_preprocess  fileExt_CDT fileExt_GLM  ];
    [stats, mfit] = loadStatsR(filename_R, 'RT',nSub,nbin2use);
    [pstring_chi,starstring] = getSignificanceStrings(stats.p, 0, 1, 'p ');
        
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
        otherwise
            betaString = [];
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

% h = suptitle('Stimulus onset component amplitude per RT bin');
% h.FontSize = fSet.Fontsize_title;
saveFigName = ['GLM_bin_componentComparison_beta' fileExt_GLM];

figSave(saveFigName, paths.fig, {'png', 'svg'})



%% check n trials rejected when VIF > threshold

nTrialsRej = NaN(nSub,2);

for isub = 1:nSub
    
    trIdx = (allSubject == isub) & allvalidtr_neg100_RT_200;
    nTrialsRej(isub, 1) = length(find(trIdx));
    nTrialsRej(isub, 2) = length(find((allGLM_pupil_Ramp_VIF(trIdx,1) >= 5)));
end

percentageRejected = nTrialsRej(:,2)./nTrialsRej(:,1) * 100;

fprintf('%1.2f %% +- %1.2f trials rejected\n', mean(percentageRejected), std(percentageRejected)/sqrt(nSub))

%% check n trials in binned GLM
nbin2use = 5;
nTrialsBin = NaN(nSub,nbin2use);

for isub = 1:nSub
    
    trIdx = (allSubject == isub) & allvalidtr_neg100_RT_200;
    
    [binEdges, binning2use, ~] = binVal(allRT(trIdx), nbin2use, 'equal');
    
    for ibin = 1:nbin2use
        nTrialsBin(isub, ibin) = length(find(binning2use==ibin));
    end
end

fprintf('%1.2f +- %1.2f trials per bin\n', mean(mean(nTrialsBin, 2),1), std(squeeze(mean(nTrialsBin, 2)), [], 1)/sqrt(nSub))



