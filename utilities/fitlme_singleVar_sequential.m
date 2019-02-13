function [stats2report, meanFit, CIFit, SEFit, STATS, model] = fitlme_singleVar_sequential(dataMat, plotFit)
% sequentially fit LME to data and compare whether linear or quadratic fit are better
% than constant or linear fit, respectively.
% 
% INPUT:
% dataMat: table with all necessary variables (see table.m)
% Y: string that indicates the response variable, this has to be a column in dataMat
% plotFit: plot the fit of the model, default is 0
%
% OUTPUT:
% stats2report: Vector with 
%     1) loglikelihood of model comparison, 
%     2) delta DF: the difference in degrees of freedom in the models that were compared
%     3) p value of model comparison
% meanFit: the mean fit of the model
% CIFit: the confidence intervals of the model fit
% CIFit: the standard error intervals of the model fit
% STATS: stats on the model, the raw beta values etc. 
%
% Example:
% datamat = table(RT, Subject, Bin, 'VariableNames', {'Y','Subject','Bin'});
% [stats2report, meanFit, CIFit, SEFit, STATS, model] = fitlme_singleVar_sequential(datamat, 0);
%
% Jochem van Kempen, 26/06/2018

if nargin < 2
    plotFit = 0;
end

nSubject = length(unique(dataMat.Subject));
nBin = length(unique(dataMat.Bin));

onlyLinearFit = 0;
if nBin < 3
    onlyLinearFit = 1;
end
% dataMat.Bin = categorical(dataMat.Bin);

lme = fitlme(dataMat,['Y ~ 1 + (1|Subject)']);% fit constant
lme2 = fitlme(dataMat,['Y ~ Bin + (1|Subject)']); % fit linear model
if ~onlyLinearFit
    lme3 = fitlme(dataMat,['Y ~ Bin^2 + (1|Subject)']); % fit quadratic model (this includes the linear model), note this doesn't have orthogonal contrasts!!
end
% compare the model fits
comp1 = compare(lme, lme2, 'CheckNesting',true); % compare linear to constant
if ~onlyLinearFit
    comp2 = compare(lme2, lme3, 'CheckNesting',true); % compare quadratic to linear
else
    comp2.pValue = 1;
end

if comp2.pValue < 0.05
    [fit, fitCI] = predict(lme3); % get the fitted values for each subject  
    [~,~,STATS] = fixedEffects(lme3);
    
    stats2report = [comp2.LRStat(2) comp2.deltaDF(2) comp2.pValue(2)];
    
    STATS = STATS(3,:);

    model = 'quadratic';
elseif comp1.pValue < 0.05
    [fit, fitCI] = predict(lme2); % get the fitted values for each subject
    [~,~,STATS] = fixedEffects(lme2);
    stats2report = [comp1.LRStat(2) comp1.deltaDF(2) comp1.pValue(2)];
    STATS = STATS(2,:);
    model = 'linear';
else
    stats2report = [comp1.LRStat(2) comp1.deltaDF(2) comp1.pValue(2)];
    meanFit = [];
    CIFit = []; SEFit = []; STATS = [];
    model = 'none';
    return
end

fit_Y = reshape(fit, [nSubject,nBin]);% reshape back into subject x Bin
fit_Y_CI1 = reshape(fitCI(:,1), [nSubject,nBin]);% reshape back into subject x Bin
fit_Y_CI2 = reshape(fitCI(:,2), [nSubject,nBin]);% reshape back into subject x Bin

meanFit = mean(fit_Y);
CIFit = [1:nBin nBin:-1:1; mean(fit_Y_CI1) flip(mean(fit_Y_CI2))];
SEFit = [1:nBin nBin:-1:1; meanFit - std(fit_Y)/sqrt(nSubject) flip(meanFit + std(fit_Y)/sqrt(nSubject))];

if plotFit
    var2plot = reshape(eval(['dataMat.' Y ]), [nSubject,nBin]);
    figure
    hold on
    h = plot(meanFit, 'linewidth', 2);
    h = patch(CIFit(1,:), CIFit(2,:), h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
    plot(mean(var2plot),'.','markersize',15)
end






