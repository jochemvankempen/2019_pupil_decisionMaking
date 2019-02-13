function [stats, fit, mfit] = loadStatsR(filename, type, nSub, nBin2use, forceMod)
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
% Load data from csv file created by R script bigDots_stats

if nargin<5
    forceMod = [];
end
SS = CSVSheets(CSVSheets,'filepath',[filename '_R_statistics.csv']);
SS = load(SS);
R_stats = convert2struct(SS);

idx     = strcmpi(R_stats.key,type);
mfit.L.df       = R_stats.Df_L(idx) - R_stats.Df_Intercept(idx);
mfit.L.LRatio   = R_stats.Chisq_L(idx);
mfit.L.p        = R_stats.P_L(idx);
mfit.L.B        = R_stats.BinL_Estimate_Trend(idx);
mfit.L.B_SE     = R_stats.BinL_SE_Trend(idx);

mfit.Q.df       = R_stats.Df_Q(idx) - R_stats.Df_L(idx);
mfit.Q.LRatio   = R_stats.Chisq_Q(idx);
mfit.Q.p        = R_stats.P_Q(idx);
mfit.Q.B        = R_stats.BinQ_Estimate_Trend(idx);
mfit.Q.B_SE     = R_stats.BinQ_SE_Trend(idx);

[~,fitLIdx] = grep(fields(R_stats), 'fitL_');
[~,fitQIdx] = grep(fields(R_stats), 'fitQ_');

SNames = fieldnames(R_stats); 
LNames = SNames(fitLIdx);
QNames = SNames(fitQIdx);
tmpFitL = zeros(numel(LNames) ,1);
tmpFitQ = zeros(numel(LNames) ,1);
for loopIndex = 1:numel(LNames) 
    tmpFitL(loopIndex) = R_stats.(LNames{loopIndex})(idx);
    tmpFitQ(loopIndex) = R_stats.(QNames{loopIndex})(idx);
end 

fitL_mat = reshape(tmpFitL, [nSub, nBin2use]);
fitQ_mat = reshape(tmpFitQ, [nSub, nBin2use]);

fitL.mean = nanmean(fitL_mat);
fitL.se = nanstd(fitL_mat,[],1)/sqrt(nSub);

fitQ.mean = nanmean(fitQ_mat);
fitQ.se = nanstd(fitQ_mat)/sqrt(nSub);

if ~isempty(forceMod)
    switch forceMod
        case 'Linear'
            stats = mfit.L;
            stats.model = 'Linear';
            fit = fitL;
        case 'Quadratic'
            stats = mfit.Q;
            stats.model = 'Quadratic';
            stats.U = R_stats.U(idx);
            fit = fitQ;
    end
    
else
    if mfit.Q.p < 0.05
        stats = mfit.Q;
        stats.model = 'Quadratic';
        stats.U = R_stats.U(idx);
        fit = fitQ;
    elseif mfit.L.p < 0.05
        stats = mfit.L;
        stats.model = 'Linear';
        fit = fitL;
    else
        stats = mfit.L;
        stats.model = 'NA';
        fit = [];
    end
end


