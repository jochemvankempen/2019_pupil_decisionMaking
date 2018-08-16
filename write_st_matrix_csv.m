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
%% grab some single trial data, and write to csv

allRT_log(allRT_log==-inf)=0; % correct for infinity result from log transform

disp('writing st matrix')

varNames = {...
    'ID', 'Subject', 'allTrials', 'allvalidtr_neg500_RT_200', 'allvalidtr_neg100_RT_200',...
    'Outcome', 'StimLoc', 'ITI', ...
    'Trial',  ...
    'blockTrial', ...
    'RT', 'RT_log', ...
    'allPupil_lp_baseline','allPupil_bp_baseline','allPupil_lp_RT_neg200_200','allPupil_bp_RT_neg200_200'};

MAT = table(...
    ID, allSubject, (1:length(allSubject))', allvalidtr_neg500_RT_200,   allvalidtr_neg100_RT_200, ... 
    allHits, allSideStimtr, allItitr, ...
    scaleVar(allTrial, 'minmax'), ...
    scaleVar(allBlockTrial, 'minmax'), ...
    scaleVar(allRT, 'minmax'), ...
    scaleVar(allRT_log, 'minmax'), ...
    scaleVar(allPupil_lp_baseline, 'minmax'), ...
    scaleVar(allPupil_bp_baseline, 'minmax'), ...
    scaleVar(allPupil_lp_RT_neg200_200, 'minmax'), ...
    scaleVar(allPupil_bp_RT_neg200_200, 'minmax'), ...
    'VariableNames',varNames);

filename_csv = [paths.pop 'allSub_singleTrial' fileExt '.csv'];
writetable(MAT,filename_csv)

