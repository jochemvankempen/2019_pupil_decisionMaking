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

%% regress out phase of pupil response
pupilResp = nan(size(allGLM_pupil_Ramp_stim_beta));
for isub = 1:nSub
    trIdx = (allSubject==isub) & allvalidtr_neg100_RT_200;
    designM = [ones(size(allGLM_pupil_Ramp_stim_beta(trIdx))) cos(allPupil_bp_baselinePhase(trIdx)) sin(allPupil_bp_baselinePhase(trIdx))];
    [b, ~, resid] = regress(allGLM_pupil_Ramp_stim_beta(trIdx), designM);
    pupilResp(trIdx) = resid;
end

%% grab some single trial data, and write to csv

allRT_log(allRT_log==-inf)=0; % correct for infinity result from log transform

disp('writing st matrix')

varNames = {...
    'Subject', 'allvalidtr_neg500_0', 'allvalidtr_neg100_RT_200',...
    'Outcome', 'StimLoc', 'ITI', ...
    'Trial',  ...
    'blockTrial', ...
    'RT', ...
    'pupilBaseline','pupilResp'};

MAT = table(...
    allSubject, allvalidtr_neg500_0,   allvalidtr_neg100_RT_200, ... 
    allHits, allSideStimtr, allItitr, ...
    scaleVar(allTrial, 'minmax'), ...
    scaleVar(allBlockTrial, 'minmax'), ...
    scaleVar(allRT, 'minmax'), ...
    scaleVar(allPupil_lp_baseline, 'minmax'), ...
    scaleVar(pupilResp, 'minmax'), ...
    'VariableNames',varNames);

filename_csv = [paths.pop 'allSub_singleTrial' [fileExt_preprocess] '.csv'];
writetable(MAT,filename_csv)

%% extract some info

% how many trials rejected?
trRej = cell(1,2);
for itrwin = 1:2
    
    switch itrwin
        case 1
            trwin2check = 'allvalidtr_neg500_0';
        case 2
            trwin2check = 'allvalidtr_neg100_RT_200';
    end
    
    for isub = 1:nSub
        trIdx = MAT.Subject == isub;
        eval(['trRej{itrwin}(isub, 1) = sum( MAT.' trwin2check ' (trIdx, :) == 0);']);
        
    end
    
    
    fprintf('number of trials removed: %1.2f +- %1.2f\n', mean(trRej{itrwin}), std(trRej{itrwin})/sqrt(nSub))
    
end




