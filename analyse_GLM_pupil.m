function analyse_GLM_pupil(isub, single_participants, dataset, fileExt_CDT, fileExt_GLM)
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
% load all relevant data and execute specific GLM analysis set below

%%
tic
switch dataset
    case 'CD'
        getSubSpecs_CD
        getfilenames_CD
    case 'bigDots'
        getSubSpecs_bigDots
        getfilenames_bigDots
end

setAnalysisSettings_bigDots

%% get filenames and drug information, if applicable

switch dataset
    case 'CD'
        filenames = [paths.s(single_participants(isub)).readbase allsubj{single_participants(isub)} fileExt_CDT '_ST.mat'];
    case 'bigDots'
        filenames = [paths.s(single_participants(isub)).readbase allsubj{single_participants(isub)} fileExt_CDT '_ST.mat'];
end
disp(['loading: ' filenames])
load([filenames], 'pupil', 'validtr', 'subRT');

%% z-score pupil
if pupilSet.zscore
    pupil.bp = zscore(pupil.bp')';
end
%% perform GLM

[st_GLMfit, conc_GLMfit, bin_GLMfit] = analysis_GLM_pupil_singleTrial(pupil.bp, pupilSet, subRT, validtr.neg100_RT_200, t_pupil, fs);

%% save
savefilename = [paths.s(single_participants(isub)).savebase allsubj{single_participants(isub)} fileExt_GLM '.mat'];

save(savefilename, ...
    'st_GLMfit','conc_GLMfit','bin_GLMfit',...
    '-v7.3')
toc
