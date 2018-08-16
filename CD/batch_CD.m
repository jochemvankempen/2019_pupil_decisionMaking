% function batch_CD(massiveSubIdx)
% These scripts reproduce the analysis in the paper: van Kempen et al.,
% (2018) 'Behavioural and neural signatures of perceptual evidence
% accumulation are modulated by pupil-linked arousal'. 
% 
% Many of these scripts are based on the original scripts for the papers
% Newman et al. (2017), Journal of Neuroscience.
% https://github.com/gerontium/big_dots, and Loughnane et al. (2018), The
% Journal of Neuroscience
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
%
% Input: serverSubIdx. This script is adapted to run on linux cluster,
% performing 

clear all
close all
clc

if isdir('/home/vjochem/Monash126/Jochem/bigDots/') % massive
    addpath(genpath('/home/vjochem/Monash126/Jochem/repositories/'))
    cd('/home/vjochem/Monash126/Jochem/toolboxes/eeglab13_5_4b/')
    eeglab
end
    
%%

Job.batch_preprocess_CD = 1;
Job.analyse_CDT = 1;

dataset = 'CD';
upDown = {'U','D'};
%%
getSubSpecs_CD
getfilenames_CD

chanlocs = readlocs('JBhead96_sym.loc'); %biosemi
plot_chans = 1:96;

single_participants = [1:2 4:17 19];

nSub = length(single_participants);

setAnalysisSettings_bigDots % use same settings as bigDots

fileExt = ['_CSD(' num2str(doCSD) ')_chAlpha(' num2str(nch_alpha) ')'];
fileExt
%% preprocess
if Job.batch_preprocess_CD
    parfor isub=1:nSub
        subIdx  = single_participants(isub); %get matching subject
        disp(allsubj{subIdx})
        
        for iside = 1:2 %up, down contrast change
            fileExtDir  = upDown{iside};
            blocks      = allblocks{subIdx}{iside};
            badchans    = allbadchans{subIdx};
            
            batch_preprocess_CD(allsubj, subIdx, paths, squeeze(files(:,iside,:)), blocks, badchans, dataset, [fileExtDir '_final'])
        end
    end
end

%% Repeat of the same code as bigdots analysis below
%%% after preprocessing, the analysis scripts are exactly the same as the
%%% scripts for bigDots dataset. 
%%%
%%% If you want to analyse the EEG data for the CD paradigm, then settings
%%% need to be changed. 
%%% For example, the channel numbers/names don't correspond between the two
%%% datasets. Here I only analyzed the behavioural data, so this is not an
%%% issue. 

%% single trial analysis
if Job.analyse_CDT
    parfor isub = 1:nSub
        analyse_CDT(isub, single_participants, dataset, [fileExt '_final']) % use same analysis script as bigDots

    end
end

%% create single trial matrix 
analyse_CDT_fileExt = ['_final'];
if 0
    collect_data
else
    filename_mat = [paths.pop 'allSub_singleTrial' fileExt analyse_CDT_fileExt '.mat'];
    load([filename_mat])
end

%% write single trial matrix to csv

write_st_matrix_csv

%% collect p_level
% sort and average data, write out to csv file for R analysis

% settings for data sorting
normaliseData   = 0; % scale variable between 0 and 100 per subject
sideInfo        = 1; % neglect which side the target was presented (regress this out from variable used for sorting)
nbin2use        = 5; % 5 is good amount for quadratic relationships
bintype         = 'equal'; % equal sized bins

clear bin2use allBin2use

allBin2use = {...
    'pupil_lp_baseline_regress_iti_side'; ...
    'pupil_lp_RT_neg200_200_regress_bl_iti_side'; ...
    };

write_participant_level_csv





