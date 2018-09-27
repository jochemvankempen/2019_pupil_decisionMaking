function batch_bigDots(serverSubIdx)
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
%
% Input: serverSubIdx. This script is adapted to run on linux cluster,
% performing 

%%
clear all % comment out if running on server
close all
clc

if isdir('/home/vjochem/Monash126/Jochem/bigDots/') % massive cluster
    addpath(genpath('/home/vjochem/Monash126/Jochem/repositories/'))
    cd('/home/vjochem/Monash126/Jochem/toolboxes/eeglab13_5_4b/')
    eeglab
end

% prevent PC from sleeping when running computations
alterPerformanceMode('highPerformance')

%% Set jobs to execute

Job.batch_preprocess    = 1;
Job.analyse_CDT         = 1;

dataset = 'bigDots'

%% get subject data, paths, and analysis settings

getSubSpecs_bigDots
getfilenames_bigDots

% chanlocs = readlocs('actiCAP65_ThetaPhi.elp','filetype','besa'); % brainproducts channel configuration
% chanlocs = readlocs('cap64.loc'); %biosemi
% [~,plot_chans, exclude_chans] = getChannelName;
% figure; topoplot([],chanlocs,'style','blank','electrodes','labelpoint','chaninfo',plot_chans);

single_participants = [2:81];
nSub = length(single_participants);

setAnalysisSettings_bigDots

% add extension to file name with settings set in setAnalysisSettings_bigDots
fileExt = ['_CSD(' num2str(doCSD) ')' ... 
    '_chAlpha(' num2str(nch_alpha) ')'];

fileExt

%% preprocess
if Job.batch_preprocess
    parfor isub = 1:nSub
%     for isub = serverSubIdx
        subIdx  = single_participants(isub); %get matching subject
        disp(allsubj{subIdx})
        
        blocks  = allblocks{subIdx};
        badchans = allbadchans{subIdx};
        if ismember(subject_folder{subIdx},TCD_bigdots)
            batch_preprocess(allsubj, subIdx, paths, files, blocks, badchans, 'TCD',dataset,'big_dots_erp_final')
        elseif ismember(subject_folder{subIdx},Monash_bigdots)
            batch_preprocess(allsubj, subIdx, paths, files, blocks, badchans, 'Monash',dataset,'big_dots_erp_final')
        end
    end
end

%% single trial analysis
if Job.analyse_CDT
    parfor isub = 1:nSub
%     for isub = serverSubIdx
        analyse_CDT(isub, single_participants, dataset, [fileExt '_final'])
        
    end
end

%% exit script if run on massive
if serverSubIdx
    return
end

%% create single trial matrix, all subjects
analyse_CDT_fileExt = ['_final'];
if 0
    collect_data
else
    filename_mat = [paths.pop 'allSub_singleTrial' fileExt analyse_CDT_fileExt '.mat'];
    load([filename_mat])
end

alterPerformanceMode('balanced')

%% write single trial matrix to csv

write_st_matrix_csv

%% collect p_level
% sort and average data, write out to csv file for R analysis

% settings for data sorting
sideInfo        = 1; % neglect which side the target was presented (regress this out from variable used for sorting)
nbin2use        = 5; % 5 is good amount for quadratic relationships
bintype         = 'equal'; % equal sized bins

clear bin2use allBin2use

allBin2use = {...
    'pupil_lp_baseline_regress_iti_side'; ...
%     'pupil_bp_baseline_regress_iti_side'; ...
%     'pupil_lp_linproj_resp_locked_neg500_200_regress_bl_iti_side'; ...
%     'pupil_bp_linproj_resp_locked_neg500_200_regress_bl_iti_side'; ...
    'pupil_lp_RT_neg200_200_regress_bl_iti_side'; ...
%     'pupil_bp_RT_neg200_200_regress_bl_iti_side'; ...
    'pretarget_alpha'; ...
    'pretarget_alpha_asym'; ...
    'N2i_amplitude_regress_iti_side'; ...

    };

write_participant_level_csv

%% plot figures (run write_participant_level_csv, and R analysis scripts first)

plot_pLevel_bigDots_manuscript_final

%% print stats on screen for easy manuscript writing, uses stats computed in R
analyse_CDT_fileExt = ['_final'];
sideInfo        = 1;
nbin2use        = 5;
bintype         = 'equal';

allBin2usePrint = {...
%     'pupil_lp_baseline_regress_iti_side'; ...
    'pupil_bp_baseline_regress_iti_side'; ...
%     'pupil_lp_linproj_resp_locked_neg500_100_regress_bl_iti_side'; ...
%     'pupil_bp_linproj_resp_locked_neg500_100_regress_bl_iti_side'; ...
%     'pupil_lp_RT_neg200_200_regress_bl_iti_side'; ...
    'pupil_bp_RT_neg200_200_regress_bl_iti_side'; ...
%     'pretarget_alpha'; ...
%     'pretarget_alpha_asym'; ...
%     'N2i_amplitude_regress_iti_side'; ...

    };

allTypes = {'RT', 'RT_CV', ...
    'alpha','alpha_asym',...
    'N2c_latency', 'N2c_amplitude','N2i_latency', 'N2i_amplitude', 'N2c_ITPC', 'N2i_ITPC', ...
    'CPP_onset', 'CPP_slope2', 'CPP_slope_var', 'CPPr_amplitude', 'CPPr_amplitude_var', 'CPP_ITPC', ...
    'preRespBeta_slope', 'preRespBeta', 'preRespBeta_base','preTargetBeta', ...
    'plevelPupil_RT200'};

for ibinType = 1:length(allBin2usePrint)
    fprintf('\n')
    clear statsTable
    clear df_lme LRatio p2use
    bin2use = allBin2usePrint{ibinType};
    
    filename = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
    
    for itype = 1:length(allTypes)
        type = allTypes{itype};
        
        [stats] = loadStatsR(filename, type, nSub, nbin2use);
        [pstring_chi] = getSignificanceStrings(stats.p, 0, 0);
        
        fprintf('TYPE:%s. Chi(%d) = %1.2f, p %s \n', type, stats.df, stats.LRatio, pstring_chi)
        switch stats.model
            case 'Linear'
                fprintf('TYPE:%s. %s. b = %1.2f, SE = %1.2f \n', type, stats.model, stats.B, stats.B_SE)
            case 'Quadratic'
                fprintf('TYPE:%s. %s. b2 = %1.2f, SE = %1.2f \n', type, stats.model, stats.B, stats.B_SE)
        end
        
        
    end
end




