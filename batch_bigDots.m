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
% Input: serverSubIdx. This script is adapted to run on linux cluster.
%%
clear all; serverSubIdx=[]; % comment out if running on server
close all
clc

if isdir('/home/vjochem/Monash126/Jochem/bigDots/') % massive cluster
    addpath(genpath('/home/vjochem/M/Jochem/repositories/'))
    cd('/home/vjochem/Monash126/Jochem/toolboxes/eeglab13_5_4b/')
    eeglab
end

% prevent PC from sleeping when running computations
alterPerformanceMode('highPerformance')

%% Set jobs to execute

Job.batch_preprocess    = 0;
Job.analyse_CDT         = 0;
Job.analyse_GLM_pupil   = 0;

dataset = 'bigDots'

readDataFromServer      = 0;

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
fileExt_preprocess = ['_pupil_0.01_6Hz'];
fileExt_CDT = ['_CSD(' num2str(doCSD) ')' ... 
    '_chAlpha(' num2str(nch_alpha) ')'...
    '_review'];

fileExt_GLM = ['_GLM_meanC(' num2str(pupilSet.meanCenter_predictors) ...
    ')_orthPr(' num2str(pupilSet.orthogonalise_predictors) ...
    ')_normPu(' num2str(pupilSet.normalise_pupil) ')'...
    '_bin(' pupilSet.bin2use ')',...
    '_bl_' num2str(abs(pupilSet.prestim_range(1)))];

fprintf('file extension preprocess: %s\n', fileExt_preprocess)
fprintf('file extension CDT: %s\n', fileExt_CDT)
fprintf('file extension GLM: %s\n', fileExt_GLM)

%% preprocess
if Job.batch_preprocess
    for isub = 1:nSub % comment out if running on server
        %     for isub = serverSubIdx
%         try
            subIdx  = single_participants(isub); %get matching subject
            disp(allsubj{subIdx})
            
            blocks  = allblocks{subIdx};
            badchans = allbadchans{subIdx};
            if ismember(subject_folder{subIdx},TCD_bigdots)
                batch_preprocess(allsubj, subIdx, paths, files, blocks, badchans, 'TCD',dataset,['big_dots_erp' fileExt_preprocess])
            elseif ismember(subject_folder{subIdx},Monash_bigdots)
                batch_preprocess(allsubj, subIdx, paths, files, blocks, badchans, 'Monash',dataset,['big_dots_erp' fileExt_preprocess])
            end
%         catch
%         end
    end
end

%% single trial analysis

if Job.analyse_CDT
    fprintf('Running single trial analysis\n')
    for isub = 1:nSub % comment out if running on server
%     for isub = serverSubIdx
        analyse_CDT(isub, single_participants, dataset, fileExt_preprocess, [fileExt_CDT])
        
    end
    if isempty(serverSubIdx)
        collect_data % do not run this script if running on server
    end
else
    filename_mat = [paths.pop 'allSub_singleTrial' fileExt_preprocess fileExt_CDT '.mat'];
    fprintf('Loading %s\n', filename_mat)
    load([filename_mat])
end

%% GLM analysis

if Job.analyse_GLM_pupil
    fprintf('Running GLM pupil analysis\n')
    for isub = 1:nSub % comment out if running on server
        %     for isub = serverSubIdx
        analyse_GLM_pupil(isub, single_participants, dataset, [fileExt_preprocess fileExt_CDT], [fileExt_preprocess fileExt_GLM])
        
    end
    if isempty(serverSubIdx)
        collect_data_GLM  % do not run this script if running on server
    end
else
    filename_mat = [paths.pop 'allSub' fileExt_preprocess fileExt_GLM '.mat'];
    fprintf('Loading %s\n', filename_mat)
    load([filename_mat])
end

%% exit script if run on massive
if isempty(serverSubIdx)
    return
end

%%
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
    %%% baseline
    'pupil_lp_baseline_regress_iti_side'; ...
%     'pupil_lp_1Hz_baseline_regress_iti_side'; ...
%     'pupil_lp_baselineSlope_regress_iti_side'; ...
%     'pupil_lp_baselineDiff_regress_iti_side'; ...
%     'pupil_bp_baseline_regress_iti_side'; ...
%     'pupil_bp_baselinePhase_regress_iti_side'; ...
%     'pupil_bp_baselinePhase_regress_iti_side_fix'; ...
    
    %%% pupil response
%     'pupil_bp_RT_neg200_200_regress_bl_iti_side';
%     'pupil_bp_RT_neg200_200_regress_bl_iti_side_RT';
%     'pupil_bp_average_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side';
%     'pupil_bp_slope_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side';
%     'pupil_bp_diff_maxDiff_pupilIRF_200_200_regress_bl_iti_side';
%     'pupil_bp_linearProjection_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side';
    
    %%% GLM
%     'GLM_pupil_Ramp_stim_regress_iti_side'; ...
%     'GLM_pupil_Ramp_stim_regress_bl_iti_side'; ...
%     'GLM_pupil_Ramp_stim_regress_blPhase_iti_side';
    'GLM_pupil_Ramp_stim_regress_bl_blPhase_iti_side'; ...
%     'GLM_pupil_Ramp_stim_VIF5_regress_bl_iti_side'; ...
%     'GLM_pupil_StimResp_stim_regress_bl_iti_side'; ...
%     'GLM_pupil_Ramp_sust_regress_iti_side'; ...
%     'GLM_pupil_Ramp_sust_regress_bl_iti_side'; ...
    
    %%% other
%     'RT_regress_iti_side'; ...
    'pretarget_alpha'; ...
    };

write_participant_level_csv




%% plot figures (run write_participant_level_csv, and R analysis scripts first)

plot_pLevel_bigDots_manuscript_final


%% plot GLM fit

plot_pLevel_bigDots_GLM


%% control plots

plot_pLevel_bigDots_extraControl

%% print stats on screen for easy manuscript writing, uses stats computed in R
sideInfo        = 1;
nbin2use        = 5;
bintype         = 'equal';

allBin2usePrint = {...
%     'pupil_lp_baseline_regress_iti_side'; ...
%     'pupil_bp_baseline_regress_iti_side'; ...

%     'RT_regress_iti_side'; ...
%     'GLM_pupil_Ramp_stim_regress_iti_side'; ...
%     'GLM_pupil_Ramp_stim_regress_bl_iti_side'; ...
    
%     'GLM_pupil_Ramp_sust_regress_iti_side'; ...
%     'GLM_pupil_Ramp_sust_regress_bl_iti_side'; ...
%     'GLM_pupil_Ramp_stim_regress_bl_blPhase_iti_side'; ...

    'pretarget_alpha'; ...
    };


allTypes = {'RT', 'RT_CV', ...
    'alpha',...
    'N2c_latency', 'N2c_amplitude','N2i_latency', 'N2i_amplitude', 'N2c_ITPC', 'N2i_ITPC', ...
    'CPP_onset', 'CPP_csd_slope2', 'CPPr_csd_amplitude', 'CPP_ITPC', ...
    'preRespBeta_slope', 'preRespBeta', 'preRespBeta_base', 'preRespBeta_baseAT', ...
    };

for ibinType = 1:length(allBin2usePrint)
    fprintf('\n')
    clear statsTable
    clear df_lme LRatio p2use
    bin2use = allBin2usePrint{ibinType};
    
    switch bin2use
        case {'pupil_lp_baseline_regress_iti_side','pupil_bp_baseline_regress_iti_side','pretarget_alpha'}
            filename = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt_preprocess fileExt_CDT  ];
        otherwise
            filename = [paths.pop 'participant_level_scaleVar(0)_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt_preprocess fileExt_CDT fileExt_GLM ];
    end
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




