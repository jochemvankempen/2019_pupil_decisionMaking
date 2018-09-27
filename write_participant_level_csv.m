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
%% participant level csv
% sort and average data, write out to csv file for R analysis

clear pl* CPP* RT* N2* alpha* beta* pupil* filename_csv

plotCPP = 0; % if 1, then onset of CPP will be plotted for each subject and bin. This is relevant for debugging, sometimes CPP onset latency is not found correctly
scaleVariables = 0; % if 1, then variables will be scaled between 0-1, this is used for final hierarchical regression analysis

% load some settings
setAnalysisSettings_bigDots
switch dataset
    case 'bigDots'
        getSubSpecs_bigDots
        getfilenames_bigDots
    case 'CD'
        getSubSpecs_CD
        getfilenames_CD
        
end

% the variables to be sorted
types2check = {...
    'RT','RT_CV','RT_window','RTcv_window',...
    'N2c_amplitude','N2c_latency','N2i_amplitude','N2i_latency',...
    'beta_base_response_amplitude_window','beta_pre_response_slope_window','CPPr_csd_slope','CPPr_amplitude','CPPr_csd_slope_window','CPPr_amplitude_window','CPP_onset',...
    'alpha_asym','alpha','alpha_preTarget_topo','alpha_asym_topo',...
    'beta_base_response', 'beta_response_amplitude', 'beta_base_response_amplitude','beta_pre_response_topo','beta_pre_response_slope',... 
    'pupil_lp_baseline', 'pupil_bp_baseline' , 'pupil_lp', 'pupil_bp', 'pupilr_lp', 'pupil_lp_RT_neg200_200',...
    'CPP','CPPr_csd','CPP_topo'...
    'N2c','N2i','N2c_topo',...
    'N2c_ITPC_bar','N2i_ITPC_bar','CPP_ITPC_bar',...
    'N2c_ITPC_band','N2i_ITPC_band','CPP_ITPC_band',...
    'N2c_ITPC','N2i_ITPC','CPP_ITPC',...
    'CPP_ITPC_window'};


for ibinType = 1:length(allBin2use)
    clear plevel*
    
    bin2use = allBin2use{ibinType};
    for itype = 1:length(types2check)
        clear type plot2make pVal normaliseData bar2plot topo2plot erp2plot asym2plot CPP_side_slope2 CPP_side_slope1 CPP_side_onsets plot2make
        
        type = types2check{itype};
        
        switch type
            case {'CPP_onset','CPPr_csd_slope','CPPr_csd_slope_var','CPPr_csd_slope_window'}
                datatype       = 'CPP';
            case {'CPP','CPPr_csd','N2c','N2i','pupil_lp','pupil_bp','pupilr_lp','beta','beta_response','beta_base_response','RT_window','RTcv_window','pupil_lp_baseline_window','pupil_lp_response_window','CPPr_amplitude_window','beta_base_response_amplitude_window'}
                datatype       = 'erp';%,'beta_base_response'
            case {'CPP_topo','alpha_preTarget_topo','alpha_asym_topo','N2c_topo','beta_base_postTarget_topo','beta_pre_response_topo'}
                datatype       = 'topo';
            case {'N2c_ITPC_bar','N2i_ITPC_bar','CPP_ITPC_bar'}
                datatype       = 'ITPC_bar'; %
            case {'N2c_ITPC_band','N2i_ITPC_band','CPP_ITPC_band'}
                datatype       = 'ITPC_band'; %
            case {'N2c_ITPC','N2i_ITPC','CPP_ITPC'}
                datatype       = 'ITPC'; %
            case 'CPP_ITPC_window'
                datatype       = 'ITPC_window';
            otherwise
                datatype       = 'bar';
        end
        
        
        binData_bigDots
        eval([types2check{itype} '{ibinType} = binnedData;'])
        
    end
    
    
    
    %% scale variables for regression analysis, or not (set in scaleVariables)

    plevel_Subject                       = reshape(plevelSub                                     , size(plevelSub,1)*size(plevelSub,2)*size(plevelSub,3), 1);
    plevel_Side                          = reshape(plevelSide                                    , length(plevel_Subject), 1);
    plevel_Bin                           = reshape(plevelBin                                     , length(plevel_Subject), 1);
    plevel_nTrial                        = reshape(plevelnTrial                                  , length(plevel_Subject), 1);

    for itype = 1:length(types2check)
        if eval(['ndims(' types2check{itype} '{ibinType} ) == 2;']) && eval(['size(' types2check{itype} '{ibinType} ,2) == nbin2use;'])
            if scaleVariables
                % yi=(xi-min xi) /(max xi-min xi). Across subjects
                scaletype = 'minmax';
                eval(['plevel_' types2check{itype} ' = scaleVar(reshape(' types2check{itype} ' {ibinType}  , length(plevel_Subject), 1), scaletype);']);
            else
                eval(['plevel_' types2check{itype} ' = reshape(' types2check{itype} ' {ibinType}  , length(plevel_Subject), 1);']);                
            end
        end
    end
    
    %     plevel_RT_window         = reshape(RT_window{ibinType}       , length(plevel_Subject), size(RT_window{ibinType},4));
    %     plevel_RTcv_window       = reshape(RTcv_window{ibinType}     , length(plevel_Subject), size(RTcv_window{ibinType},4));
    %     plevel_CPP               = reshape(CPP{ibinType}             , length(plevel_Subject), size(CPP{ibinType},4));
    %     plevel_CPPr              = reshape(CPPr_csd{ibinType}        , length(plevel_Subject), size(CPPr_csd{ibinType},4));
    %     plevel_CPP_ITPC_band     = reshape(CPP_ITPC_band{ibinType}   , length(plevel_Subject), size(CPP_ITPC_band{ibinType},4));
    %     plevel_N2c               = reshape(N2c{ibinType}             , length(plevel_Subject), size(N2c{ibinType},4));
    %     plevel_N2c_ITPC_band     = reshape(N2c_ITPC_band{ibinType}   , length(plevel_Subject), size(N2c_ITPC_band{ibinType},4));
    %     plevel_N2i_ITPC_band     = reshape(N2i_ITPC_band{ibinType}   , length(plevel_Subject), size(N2i_ITPC_band{ibinType},4));
    %     plevel_Pupil             = reshape(pupil_lp{ibinType}        , length(plevel_Subject), size(pupil_lp{ibinType},4));
    
    
    
    %% write out to csv
    varNames = {...
        'Subject'       , 'Side'            , 'Bin'       , 'nTrial', ...
        'RT'            , 'RT_CV'           ,   ...
        'CPP_onset'     , 'CPP_slope2'      , 'CPPr_amplitude', ...
        'N2c_amplitude' , 'N2c_latency'     , 'N2i_amplitude' , 'N2i_latency',...
        'alpha_asym'    , 'alpha'           , ...
        'preRespBeta'   , 'preRespBeta_base', 'preRespBeta_slope', ...
        'pupil_bl_lp'   , 'pupil_bl_bp'     , 'plevelPupil_RT200' ,...
        'N2c_ITPC'      , 'N2i_ITPC'        , 'CPP_ITPC'...
        };
    
    %         plevelRTmin5, plevelRTmin4, plevelRTmin3, plevelRTmin2, plevelRTmin1, plevelRT0, plevelRTplus1, plevelRTplus2, plevelRTplus3, plevelRTplus4, plevelRTplus5,
    MAT = table(...
        plevel_Subject               , plevel_Side                    , plevel_Bin                         , plevel_nTrial ,...
        plevel_RT                    , plevel_RT_CV                   , ...
        plevel_CPP_onset             , plevel_CPPr_csd_slope          , plevel_CPPr_amplitude,  ...
        plevel_N2c_amplitude         , plevel_N2c_latency             , plevel_N2i_amplitude               , plevel_N2i_latency,           ...
        plevel_alpha_asym            , plevel_alpha                   , ...
        plevel_beta_response_amplitude , plevel_beta_base_response_amplitude , plevel_beta_pre_response_slope , ...
        plevel_pupil_lp_baseline     , plevel_pupil_bp_baseline       , plevel_pupil_lp_RT_neg200_200, ...
        plevel_N2c_ITPC_bar          , plevel_N2i_ITPC_bar            , plevel_CPP_ITPC_bar            , ...
        'VariableNames',varNames);
    
    filename_csv{ibinType} = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype fileExt analyse_CDT_fileExt  ];
    % keyboard
    writetable(MAT,[filename_csv{ibinType} '.csv'])
    %
    
    
    %% Prepare some variables for the window analysis, do this separately because has different number of rows in table
    
    windowSize = -RTwinsize:RTwinsize;
    window2extract = -3:2;
    windowIdx = ismember(windowSize, window2extract);
    
    nwin = length(window2extract);
    plevelWin_Subject   = reshape(repmat((1:nSub)', 1, nwin), nSub*nwin, 1);
    plevelWin_Win       = reshape(repmat(1:nwin, nSub, 1)   , nSub*nwin, 1);
    
    windowTypes = find(~cellfun('isempty',(grep(types2check,'_window'))));
    for itype = 1:length(windowTypes)
        for ibin = 1:nbin2use
            eval(['plevelWin' num2str(ibin) '_' types2check{windowTypes(itype)} ' = reshape(' types2check{windowTypes(itype)} ' {ibinType}(:,ibin,:,windowIdx)  , nSub*nwin, 1);']);
        end
    end
    
    %% write out to csv
    varNames = {...
        'Subject'           , 'Win'                     , ...
        'RT_win5'           , 'RT_CV_win5'              , ...
        'CPP_slope_win5'    , 'CPPr_amplitude_win5'     , 'CPP_ITPC_win5' , ...
        'LHB_slope_win5'    , 'LHB_amplitude_win5'      , ...
        };
    
    %         plevelRTmin5, plevelRTmin4, plevelRTmin3, plevelRTmin2, plevelRTmin1, plevelRT0, plevelRTplus1, plevelRTplus2, plevelRTplus3, plevelRTplus4, plevelRTplus5,
    MAT = table(...
        plevelWin_Subject                           , plevelWin_Win             , ...
        plevelWin5_RT_window                        , plevelWin5_RTcv_window    , ...
        plevelWin5_CPPr_csd_slope_window            , plevelWin5_CPPr_amplitude_window                  , plevelWin5_CPP_ITPC_window, ...
        plevelWin5_beta_pre_response_slope_window   , plevelWin5_beta_base_response_amplitude_window    , ...
        'VariableNames',varNames);
    
    filename_csv{ibinType} = [paths.pop 'participant_level_side(' num2str(sideInfo) ')_bin(' num2str(nbin2use) ')_' bin2use '_' bintype '_window(' num2str(nwin) ')'  fileExt analyse_CDT_fileExt  ];
    % keyboard
    writetable(MAT,[filename_csv{ibinType} '.csv'])
    
    

end
















