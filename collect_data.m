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

%% append all subject files, collect in to one big matrix

ID=[]; allTrial=[]; allBlock=[]; allBlockTrial=[]; allSideStimtr=[]; allSideStimtr_window=[]; allMotiontr=[]; allItitr=[]; allHits=[]; allMisses=[]; allValidrlock=[]; allSubject=[];  %idx

allvalidtr=[]; allvalidtr_neg100_0=[]; allvalidtr_neg500_0=[];
allvalidtr_neg100_RT_200=[]; allvalidtr_neg500_RT_200=[];

%RT etc
allRT=[]; allRT_log=[]; allRT_zscore=[]; allValidRT=[]; 

%N2
allN2c=[]; allN2i=[]; allN2c_topo=[];
%CPP
allCPP=[]; allCPP_csd=[]; allCPPr=[]; allCPPr_csd=[]; allCPP_topo=[]; %

% alpha
allAlphaRh_preTarget=[]; allAlphaLh_preTarget=[];
allAlpha_preTarget=[]; allAlphaAsym_preTarget=[]; allAlpha_preTarget_topo=[]; allAlphaAsym_preTarget_topo=[];

% beta
allBeta_postTarget=[]; allBeta_base_postTarget=[]; allBeta_preResponse=[]; allBeta_base_preResponse=[]; allBeta_preResponse_topo=[]; allBeta_base_preResponse_topo=[];

% SPG
allN2c_power = []; allN2c_phase = []; allN2i_power = []; allN2i_phase = [];
allN2c_256_power = []; allN2c_256_phase = []; allN2i_256_power = []; allN2i_256_phase = [];
allCPP_power = []; allCPP_phase = [];

%pupil.
allPupil_lp=[]; allPupilr_lp=[];
allPupil_bp=[]; allPupilr_bp=[];
allPupil_lp_baseline=[]; allPupil_bp_baseline=[];
allPupil_lp_RT_neg200_200=[]; allPupil_bp_RT_neg200_200=[];
allMinRT = []; allMeanRT = [];
allPupil_lp_minRT_neg200_200=[]; allPupil_bp_minRT_neg200_200=[];
allPupil_lp_meanRT_neg200_200=[]; allPupil_bp_meanRT_neg200_200=[];
allPupil_lp_linproj_resp_locked_neg500_200=[];allPupil_bp_linproj_resp_locked_neg500_200=[];
allPupil_lp_baseline_zscore = [];
% trial window index
allTrial_window=[]; allTrial_window_valid_RT=[]; allTrial_window_valid_neg100_RT_200=[]; allTrial_window_valid_neg500_0=[];

clear validtrials filenames
tmpsub=0;
for isub = single_participants
    tmpsub=tmpsub+1;
    
    switch dataset
        case 'CD'
            filenames = [paths.s(isub).base allsubj{isub} fileExt analyse_CDT_fileExt '_ST.mat'];
        case 'bigDots'
            filenames = [paths.s(isub).base allsubj{isub} fileExt analyse_CDT_fileExt '_ST.mat'];
    end
    disp(['loading: ' filenames])
    clear DAT
    DAT = load([filenames]);
    
    % idx
    allTrial        = [allTrial         ; DAT.trialIdx];
    allBlock        = [allBlock         ; DAT.blockIdx];
    allBlockTrial   = [allBlockTrial    ; DAT.blockTrialIdx];
    allSideStimtr   = [allSideStimtr    ; DAT.sideStimtr];
    if isfield(DAT,'motiontr')
        allMotiontr   = [allMotiontr        ; DAT.motiontr];
    end
    allItitr        = [allItitr         ; DAT.ititr];
    allHits         = [allHits          ; DAT.hits];
    allMisses       = [allMisses        ; DAT.misses];
    allValidrlock   = [allValidrlock    ; DAT.validrlock];
    allSubject      = [allSubject       ; ones(length(DAT.subject),1)*tmpsub];
    ID              = [ID ; repmat(allsubj(isub), length(DAT.subject),1)];
    
    
    % valid tr ind
    allvalidtr                  = [allvalidtr               ; DAT.validtr.validRT];
    allvalidtr_neg100_0         = [allvalidtr_neg100_0      ; DAT.validtr.neg100_0];
    allvalidtr_neg500_0         = [allvalidtr_neg500_0      ; DAT.validtr.neg500_0];
    allvalidtr_neg100_RT_200    = [allvalidtr_neg100_RT_200 ; DAT.validtr.neg100_RT_200];
    allvalidtr_neg500_RT_200    = [allvalidtr_neg500_RT_200 ; DAT.validtr.neg500_RT_200];
    
    %RT etc
    allRT               = [allRT        ; DAT.subRT];
    allRT_log           = [allRT_log    ; DAT.subRT_log];
    allRT_zscore        = [allRT_zscore ; DAT.subRT_zscore];
    allValidRT          = [allValidRT   ; DAT.validRT];
        
    %     % erp
    %     allERP              = cat(3, allERP,        ERP);
    %     allERPr             = cat(3, allERPr,       ERPr);
    %     allERP_csd          = cat(3, allERP_csd,    ERP_csd);
    %     allERPr_csd         = cat(3, allERPr_csd,   ERPr_csd);
    
    %N2
    allN2c          = [allN2c DAT.N2c];
    allN2i          = [allN2i DAT.N2i];
    allN2c_topo     = [allN2c_topo    DAT.N2c_topo];
    %     allN2i_topo     = [allN2i_topo    DAT.N2i_topo];
    
    %CPP
    allCPP          = [allCPP       DAT.CPP];
    allCPPr         = [allCPPr      DAT.CPPr];
    allCPP_csd      = [allCPP_csd   DAT.CPP_csd];
    allCPPr_csd     = [allCPPr_csd  DAT.CPPr_csd];
    allCPP_topo     = [allCPP_topo  DAT.CPP_topo];
    
    %alpha
    allAlphaRh_preTarget             = [allAlphaRh_preTarget            ; DAT.alphaRh_preTarget];
    allAlphaLh_preTarget             = [allAlphaLh_preTarget            ; DAT.alphaLh_preTarget];
    allAlpha_preTarget               = [allAlpha_preTarget              ; DAT.alpha_preTarget];
    allAlphaAsym_preTarget           = [allAlphaAsym_preTarget          ; DAT.alphaAsym_preTarget];
    allAlpha_preTarget_topo          = [allAlpha_preTarget_topo         DAT.alpha_preTarget_topo];
    allAlphaAsym_preTarget_topo      = [allAlphaAsym_preTarget_topo     DAT.alphaAsym_preTarget_topo];
    
    %beta
    allBeta_postTarget              = [allBeta_postTarget           DAT.beta_postTarget];
    allBeta_base_postTarget         = [allBeta_base_postTarget      DAT.beta_base_postTarget];
    allBeta_preResponse             = [allBeta_preResponse          DAT.beta_preResponse];
    allBeta_base_preResponse        = [allBeta_base_preResponse         DAT.beta_base_preResponse];
    allBeta_preResponse_topo        = [allBeta_preResponse_topo         DAT.beta_preResponse_topo];
    allBeta_base_preResponse_topo   = [allBeta_base_preResponse_topo    DAT.beta_base_preResponse_topo];
    
    allN2c_power        = cat(3,allN2c_power, DAT.SPG.N2c_power);
    allN2c_phase        = cat(3,allN2c_phase, DAT.SPG.N2c_phase);
    allN2c_256_power    = cat(3,allN2c_256_power, DAT.SPG.N2c_256_power);
    allN2c_256_phase    = cat(3,allN2c_256_phase, DAT.SPG.N2c_256_phase);
    allN2i_power        = cat(3,allN2i_power, DAT.SPG.N2i_power);
    allN2i_phase        = cat(3,allN2i_phase, DAT.SPG.N2i_phase);
    allN2i_256_power    = cat(3,allN2i_256_power, DAT.SPG.N2i_256_power);
    allN2i_256_phase    = cat(3,allN2i_256_phase, DAT.SPG.N2i_256_phase);
    allCPP_power        = cat(3,allCPP_power, DAT.SPG.CPP_power);
    allCPP_phase        = cat(3,allCPP_phase, DAT.SPG.CPP_phase);
    
    stft_times      = DAT.spectral_t;
    stft_timesr     = DAT.spectral_tr;
    
    SPG_times       = DAT.SPG.tt;
    SPG_times256    = DAT.SPG.tt_256;
    SPG_freq        = DAT.SPG.ff;
    SPG_freq256     = DAT.SPG.ff_256;
    
    %pupil
    allPupil_bp                = [allPupil_bp  DAT.pupil.bp];
    allPupilr_bp               = [allPupilr_bp DAT.pupilr.bp];
    allPupil_lp                = [allPupil_lp  DAT.pupil.lp];
    allPupilr_lp               = [allPupilr_lp DAT.pupilr.lp];
    
    allPupil_lp_baseline       = [allPupil_lp_baseline ; DAT.pupil.baseline.lp];
    allPupil_bp_baseline       = [allPupil_bp_baseline ; DAT.pupil.baseline.bp];

    allPupil_lp_baseline_zscore       = [allPupil_lp_baseline_zscore ; DAT.pupil.baseline_zscore.lp];

    allPupil_bp_RT_neg200_200   = [allPupil_bp_RT_neg200_200 ; DAT.pupil.RT.bp.neg200_200];
    allPupil_lp_RT_neg200_200   = [allPupil_lp_RT_neg200_200 ; DAT.pupil.RT.lp.neg200_200];

    allMinRT                       = [allMinRT ; DAT.pupil.minRT.windowCentre];
    allPupil_bp_minRT_neg200_200   = [allPupil_bp_minRT_neg200_200 ; DAT.pupil.minRT.bp.neg200_200];
    allPupil_lp_minRT_neg200_200   = [allPupil_lp_minRT_neg200_200 ; DAT.pupil.minRT.lp.neg200_200];

    allMeanRT                       = [allMeanRT ; DAT.pupil.meanRT.windowCentre];
    allPupil_bp_meanRT_neg200_200   = [allPupil_bp_meanRT_neg200_200 ; DAT.pupil.meanRT.bp.neg200_200];
    allPupil_lp_meanRT_neg200_200   = [allPupil_lp_meanRT_neg200_200 ; DAT.pupil.meanRT.lp.neg200_200];

    allPupil_lp_linproj_resp_locked_neg500_200      = [allPupil_lp_linproj_resp_locked_neg500_200     ; DAT.lpPupil.lp.resp_locked_neg500_200];
    allPupil_bp_linproj_resp_locked_neg500_200      = [allPupil_bp_linproj_resp_locked_neg500_200     ; DAT.lpPupil.bp.resp_locked_neg500_200];
    
    if isub == single_participants(1)
        trialadd = 0;
    else
        trialadd = max(allTrial_window(:));
    end
    allTrial_window                         = [allTrial_window                      (DAT.trialWindow_idx' + trialadd)];
    allTrial_window_valid_RT                = [allTrial_window_valid_RT             DAT.trialWindow_valid_RT'];
    allTrial_window_valid_neg100_RT_200     = [allTrial_window_valid_neg100_RT_200  DAT.trialWindow_valid_neg100_RT_200'];
    allTrial_window_valid_neg500_0          = [allTrial_window_valid_neg500_0       DAT.trialWindow_valid_neg500_0'];

    clear DAT
end

%% save
filename_mat = [paths.pop 'allSub_singleTrial' fileExt analyse_CDT_fileExt '.mat'];
disp(['saving ' filename_mat])
save([filename_mat],'-v7.3')

