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

% sort and bin data

switch type
    case {'CPP_onset','CPPr_slope','CPPr_csd_slope','CPPr_slope_var', 'CPPr_csd_slope_var'}
        datatype       = 'CPP';
        
    case {'CPP','CPPr','CPPr_csd',...
            'N2c','N2i',...
            'pupil_lp','pupil_bp','pupil_lp_1Hz','pupilr_lp',...
            'beta','beta_response','beta_base_response','beta_baseAT_response',...
            'pupil_GLM_Ramp_yhat', 'pupil_conc_GLM_Ramp_yhat'}
        datatype       = 'erp';
        
    case {'CPP_topo','alpha_preTarget_topo','alpha_asym_topo','N2c_topo','beta_pre_response_topo','beta_base_pre_response_topo'}
        datatype       = 'topo';
        
    case {'N2c_ITPC_bar','N2i_ITPC_bar','N2c_256_ITPC_bar','N2i_256_ITPC_bar',...
            'CPP_ITPC_bar','CPPr_ITPC_bar','CPP_csd_ITPC_bar','CPPr_csd_ITPC_bar'}
        datatype       = 'ITPC_bar'; %
        
    case {'N2c_ITPC_band','N2i_ITPC_band','N2c_256_ITPC_band','N2i_256_ITPC_band',...
            'CPP_ITPC_band','CPPr_ITPC_band','CPP_csd_ITPC_band','CPPr_csd_ITPC_band'}
        datatype       = 'ITPC_band'; %
        
    case {'N2c_ITPC','N2i_ITPC','N2c_256_ITPC','N2i_256_ITPC',...
            'CPP_ITPC','CPPr_ITPC','CPP_csd_ITPC','CPPr_csd_ITPC'}
        datatype       = 'ITPC'; %
        
    case 'none'
        datatype = 'none';
    otherwise
        datatype       = 'bar';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select data to bin/sort %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preTarget_validtrials = 0; % If 1, trials with artifacts between -500 and -100 ms before target onset will be excluded, e.g. for alpha analysis. Default 0, changed below if necessary
clear data2use tt t2test twin_bar
switch type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% RT etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'RT'
        data2use    = allRT;
    case 'RT_CV'
        data2use    = allRT;
    case 'RT_zscore'
        data2use    = allRT_zscore;
    case {'RT_window'}
        data2use    = NaN(size(allTrial_window));        
        for irow = 1:size(allTrial_window,1)
            trIdx       = allTrial_window(irow,:);
            trIdx_get   = trIdx(find(~isnan(trIdx)));
            trIdx_save  = find(~isnan(trIdx));
            data2use(irow,trIdx_save) = allRT_zscore(trIdx_get);
        end       
        data2use(~allTrial_window_valid_RT) = NaN;

        tt2use      = -5:5;
        t2test      = -5:1:5;
        twin_bar    = [-1 -1];
    case 'RTcv_window'        
        data2use    = NaN(size(allTrial_window));
        for irow = 1:size(allTrial_window,1)
            trIdx       = allTrial_window(irow,:);
            trIdx_get   = trIdx(find(~isnan(trIdx)));
            trIdx_save  = find(~isnan(trIdx));
            data2use(irow,trIdx_save) = allRT(trIdx_get);
        end
        data2use(~allTrial_window_valid_RT) = NaN;
                
        tt2use      = -5:5;
        t2test      = -5:1:5;
        twin_bar    = [-1 -1];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% pupil
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'pupil_bp_baseline'}
        data2use    = allPupil_bp_baseline;    
    case {'pupil_lp_baseline'}
        data2use    = allPupil_lp_baseline;

    case 'pupil_lp'
        data2use    = allPupil_lp;
        tt2use      = t_pupil;
        t2test      = -500:50:1500;
        twin_bar    = [500 800];

    case 'pupil_bp'
        data2use    = allPupil_bp;
        tt2use      = t_pupil;
        t2test      = -500:50:1500;
        twin_bar    = [500 800];
        
    case 'pupil_lp_1Hz'
        data2use    = allPupil_lp_1Hz;
        tt2use      = t_pupil;
        t2test      = -500:50:1500;
        twin_bar    = [500 800];

    case {'pupilr_lp'}
        data2use    = allPupilr_lp;
        tt2use      = tr_pupil;
        t2test      = -1000:50:400;
        twin_bar    = [-200 0];
    case {'pupilr_bp'}
        data2use    = allPupilr_bp;
        tt2use      = tr_pupil;
        t2test      = -1000:50:400;
        twin_bar    = [-200 0];
        
    case 'pupil_lp_RT_neg200_200'
        data2use    = allPupil_lp_RT_neg200_200;
    case 'pupil_bp_RT_neg200_200'
        data2use    = allPupil_bp_RT_neg200_200;
        
    case {'pupil_GLM_Ramp_yhat', 'pupil_conc_GLM_Ramp_yhat'}
        
        switch type
            case 'pupil_GLM_Ramp_yhat'
                tmp    = allGLM_Ramp_yhat;
            case 'pupil_conc_GLM_Ramp_yhat'
                tmp    = allGLM_conc_Ramp_yhat;
        end
        t2test      = -500:50:1500;

        minLength = min(cellfun(@(X) (length(X)), tmp, 'UniformOutput', true));

        tt2use = t_pupil(t_pupil>=pupilSet.prestim_range(1));
        tt2use = tt2use(1:minLength);

        % put it in same size matrix as other erp variables
        data2use    = NaN(minLength, length(allSubject));
        for isub = 1:nSub
            trIdx = allSubject==isub;
            data2use(:, trIdx) = tmp{isub}(1:minLength, :);
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CPP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    case {'CPP','CPP_onset','CPP_csd','CPP_csd_onset'}
        % stim locked CPP
        switch type
            case {'CPP','CPP_onset'}
                data2use    = allCPP;
            case {'CPP_csd','CPP_csd_onset'}
                data2use    = allCPP_csd;
        end
        tt2use      = t;
        t2test     = -100:10:900;
        twin_bar    = [300 400];
    case 'CPP_topo'
        data2use    = allCPP_topo;

    case {'CPPr','CPPr_amplitude','CPPr_slope','CPPr_slope_var','CPPr_csd','CPPr_csd_amplitude','CPPr_csd_slope','CPPr_csd_slope_var','CPPr_8Hz_slope_var','CPPr_csd_8Hz_slope_var'}
        % response locked CPP and CPP_csd
        switch type
            case {'CPPr','CPPr_amplitude','CPPr_slope','CPPr_slope_var'}
                % non-csd
                data2use    = allCPPr;
            case {'CPPr_8Hz_slope_var'}
                data2use    = allCPPr_8Hz;
            case {'CPPr_csd','CPPr_csd_amplitude','CPPr_csd_slope','CPPr_csd_slope_var'}
                % csd
                data2use    = allCPPr_csd;
            case {'CPPr_csd_8Hz_slope_var'}
                data2use    = allCPPr_csd_8Hz;
        end
        tt2use      = tr;
        t2test      = -400:10:100;
        switch type
            case {'CPPr','CPPr_amplitude','CPPr_csd','CPPr_csd_amplitude'}
                % CPPr amplitude
                twin_bar    = [-50 50];
            case {'CPPr_slope','CPPr_slope_var','CPPr_csd_slope','CPPr_csd_slope_var','CPPr_8Hz_slope_var','CPPr_csd_8Hz_slope_var'}
                % CPP build-up rate was defined as the slope of a straight line fitted to the response-locked waveform (O’Connell et al., 2012; Kelly and O’Connell, 2013; Loughnane et al., 2016) with the time window defined individually for each participant as the 100ms prior to the maximum CPP amplitude pre-response.
                twin_bar    = [-250 -50];
        end
        
    case {'CPP_ITPC','CPP_ITPC_bar','CPP_ITPC_band','CPP_csd_ITPC','CPP_csd_ITPC_bar','CPP_csd_ITPC_band'}
        switch type
            case {'CPP_ITPC','CPP_ITPC_bar','CPP_ITPC_band'}
                data2use    = allCPP_phase;
                tt2use      = (allSPG_times*1000) + t(1);
                ff2use      = allSPG_freq;
            case {'CPP_csd_ITPC','CPP_csd_ITPC_bar','CPP_csd_ITPC_band'}
                data2use    = allCPP_csd_phase;
                tt2use      = (allSPG_times*1000) + t(1);
                ff2use      = allSPG_freq;
        end
        switch datatype
            case {'ITPC','ITPC_bar'}
                switch dataset
                    case 'bigDots'
                        t2test      = [300 550];
                    case 'CD'
                        t2test      = [750 1200];
                end
                
            case 'ITPC_band'
                t2test      = [-200 1000];
        end
        f2test      = [0 4];
        
    case {'CPPr_ITPC','CPPr_ITPC_bar','CPPr_ITPC_band','CPPr_csd_ITPC','CPPr_csd_ITPC_bar','CPPr_csd_ITPC_band'}
        switch type
            case {'CPPr_ITPC','CPPr_ITPC_bar','CPPr_ITPC_band'}
                data2use    = allCPPr_phase;
                tt2use      = (allSPG_timesr*1000) + tr(1);
                ff2use      = allSPG_freq;
            case {'CPPr_csd_ITPC','CPPr_csd_ITPC_bar','CPPr_csd_ITPC_band'}
                data2use    = allCPPr_csd_phase;
                tt2use      = (allSPG_timesr*1000) + tr(1);
                ff2use      = allSPG_freq;
        end
        switch datatype
            case {'ITPC','ITPC_bar'}
                switch dataset
                    case 'bigDots'
                        t2test      = [-300 -50];
                    case 'CD'
                        t2test      = [-300 -50];
                end
                
            case 'ITPC_band'
                t2test      = [-200 1000];
        end
        f2test      = [0 4];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% N2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'N2c'
        data2use    = allN2c;
        tt2use      = t;
        t2test      = -100:10:500;
        twin_bar    = [200 300];
    case 'N2i'
        data2use    = allN2i;
        tt2use      = t;
        t2test      = -100:10:500;
        twin_bar    = [300 400];   
    case 'N2c_amplitude'
        % N2-amplitude was measured as the mean amplitude inside a 100ms window centred on the stimulus-locked grand average peak (N2c: 266ms; N2i: 340ms, Loughnane et al., 2016)
        data2use    = allN2c;
        tt2use      = t;
        t2test      = 266 + [-50 50];
    case 'N2c_latency'
        % N2-latency was then identified as the time-point with the most negative amplitude value in the stimulus-locked waveform between 150-400ms for the N2c and 200-450ms for N2i
        data2use    = allN2c;
        tt2use      = t;
        t2test      = [150 400];

    case {'N2c_ITPC','N2c_ITPC_bar','N2c_ITPC_band'}
        data2use    = allN2c_phase;
        tt2use      = (allSPG_times*1000) + t(1);
        ff2use      = allSPG_freq;
        switch datatype
            case 'ITPC_bar'
                t2test      = [200 400];
            case 'ITPC_band'
                t2test      = [-200 1000];
        end
        f2test      = [0 4];

    case {'N2c_256_ITPC','N2c_256_ITPC_bar','N2c_256_ITPC_band'}
        data2use    = allN2c_256_phase;
        tt2use      = (allSPG_times256*1000) + t(1);
        ff2use      = allSPG_freq256;
        switch datatype
            case 'ITPC_bar'
                t2test      = [200 400];
            case 'ITPC_band'
                t2test      = [-200 1000];
        end
        f2test      = [0 4];        
    case 'N2i_amplitude'
        data2use    = allN2i;
        tt2use      = t;
        t2test      = 340 + [-50 50];
    case 'N2i_latency'
        data2use    = allN2i;
        tt2use      = t;
        t2test      = [200 450];
    case {'N2i_ITPC','N2i_ITPC_bar','N2i_ITPC_band'}
        data2use    = allN2i_phase;
        tt2use      = (allSPG_times*1000) + t(1);
        ff2use      = allSPG_freq;
        switch datatype
            case 'ITPC_bar'
                t2test      = [250 450];
            case 'ITPC_band'
                t2test      = [-200 1000];
        end
        f2test      = [0 4];
    case {'N2i_256_ITPC','N2i_256_ITPC_bar','N2i_256_ITPC_band'}
        data2use    = allN2i_256_phase;
        tt2use      = (allSPG_times256*1000) + t(1);
        ff2use      = allSPG_freq256;
        switch datatype
            case 'ITPC_bar'
                t2test      = [200 400];
            case 'ITPC_band'
                t2test      = [-200 1000];
        end
        f2test      = [0 4];    
    case 'N2c_topo'
        data2use    = allN2c_topo;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% alpha
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'alpha'
        preTarget_validtrials = 1;
        data2use    = allAlpha_preTarget;
    case 'alpha_preTarget_topo'
        preTarget_validtrials = 1;
        data2use    = allAlpha_preTarget_topo;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% beta
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case {'beta','beta_pretarget_amplitude','beta_baseAT_pretarget_amplitude'}
        data2use    = allBeta_postTarget;
        tt2use      = allStft_times;
        t2test      = -500:50:700;
        switch type
            case {'beta_pretarget_amplitude','beta_baseAT_pretarget_amplitude'}
                twin_bar    = [-100 0];
            case 'beta'
                twin_bar    = [400 600];
        end
    case {'beta_response','beta_response_amplitude',...
            'beta_baseAT_response','beta_baseAT_response_amplitude',...
            'beta_base_response','beta_base_response_amplitude',...
            'beta_pre_response_slope'}
        switch type
            case {'beta_response','beta_baseAT_response','beta_response_amplitude','beta_baseAT_response_amplitude','beta_pre_response_slope', 'beta_baseAT_pre_response_slope'}
                data2use    = allBeta_preResponse;
            case {'beta_base_response','beta_base_response_amplitude'}
                data2use    = allBeta_base_preResponse;
        end
        tt2use      = allStft_timesr;
        t2test      = -500:10:100;
        
        % different time window for slope
        switch type
            case {'beta_pre_response_slope'}
                t2test      = [-300 0];
            otherwise
                twin_bar    = [-130 -70];
        end
        
    case 'beta_pre_response_topo'
        data2use    = allBeta_preResponse_topo;
    case 'beta_base_pre_response_topo'
        data2use    = allBeta_base_preResponse_topo;
end

try
    disp(['sorting ' type ' data accoriding to ' bin2use ' in ' num2str(nbin2use) ' bins, for ' num2str(sideInfo) ' stimulus sides'])
catch
    disp(['sorting ' type ' data accoriding to ' bin2use ' in ' num2str(nbin2use) ' bins'])
end

switch dataset
    case 'CD'
        iti2change = [4 7 10];
        for  iti = 1:3
            allItitr(allItitr == iti2change(iti)) = iti;
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop over subjects and bin data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear tIdx trIdx binningVariable binningData binning2use binnedData asym2plot img2plot
%%
for isub = 1:nSub
    itmpsub         = single_participants(isub);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Define the binning trial index %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % note, this is overwritten below when sorting by e.g. alpha
    if (preTarget_validtrials == 1)
        validtr2use = allvalidtr_neg500_RT_200;
    else
        validtr2use = allvalidtr_neg100_RT_200;
    end
    
    trIdx = (allSubject==isub) & validtr2use;
    
    switch bin2use
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Baseline pupil 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        case 'pupil_bp_baseline_regress_iti_side'
            
            designM = [ones(size(allPupil_bp_baseline(trIdx))) allItitr(trIdx) allSideStimtr(trIdx)]; %
            [b, ~, resid] = regress(allPupil_bp_baseline(trIdx), designM);
            binningVariable = resid;
            
            if 0 
                % Different way of regression, where we specifically
                % identify categorical variables, leads to exactly the same
                % result as using the function regress.m 
                LM = fitlm([allItitr(trIdx) allSideStimtr(trIdx)], allPupil_bp_baseline(trIdx), 'y ~ x1+x2', 'CategoricalVars', [1 2]);
                [resid LM.Residuals.Raw]
            end

        case 'pupil_lp_baseline_regress_iti_side'
            
            designM = [ones(size(allPupil_lp_baseline(trIdx))) allItitr(trIdx) allSideStimtr(trIdx)]; %
            [b, ~, resid] = regress(allPupil_lp_baseline(trIdx), designM);
            binningVariable = resid;
            
        case 'pupil_lp_baselineSlope_regress_iti_side'
            
            validtr2use = allvalidtr_neg500_RT_200;
            
            tIdx = find(t>=-500 & ts<=0);
            
            tmp_pupil = allPupil_lp(tIdx, trIdx);
            pupilSlope = NaN(length(find(trIdx)), 1);
            for itrial = 1:length(find(trIdx))
                P = polyfit(1:length(tIdx), squeeze(tmp_pupil(:, itrial))', 1);
                pupilSlope(itrial) = P(1);
            end
            
            designM = [ones(size(pupilSlope)) allItitr(trIdx) allSideStimtr(trIdx)];
            [~, ~, resid] = regress(pupilSlope, designM);
            binningVariable = resid;
            
        case 'pupil_lp_baselineDiff_regress_iti_side'
            validtr2use = allvalidtr_neg500_RT_200;
            
            tIdx = find(t>=-500 & ts<=0);
            
            tmp_pupil = allPupil_lp(tIdx, trIdx);
            pupilDiff = mean(diff(tmp_pupil,1,1))';
            
            designM = [ones(size(pupilDiff)) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(pupilDiff, designM);
            binningVariable = resid;
            
        case 'pupil_lp_1Hz_baseline_regress_iti_side'
            
            designM = [ones(size(allPupil_lp_1Hz_baseline(trIdx))) allItitr(trIdx) allSideStimtr(trIdx)]; %
            [b, ~, resid] = regress(allPupil_lp_1Hz_baseline(trIdx), designM);
            binningVariable = resid;
            
        case {'pupil_bp_baselinePhase_regress_iti_side','pupil_bp_baselinePhase_regress_iti_side_fix'}
            
            designM = [ones(size(allPupil_bp_baselinePhase(trIdx))) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allPupil_bp_baselinePhase(trIdx), designM);
            binningVariable = resid;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% pupil, average around RT
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 'pupil_bp_RT_neg200_200_regress_bl_iti_side'

            designM = [ones(size(allPupil_bp_RT_neg200_200(trIdx))) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allPupil_bp_RT_neg200_200(trIdx), designM);
            binningVariable = resid;
            
        case 'pupil_bp_RT_neg200_200_regress_bl_iti_side_RT'

            tmpRT = allRT(trIdx)/norm(allRT(trIdx));
            designM = [ones(size(allPupil_bp_RT_neg200_200(trIdx))) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx) tmpRT];
            [b, ~, resid] = regress(allPupil_bp_RT_neg200_200(trIdx), designM);
            binningVariable = resid;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% pupil, slope/linear projection
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         case 'pupil_bp_average_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side'
           
             % around maximum of derivative of pupil IRF
             ts_maxDiff_pupilIRF = 318; % time index of maximum of derivative of pupil_IRF
             tIdx = find(ts>=ts_maxDiff_pupilIRF-100 & ts<=ts_maxDiff_pupilIRF+100);
             
             pupilMean = mean(allPupil_bp(tIdx, trIdx), 1)';
             
             designM = [ones(size(pupilMean)) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
             [b, ~, resid] = regress(pupilMean, designM);
             binningVariable = resid;

        case 'pupil_bp_slope_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side'
           
             % around maximum of derivative of pupil IRF
             ts_maxDiff_pupilIRF = 318; % time index of maximum of derivative of pupil_IRF
             tIdx = find(ts>=ts_maxDiff_pupilIRF-100 & ts<=ts_maxDiff_pupilIRF+100);
             
             pupilSlope = NaN(length(find(trIdx)), 1);
             tmp_pupil = allPupil_bp(tIdx, trIdx);
             for itrial = 1:length(find(trIdx))
                 P = polyfit(1:length(tIdx), squeeze(tmp_pupil(:, itrial))', 1);
                 pupilSlope(itrial) = P(1);
             end
             
             designM = [ones(size(pupilSlope)) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
             [b, ~, resid] = regress(pupilSlope, designM);
             binningVariable = resid;
        case 'pupil_bp_slope_maxDiff_pupilIRF_neg250_250_regress_bl_iti_side'
            
            % around maximum of derivative of pupil IRF
            ts_maxDiff_pupilIRF = 318; % time index of maximum of derivative of pupil_IRF
            tIdx = find(ts>=ts_maxDiff_pupilIRF-125 & ts<=ts_maxDiff_pupilIRF+125);
            
            pupilSlope = NaN(length(find(trIdx)), 1);
            tmp_pupil = allPupil_bp(tIdx, trIdx);
            for itrial = 1:length(find(trIdx))
                P = polyfit(1:length(tIdx), squeeze(tmp_pupil(:, itrial))', 1);
                pupilSlope(itrial) = P(1);
            end
            
            designM = [ones(size(pupilSlope)) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
             [b, ~, resid] = regress(pupilSlope, designM);
             binningVariable = resid;
             
             
        case 'pupil_bp_slope_maxDiff_pupilIRF_300_200_regress_bl_iti_side'
           
             % around maximum of derivative of pupil IRF
             ts_maxDiff_pupilIRF = 318; % time index of maximum of derivative of pupil_IRF
             tIdx = find(ts>=150 & ts<=ts_maxDiff_pupilIRF+100);
             
             pupilSlope = NaN(length(find(trIdx)), 1);
             tmp_pupil = allPupil_bp(tIdx, trIdx);
             for itrial = 1:length(find(trIdx))
                 P = polyfit(1:length(tIdx), squeeze(tmp_pupil(:, itrial))', 1);
                 pupilSlope(itrial) = P(1);
             end
             
             designM = [ones(size(pupilSlope)) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
             [b, ~, resid] = regress(pupilSlope, designM);
             binningVariable = resid;
             
        case 'pupil_bp_diff_maxDiff_pupilIRF_200_200_regress_bl_iti_side'
           
             % around maximum of derivative of pupil IRF
             ts_maxDiff_pupilIRF = 318; % time index of maximum of derivative of pupil_IRF
             tIdx = find(ts>=100 & ts<=ts_maxDiff_pupilIRF+101);
             
             tmp_pupil = allPupil_bp(tIdx, trIdx);
             pupilDiff = mean(diff(tmp_pupil,1,1))';
             
             designM = [ones(size(pupilDiff)) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
             [b, ~, resid] = regress(pupilDiff, designM);
             binningVariable = resid;
             
             
         case 'pupil_bp_linearProjection_maxDiff_pupilIRF_neg200_200_regress_bl_iti_side'
           
             % around maximum of derivative of pupil IRF
             ts_maxDiff_pupilIRF = 318; % time index of maximum of derivative of pupil_IRF
             tIdx = find(ts>=ts_maxDiff_pupilIRF-100 & ts<=ts_maxDiff_pupilIRF+100);
             
             pupilLinProj = linearProjection(allPupil_bp(tIdx, trIdx), 1:length(tIdx));
             designM = [ones(size(pupilLinProj)) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
             [b, ~, resid] = regress(pupilLinProj, designM);
             binningVariable = resid;
             
         case 'pupil_bp_linearProjection_maxDiff_pupilIRF_0_200_regress_bl_iti_side'
           
             % around maximum of derivative of pupil IRF
             ts_maxDiff_pupilIRF = 318; % time index of maximum of derivative of pupil_IRF
             tIdx = find(ts>=0 & ts<=ts_maxDiff_pupilIRF+100);
             
             pupilLinProj = linearProjection(allPupil_bp(tIdx, trIdx), t(tIdx));
             designM = [ones(size(pupilLinProj)) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
             [b, ~, resid] = regress(pupilLinProj, designM);
             binningVariable = resid;
             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% GLM, pupil
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 'GLM_pupil_StimResp_stim_regress_bl_iti_side'
            designM = [ones(size(allGLM_pupil_StimResp_stim_beta(trIdx))) allPupil_bp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allGLM_pupil_StimResp_stim_beta(trIdx), designM);
            binningVariable = resid;
            
            
        case 'GLM_pupil_Ramp_stim_regress_iti_side'
            designM = [ones(size(allGLM_pupil_Ramp_stim_beta(trIdx))) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allGLM_pupil_Ramp_stim_beta(trIdx), designM);
            binningVariable = resid;
            
        case 'GLM_pupil_Ramp_stim_regress_bl_iti_side'
            designM = [ones(size(allGLM_pupil_Ramp_stim_beta(trIdx))) allPupil_bp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allGLM_pupil_Ramp_stim_beta(trIdx), designM);
            binningVariable = resid;
            
        case 'GLM_pupil_Ramp_stim_regress_bl_blPhase_iti_side'
            designM = [ones(size(allGLM_pupil_Ramp_stim_beta(trIdx))) allPupil_bp_baseline(trIdx) cos(allPupil_bp_baselinePhase(trIdx)) sin(allPupil_bp_baselinePhase(trIdx)) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allGLM_pupil_Ramp_stim_beta(trIdx), designM);
            binningVariable = resid;
            
        case 'GLM_pupil_Ramp_stim_regress_blPhase_iti_side'
            % use the cosine and sine of the phase as independent
            % orthogonal predictors of pupil response
            designM = [ones(size(allGLM_pupil_Ramp_stim_beta(trIdx))) cos(allPupil_bp_baselinePhase(trIdx)) sin(allPupil_bp_baselinePhase(trIdx)) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allGLM_pupil_Ramp_stim_beta(trIdx), designM);
            binningVariable = resid;
            
            if 0
                subplotgap = [0.06 0.08];

                % confirm there is no correlation between phase and
                % response after regression
                [rho1 pval1 ] = circ_corrcl(allPupil_bp_baselinePhase(trIdx), allGLM_pupil_Ramp_stim_beta(trIdx));
                [rho2 pval2 ] = circ_corrcl(allPupil_bp_baselinePhase(trIdx), binningVariable);
                
                [figHandle, fSet] = figInit('fig',1, {'height', 10; 'width', 22});
                axH = subtightplot(1, 3, 2, subplotgap);
                figInit('ax');
                plot_subplot_label(axH, fSet, 'B')
                
                hold on
                x1 = linspace(-pi,pi);
                p = polyfit(allPupil_bp_baselinePhase(trIdx), allGLM_pupil_Ramp_stim_beta(trIdx),3);
                y1 = polyval(p,x1);
                plot(allPupil_bp_baselinePhase(trIdx), allGLM_pupil_Ramp_stim_beta(trIdx),'.', 'color', fSet.colors(1,:),'markersize',fSet.MarkerSize)
                plot(x1, y1,'color', fSet.colors(1,:),'linewidth',fSet.LineWidth)
                
                p = polyfit(allPupil_bp_baselinePhase(trIdx), binningVariable,3);
                y1 = polyval(p,x1);
                plot(allPupil_bp_baselinePhase(trIdx), binningVariable,'.', 'color', fSet.colors(2,:),'markersize',fSet.MarkerSize)
                plot(x1, y1,'color', fSet.colors(2,:),'linewidth',fSet.LineWidth)
                
                tmp = get(gca,'ylim');
                
                text(-0.5, tmp(1) * 0.78, ['r = ' num2str(rho1, '%1.2f') ', p = ' num2str(pval1, '%1.2e\n')],'color',fSet.colors(1,:),'fontsize',fSet.Fontsize_text_in)
                text(-0.5, tmp(1) * 0.9, ['r = ' num2str(rho2, '%1.2f') ', p = ' num2str(pval2, '%1.2f')],'color',fSet.colors(2,:),'fontsize',fSet.Fontsize_text_in)
                
                xlim([-pi*1.1 pi*1.1])
                %                 ylim([min(allPupil_lp_resp_locked_neg500_100(trIdx)) * 1.1 max(allPupil_lp_resp_locked_neg500_100(trIdx))*1.1])
                
                xlabel('Pre-target pupil phase (Rad)','fontsize',fSet.Fontsize_text)
                ylabel('Pupil response amplitude','fontsize',fSet.Fontsize_text)
                
                axis square
                saveFigName = 'circularRegression';
                saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.svg'], 'svg');
                saveas(gcf, [paths.pop 'fig' filesep 'manuscript' filesep saveFigName '.png'], 'png');
                
            end
            
        case 'GLM_pupil_Ramp_sust_regress_iti_side'
            designM = [ones(size(allGLM_pupil_Ramp_sust_beta(trIdx))) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allGLM_pupil_Ramp_sust_beta(trIdx), designM);
            binningVariable = resid;

        case 'GLM_pupil_Ramp_sust_regress_bl_iti_side'
            designM = [ones(size(allGLM_pupil_Ramp_sust_beta(trIdx))) allPupil_bp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allGLM_pupil_Ramp_sust_beta(trIdx), designM);
            binningVariable = resid;
            
            
        case 'GLM_pupil_Ramp_stim_VIF5_regress_iti_side'
            
            trIdx = trIdx & (allGLM_pupil_Ramp_VIF(:,1) < 5);
            
            designM = [ones(size(allGLM_pupil_Ramp_stim_beta(trIdx))) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allGLM_pupil_Ramp_stim_beta(trIdx), designM);
            binningVariable = resid;
            
        case 'GLM_pupil_Ramp_stim_VIF5_regress_bl_iti_side'
            
            trIdx = trIdx & (allGLM_pupil_Ramp_VIF(:,1) < 5);
            
            designM = [ones(size(allGLM_pupil_Ramp_stim_beta(trIdx))) allPupil_bp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allGLM_pupil_Ramp_stim_beta(trIdx), designM);
            binningVariable = resid;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Other
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'RT_regress_iti_side'
            designM = [ones(size(allRT(trIdx))) allPupil_lp_baseline(trIdx) allItitr(trIdx) allSideStimtr(trIdx)];
            [b, ~, resid] = regress(allRT(trIdx), designM);
            binningVariable = resid;
            
        case 'pretarget_alpha'
            validtr2use = allvalidtr_neg500_0;
            trIdx = allSubject==isub & validtr2use;
            binningVariable = allAlpha_preTarget(trIdx,1);

        case 'none'
            binningVariable = (allSubject==isub & validtr2use );
            
        otherwise
            clear binningVariable
            warning('No binning variable defined')
            
    end
    
                            
    
    %%% bin the data
    binningData = nan(length(allSubject),1);
    binningData(trIdx) = binningVariable;
    
    binning2use = nan(length(allSubject),1);
    
    switch bin2use
        case 'pupil_bp_baselinePhase_regress_iti_side_fix'
            [binEdges, binning2use(trIdx), ~] = binVal(binningData(trIdx,1), nbin2use, 'manual', linspace(-pi,pi,nbin2use+1));
        otherwise
            [binEdges, binning2use(trIdx), ~] = binVal(binningData(trIdx,1), nbin2use, bintype);
    end

    binning2use(isnan(binning2use)) = 0;
    

    %% some processing specific for certain analyses

    switch type
        case {'beta_baseAT_response','beta_baseAT_response_amplitude','beta_baseAT_pretarget_amplitude'}
            % implement across trial (AT) baseline for LHB, per reviewers request
            baseline_beta = mean(mean(allBeta_postTarget(find(allStft_times>=BL_spectrum(1) & allStft_times<=BL_spectrum(2)), trIdx),1),2);
            data2use(:,trIdx) = data2use(:,trIdx)-repmat(baseline_beta,[size(data2use,1),length(find(trIdx))]);
    end
    
    clear trIdx baseline_beta

    %% Get matrix with subject number, trial number, stimulus side and bin
    % number to write to csv file
    for ibin = 1:nbin2use
        for isideStim = 1:sideInfo
            if sideInfo > 1
                trIdx = (allSideStimtr==isideStim);
                plevelSub(isub, ibin, isideStim) = isub;
                plevelnTrial(isub, ibin, isideStim) = length(find((binning2use==ibin) & trIdx));
                plevelSide(isub, ibin, isideStim) = isideStim;
                plevelBin(isub, ibin, isideStim) = ibin;
            else
                plevelSub(isub, ibin, isideStim) = isub;
                plevelnTrial(isub, ibin, isideStim) = length(find((binning2use==ibin)));
                plevelSide(isub, ibin, isideStim) = isideStim;
                plevelBin(isub, ibin, isideStim) = ibin;
            end
        end
    end
    if strcmpi(type, 'none') && strcmpi(datatype, 'none') && (isub == nSub)
        return
    end
    clear trIdx

    %%
    %%%%%%%%%%%%%%%%%%%
    %%% Bin the data %%
    %%%%%%%%%%%%%%%%%%%
    switch datatype  
        %%
        case 'erp'
            
            if isub == 1
                binnedData = NaN(nSub, nbin2use, sideInfo, size(data2use,1));
            end
            switch type
                case 'RTcv_window'                   
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                binnedData(isub, ibin, isideStim, :) = squeeze(nanstd(data2use(:,(binning2use==ibin)& trIdx),[],2)) ./ squeeze(nanmean(data2use(:,(binning2use==ibin)& trIdx),2));
                            else
                                binnedData(isub, ibin, isideStim, :) = squeeze(nanstd(data2use(:,(binning2use==ibin)),[],2)) ./ squeeze(nanmean(data2use(:,(binning2use==ibin)),2));
                            end
                        end
                    end
                    
                otherwise
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                binnedData(isub, ibin, isideStim, :) = squeeze(mean(data2use(:,(binning2use==ibin)& trIdx),2));
                            else
                                binnedData(isub, ibin, isideStim, :) = squeeze(nanmean(data2use(:,(binning2use==ibin)),2));
                            end
                        end
                    end
            end
            
        case 'bar'
            %%
            if isub == 1
                binnedData = NaN(nSub, nbin2use, sideInfo);
            end
            switch type
                case {'RT_CV'}
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                binnedData(isub, ibin, isideStim) = squeeze(std(data2use((binning2use==ibin) & trIdx))) / squeeze(mean(data2use((binning2use==ibin) & trIdx)));
                            else
                                binnedData(isub, ibin, isideStim) = squeeze(std(data2use(binning2use==ibin))) / squeeze(mean(data2use(binning2use==ibin)));
                            end
                        end
                    end
                case {'beta_pre_response_slope'}

                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (binning2use==ibin) & (allSideStimtr==isideStim);
                            else
                                trIdx = (binning2use==ibin);
                            end
                            
                            [~,slope_timeframe_index(1)] = min(abs(tt2use - t2test(1)));%
                            [~,slope_timeframe_index(2)] = min(abs(tt2use - t2test(2)));%
                            
                            tmp = squeeze(mean(data2use(:,trIdx),2));
                            
                            coef = polyfit(tt2use(slope_timeframe_index(1):slope_timeframe_index(2)),squeeze(tmp(slope_timeframe_index(1):slope_timeframe_index(2)))',1);% coef gives 2 coefficients fitting r = slope * x + intercept
                            binnedData(isub, ibin, isideStim)=coef(1);
                            
                        end
                    end

                    
                case {'N2c_latency','N2i_latency'}
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            tIdx    = find(tt2use>t2test(1) & tt2use<t2test(2));
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                tmpdat  = squeeze(mean(data2use(tIdx,(binning2use==ibin) & trIdx),2));
                            else
                                tmpdat  = squeeze(mean(data2use(tIdx,(binning2use==ibin)),2));
                            end
                            [minN2, minTime]   = min(tmpdat);
                            binnedData(isub, ibin, isideStim) = tt2use(tIdx(minTime));

                        end
                    end
                                        
                case {'N2c_amplitude','N2i_amplitude'}
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            
                            tIdx    = find(tt2use>t2test(1) & tt2use<t2test(2));
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                binnedData(isub, ibin, isideStim)  = squeeze(mean(mean(data2use(tIdx,(binning2use==ibin)& trIdx),1),2));
                            else
                                binnedData(isub, ibin, isideStim)  = squeeze(mean(mean(data2use(tIdx,(binning2use==ibin)),1),2));
                            end
                        end
                    end
                    
                case {'N2c','N2i',...
                        'CPP','CPPr_csd','CPPr_amplitude','CPPr_csd_amplitude', ...
                        'beta_response_amplitude', 'beta_base_response_amplitude', 'beta_baseAT_response_amplitude', 'beta_pretarget_amplitude', 'beta_baseAT_pretarget_amplitude', ...
                        'pupil_lp','pupilr_lp','pupil_bp','pupilr_bp'}
                    % get average data in time window from time series data
                    tIdx = (tt2use>=twin_bar(1) & tt2use<=twin_bar(2));
                    
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                binnedData(isub, ibin, isideStim) = squeeze(mean(mean(data2use(tIdx,(binning2use==ibin) & trIdx),1),2));
                            else
                                binnedData(isub, ibin, isideStim) = squeeze(nanmean(nanmean(data2use(tIdx,binning2use==ibin),1),2));
                            end
                        end
                    end                 
                    
                case 'CPPr_amplitude_var'
                    % get average data in time window from time series data
                    tIdx = (tt2use>=twin_bar(1) & tt2use<=twin_bar(2));
                    
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                binnedData(isub, ibin, isideStim) = var(squeeze(mean(data2use(tIdx,(binning2use==ibin) & trIdx),1)));
                            else
                                binnedData(isub, ibin, isideStim) = var(squeeze(mean(data2use(tIdx,binning2use==ibin),1)));
                            end
                        end
                    end
                    
                otherwise
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                binnedData(isub, ibin, isideStim) = squeeze(mean(data2use((binning2use==ibin) & trIdx),1));
                            else
                                binnedData(isub, ibin, isideStim) = squeeze(nanmean(data2use(binning2use==ibin),1));
                            end
                        end
                    end
                    
            end
            
            %% TOPO
        case 'topo'
            if isub == 1
                binnedData = NaN(nSub, nbin2use, size(data2use,1));
            end
            switch type
                case {'N2c_topo'}
                    if isub == 1
                        binnedData = NaN(nSub, nbin2use, 2, size(data2use,1));
                    end
                    for ibin = 1:nbin2use
                        trIdx1 = (allSideStimtr==1);
                        trIdx2 = (allSideStimtr==2);
                        binnedData(isub, ibin, 1, :) = squeeze(mean(data2use(:,(binning2use==ibin) & trIdx1),2));
                        binnedData(isub, ibin, 2, :) = squeeze(mean(data2use(:,(binning2use==ibin) & trIdx2),2));
                    end
                otherwise
                    for ibin = 1:nbin2use
                        binnedData(isub, ibin, :) = squeeze(mean(data2use(:,binning2use==ibin),2));
                    end
            end
            
            %% CPP
        case 'CPP'
            switch type
                case 'CPP_onset'
                    if isub == 1
                        binnedData = NaN(nSub, nbin2use, sideInfo);
                    end
                    tmpbinnedData = cell(nbin2use, sideInfo);
                    
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                tmpbinnedData{ibin, isideStim} = squeeze(data2use(:,(binning2use==ibin)& trIdx));
                            else
                                tmpbinnedData{ibin, isideStim} = squeeze(data2use(:,(binning2use==ibin)));
                            end
                            
                            % get specific CPP settings
                            % Because of binning procedure, we don't find
                            % an onset latency for each subject/bin
                            switch dataset
                                case 'bigDots'
                                    CPPonset_settings_bigDots
                                case 'CD'
                                    % CPP_settings_CD
                                    win_mean_change = 0;
                                    consecutive_windows=15;%15 works well for everybody else
                            end
                            
                            if plotCPP
                                disp([num2str(isub) ' ' subject_folder{itmpsub}])
                            end
                            
                            [binnedData(isub,ibin,isideStim), starttime] = getOnsetCPP(tmpbinnedData{ibin,isideStim}, subject_folder{itmpsub}, t, CPP_search_t, max_search_window, win_mean_change, consecutive_windows, plotCPP,'CPP');
                            
                        end
                    end
                case {'CPPr_slope', 'CPPr_csd_slope'}
                                        
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (binning2use==ibin) & (allSideStimtr==isideStim);
                            else
                                trIdx = (binning2use==ibin);
                            end
                            
                            
                            if isub == 1 && ibin == 1 && isideStim == 1
                                binnedData = NaN(nSub, nbin2use, sideInfo);
                            end
                            
                            tmpbinnedData               = squeeze(mean(data2use(:, trIdx),2));
                            [binnedData(isub,ibin,:)]   = getSlopeCPP(tmpbinnedData, sideInfo, subject_folder{itmpsub}, tr,side_tags, plotCPP, twin_bar,2);
                        end                        
                    end
                    

                case {'CPPr_slope_var','CPPr_csd_slope_var'}
                    
                    if isub == 1
                        binnedData = NaN(nSub, nbin2use, sideInfo);
                    end
                    
                    allSlopes = NaN(length(binning2use),1);
                    trIdx = find(binning2use);
                    for itrial = 1:length(trIdx)
                        allSlopes(trIdx(itrial)) = getSlopeCPP(squeeze(data2use(:,trIdx(itrial))), 1, subject_folder{itmpsub}, tr, side_tags, plotCPP, twin_bar,2);
                    end
                        
                    for ibin = 1:nbin2use
                        for isideStim = 1:sideInfo
                            if sideInfo > 1
                                trIdx = (allSideStimtr==isideStim);
                                binnedData(isub,ibin,isideStim) = var(allSlopes(binning2use==ibin)& trIdx);
                            else
                                binnedData(isub,ibin,isideStim) = var(allSlopes(binning2use==ibin));
                            end
                        end                        
                    end
            end
            
        case {'ITPC'}
            if isub == 1
                binnedData = NaN(nSub, nbin2use, sideInfo, size(data2use,1), size(data2use,2));
            end
            
            for ibin = 1:nbin2use
                for isideStim = 1:sideInfo
                    if sideInfo > 1
                        trIdx = (binning2use==ibin) & (allSideStimtr==isideStim);
                    else  
                        trIdx = (binning2use==ibin);
                    end

                    binnedData(isub,ibin,isideStim,:,:) = squeeze(abs(mean(exp(1i*(data2use(:, :, trIdx))),3)));
                    
                    if ~isempty(find(isnan(binnedData(isub,ibin,isideStim,:,:))))
                        keyboard
                    end
                end
            end
        case {'ITPC_band'}
            if isub == 1
                binnedData = NaN(nSub, nbin2use, sideInfo, size(data2use,1));
            end
            
            clear tmp
            for ibin = 1:nbin2use
                for isideStim = 1:sideInfo
                    if sideInfo > 1
                        trIdx = (binning2use==ibin) & (allSideStimtr==isideStim);
                    else  
                        trIdx = (binning2use==ibin);
                    end
                    
                    tmp = squeeze(abs(mean(exp(1i*(data2use(:, :, trIdx))),3)));
                    [~,f2plot(1)] = min(abs(f2test(1) - allSPG_freq));
                    [~,f2plot(2)] = min(abs(f2test(2) - allSPG_freq));
                    
                    binnedData(isub,ibin,isideStim,:) = mean(tmp(:,f2plot(1):f2plot(2)),2);
                    
                    if ~isempty(find(isnan(binnedData(isub,ibin,isideStim,:))))
                        keyboard
                    end

                end
            end
            
        case {'ITPC_bar'}
            if isub == 1
                binnedData = NaN(nSub, nbin2use, sideInfo);
            end
            
            clear tmp
            for ibin = 1:nbin2use
                for isideStim = 1:sideInfo
                    
                    if sideInfo > 1
                        trIdx = (binning2use==ibin) & (allSideStimtr==isideStim);
                    else
                        trIdx = (binning2use==ibin);
                    end
                    
                    tmp = squeeze(abs(mean(exp(1i*(data2use(:, :, trIdx))),3)));
                    [~,f2plot(1)] = min(abs(f2test(1) - allSPG_freq));
                    [~,f2plot(2)] = min(abs(f2test(2) - allSPG_freq));
                    
                    [~,t2plot(1)] = min(abs(t2test(1) - tt2use));
                    [~,t2plot(2)] = min(abs(t2test(2) - tt2use));
                    
                    binnedData(isub,ibin,isideStim) = mean(mean(tmp(t2plot(1):t2plot(2),f2plot(1):f2plot(2)),1),2);
                    
                    if ~isempty(find(isnan(binnedData(isub,ibin,isideStim))))
                        keyboard
                    end
                end
            end
        case 'none'
        otherwise
            error('data type not defined')
    end
    
    
end
