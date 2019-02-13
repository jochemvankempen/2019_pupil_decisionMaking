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


% concatenated GLM
allGLM_conc_pupil_BIC = [];
allGLM_conc_pupil_StimResp_VIF=[];
allGLM_conc_pupil_Boxc_tstat=[]; allGLM_conc_pupil_Ramp_tstat=[]; allGLM_conc_pupil_Ramp2thresh_tstat=[]; 
allGLM_conc_pupil_Boxc_beta=[]; allGLM_conc_pupil_Ramp_beta=[]; allGLM_conc_pupil_Ramp2thresh_beta=[]; 
allGLM_conc_pupil_Boxc_VIF=[]; allGLM_conc_pupil_Ramp_VIF=[]; allGLM_conc_pupil_Ramp2thresh_VIF=[]; 

% single trial
allGLM_pupil_StimResp_stim_beta=[]; allGLM_pupil_StimResp_resp_beta=[]; allGLM_pupil_StimResp_VIF=[]; allGLM_pupil_StimResp_rsquare=[];
allGLM_pupil_Boxc_stim_beta=[]; allGLM_pupil_Boxc_sust_beta=[]; allGLM_pupil_Boxc_resp_beta=[]; allGLM_pupil_Boxc_VIF=[]; allGLM_pupil_Boxc_rsquare=[];
allGLM_pupil_Ramp_stim_beta=[]; allGLM_pupil_Ramp_sust_beta=[]; allGLM_pupil_Ramp_resp_beta=[];  allGLM_pupil_Ramp_VIF=[]; allGLM_pupil_Ramp_rsquare=[];


allGLM_pupil        = cell(length(single_participants),1);% pupil time course used to fit GLM
allGLM_Boxc_yhat 	= cell(length(single_participants),1);% fit of Boxc (single trial)
allGLM_Boxc_r       = cell(length(single_participants),1);% residuals of Boxc
allGLM_Ramp_yhat    = cell(length(single_participants),1);
allGLM_Ramp_r       = cell(length(single_participants),1);

allGLM_conc_Boxc_yhat   = cell(length(single_participants),1);% fit of Boxc (concatenated)
allGLM_conc_Boxc_r      = cell(length(single_participants),1);% residuals of Boxc (concatenated)
allGLM_conc_Ramp_yhat   = cell(length(single_participants),1);
allGLM_conc_Ramp_r      = cell(length(single_participants),1);
    

clear validtrials filenames
tmpsub=0;
for isub = single_participants
    tmpsub=tmpsub+1;
    
    switch dataset
        case 'CD'
            filenames = [paths.s(isub).readbase allsubj{isub} fileExt_preprocess fileExt_GLM '.mat'];
        case 'bigDots'
            filenames = [paths.s(isub).readbase allsubj{isub} fileExt_preprocess fileExt_GLM '.mat'];
    end
    disp(['loading: ' filenames])
    clear DAT
    DAT = load([filenames]);
    
    nTrials = length(DAT.st_GLMfit.data);
    
    % GLMfit - conc trials
    allGLM_conc_pupil_BIC   = [allGLM_conc_pupil_BIC ; DAT.conc_GLMfit.BIC];
    
    allGLM_conc_pupil_Boxc_tstat        = [allGLM_conc_pupil_Boxc_tstat ; DAT.conc_GLMfit.Boxc.tstat.t(2:4)'];
    allGLM_conc_pupil_Ramp_tstat        = [allGLM_conc_pupil_Ramp_tstat ; DAT.conc_GLMfit.Ramp.tstat.t(2:4)'];
    allGLM_conc_pupil_Ramp2thresh_tstat = [allGLM_conc_pupil_Ramp2thresh_tstat ; DAT.conc_GLMfit.Ramp2thresh.tstat.t(2:4)'];
    
    allGLM_conc_pupil_Boxc_beta        = [allGLM_conc_pupil_Boxc_beta ; DAT.conc_GLMfit.Boxc.beta(2:4)'];
    allGLM_conc_pupil_Ramp_beta        = [allGLM_conc_pupil_Ramp_beta ; DAT.conc_GLMfit.Ramp.beta(2:4)'];
    allGLM_conc_pupil_Ramp2thresh_beta = [allGLM_conc_pupil_Ramp2thresh_beta ; DAT.conc_GLMfit.Ramp2thresh.beta(2:4)'];
    
    allGLM_conc_pupil_StimResp_VIF    = [allGLM_conc_pupil_StimResp_VIF ; DAT.conc_GLMfit.VIF{1}];
    allGLM_conc_pupil_Boxc_VIF        = [allGLM_conc_pupil_Boxc_VIF ; DAT.conc_GLMfit.VIF{2}];
    allGLM_conc_pupil_Ramp_VIF        = [allGLM_conc_pupil_Ramp_VIF ; DAT.conc_GLMfit.VIF{3}];
    allGLM_conc_pupil_Ramp2thresh_VIF = [allGLM_conc_pupil_Ramp2thresh_VIF ; DAT.conc_GLMfit.VIF{4}];
    
    if tmpsub==1
        % initialise binned GLM data
        allGLM_bin_pupil_Boxc_tstat = cell(length(DAT.bin_GLMfit), 1);
        allGLM_bin_pupil_Ramp_tstat = cell(length(DAT.bin_GLMfit), 1);
        allGLM_bin_pupil_Boxc_beta = cell(length(DAT.bin_GLMfit), 1);
        allGLM_bin_pupil_Ramp_beta = cell(length(DAT.bin_GLMfit), 1);
        allGLM_bin_pupil_StimResp_VIF = cell(length(DAT.bin_GLMfit), 1);
        allGLM_bin_pupil_Boxc_VIF = cell(length(DAT.bin_GLMfit), 1);
        allGLM_bin_pupil_Ramp_VIF = cell(length(DAT.bin_GLMfit), 1);
        allGLM_bin_pupil_BIC = cell(length(DAT.bin_GLMfit), 1);
    end
    
    % GLMfit - bin
    for ibin = 1:length(DAT.bin_GLMfit)
        for ibinsize = 1:length(DAT.bin_GLMfit(ibin).GLM)
            allGLM_bin_pupil_Boxc_tstat{ibin}(tmpsub, ibinsize, :)   = [ DAT.bin_GLMfit(ibin).GLM(ibinsize).Boxc.tstat.t(2:end)'];
            allGLM_bin_pupil_Ramp_tstat{ibin}(tmpsub, ibinsize, :)   = [ DAT.bin_GLMfit(ibin).GLM(ibinsize).Ramp.tstat.t(2:end)'];
            
            allGLM_bin_pupil_Boxc_beta{ibin}(tmpsub, ibinsize, :)   = [ DAT.bin_GLMfit(ibin).GLM(ibinsize).Boxc.beta(2:end)'];
            allGLM_bin_pupil_Ramp_beta{ibin}(tmpsub, ibinsize, :)   = [ DAT.bin_GLMfit(ibin).GLM(ibinsize).Ramp.beta(2:end)'];
            
            allGLM_bin_pupil_StimResp_VIF{ibin}(tmpsub, ibinsize, :)   = [ DAT.bin_GLMfit(ibin).GLM(ibinsize).VIF{1}'];
            allGLM_bin_pupil_Boxc_VIF{ibin}(tmpsub, ibinsize, :)   = [ DAT.bin_GLMfit(ibin).GLM(ibinsize).VIF{2}'];
            allGLM_bin_pupil_Ramp_VIF{ibin}(tmpsub, ibinsize, :)   = [ DAT.bin_GLMfit(ibin).GLM(ibinsize).VIF{3}'];
            
            allGLM_bin_pupil_BIC{ibin}(tmpsub, ibinsize, :)   = DAT.bin_GLMfit(ibin).GLM(ibinsize).BIC;
        end
    end
    
    % GLMfit - single trials    
    GLM_pupil_StimResp_rsquare = NaN(nTrials,1);
    GLM_pupil_StimResp_stim_beta = NaN(nTrials,1);
    GLM_pupil_StimResp_resp_beta = NaN(nTrials,1);

    GLM_pupil_Boxc_rsquare = NaN(nTrials,1);
    GLM_pupil_Boxc_stim_beta = NaN(nTrials,1);
    GLM_pupil_Boxc_sust_beta = NaN(nTrials,1);
    GLM_pupil_Boxc_resp_beta = NaN(nTrials,1);

    GLM_pupil_Ramp_rsquare = NaN(nTrials,1);
    GLM_pupil_Ramp_stim_beta = NaN(nTrials,1);
    GLM_pupil_Ramp_sust_beta = NaN(nTrials,1);
    GLM_pupil_Ramp_resp_beta = NaN(nTrials,1);

    GLM_pupil               = cell(nTrials,1);
    GLM_Boxc_yhat           = cell(nTrials,1);
    GLM_Boxc_r              = cell(nTrials,1);
    GLM_Ramp_yhat           = cell(nTrials,1);
    GLM_Ramp_r              = cell(nTrials,1);
    
    GLM_conc_Boxc_yhat      = cell(nTrials,1);
    GLM_conc_Boxc_r         = cell(nTrials,1);
    GLM_conc_Ramp_yhat      = cell(nTrials,1);
    GLM_conc_Ramp_r         = cell(nTrials,1);
    
    idx = 1;
    for itrial = 1:nTrials
        
        if ~isempty(DAT.st_GLMfit.fit(itrial).StimResp)
            GLM_pupil_StimResp_rsquare(itrial)       = DAT.st_GLMfit.fit(itrial).StimResp.rsquare;
            GLM_pupil_StimResp_stim_beta(itrial) = DAT.st_GLMfit.fit(itrial).StimResp.beta(2);
            GLM_pupil_StimResp_resp_beta(itrial) = DAT.st_GLMfit.fit(itrial).StimResp.beta(3);
        end
        if ~isempty(DAT.st_GLMfit.fit(itrial).Boxc)
            GLM_pupil_Boxc_rsquare(itrial) = DAT.st_GLMfit.fit(itrial).Boxc.rsquare;
            GLM_pupil_Boxc_stim_beta(itrial) = DAT.st_GLMfit.fit(itrial).Boxc.beta(2);
            GLM_pupil_Boxc_sust_beta(itrial) = DAT.st_GLMfit.fit(itrial).Boxc.beta(3);
            GLM_pupil_Boxc_resp_beta(itrial) = DAT.st_GLMfit.fit(itrial).Boxc.beta(4);
        end
        if ~isempty(DAT.st_GLMfit.fit(itrial).Ramp)
            GLM_pupil_Ramp_rsquare(itrial) = DAT.st_GLMfit.fit(itrial).Ramp.rsquare;
            GLM_pupil_Ramp_stim_beta(itrial) = DAT.st_GLMfit.fit(itrial).Ramp.beta(2);
            GLM_pupil_Ramp_sust_beta(itrial) = DAT.st_GLMfit.fit(itrial).Ramp.beta(3);
            GLM_pupil_Ramp_resp_beta(itrial) = DAT.st_GLMfit.fit(itrial).Ramp.beta(4);
        end
        
        if ~isempty(DAT.st_GLMfit.data(itrial).Pupil)
            % st
            GLM_pupil{itrial}       = DAT.st_GLMfit.data(itrial).Pupil;
            GLM_Boxc_yhat{itrial}   = DAT.st_GLMfit.fit(itrial).Boxc.yhat;
            GLM_Boxc_r{itrial}      = DAT.st_GLMfit.fit(itrial).Boxc.r;
            
            GLM_Ramp_yhat{itrial}   = DAT.st_GLMfit.fit(itrial).Ramp.yhat;
            GLM_Ramp_r{itrial}      = DAT.st_GLMfit.fit(itrial).Ramp.r;
            
            % conc
            GLM_conc_Boxc_yhat{itrial}  = DAT.conc_GLMfit.Boxc.yhat(idx:(idx+length(GLM_pupil{itrial})-1));
            GLM_conc_Boxc_r{itrial}     = DAT.conc_GLMfit.Boxc.r(idx:(idx+length(GLM_pupil{itrial})-1));
            
            GLM_conc_Ramp_yhat{itrial}  = DAT.conc_GLMfit.Ramp.yhat(idx:(idx+length(GLM_pupil{itrial})-1));
            GLM_conc_Ramp_r{itrial}     = DAT.conc_GLMfit.Ramp.r(idx:(idx+length(GLM_pupil{itrial})-1));
            
            idx = idx + length(GLM_pupil{itrial});
        end
    end
    
    allGLM_pupil_StimResp_rsquare = [allGLM_pupil_StimResp_rsquare     ; GLM_pupil_StimResp_rsquare];
    allGLM_pupil_StimResp_stim_beta = [allGLM_pupil_StimResp_stim_beta     ; GLM_pupil_StimResp_stim_beta];
    allGLM_pupil_StimResp_resp_beta = [allGLM_pupil_StimResp_resp_beta     ; GLM_pupil_StimResp_resp_beta];
    
    allGLM_pupil_Boxc_rsquare = [allGLM_pupil_Boxc_rsquare     ; GLM_pupil_Boxc_rsquare];
    allGLM_pupil_Boxc_stim_beta = [allGLM_pupil_Boxc_stim_beta     ; GLM_pupil_Boxc_stim_beta];
    allGLM_pupil_Boxc_sust_beta = [allGLM_pupil_Boxc_sust_beta     ; GLM_pupil_Boxc_sust_beta];
    allGLM_pupil_Boxc_resp_beta = [allGLM_pupil_Boxc_resp_beta     ; GLM_pupil_Boxc_resp_beta];
    
    allGLM_pupil_Ramp_rsquare = [allGLM_pupil_Ramp_rsquare     ; GLM_pupil_Ramp_rsquare];
    allGLM_pupil_Ramp_stim_beta = [allGLM_pupil_Ramp_stim_beta     ; GLM_pupil_Ramp_stim_beta];
    allGLM_pupil_Ramp_sust_beta = [allGLM_pupil_Ramp_sust_beta     ; GLM_pupil_Ramp_sust_beta];
    allGLM_pupil_Ramp_resp_beta = [allGLM_pupil_Ramp_resp_beta     ; GLM_pupil_Ramp_resp_beta];
    
    allGLM_pupil_StimResp_VIF = [allGLM_pupil_StimResp_VIF     ; DAT.st_GLMfit.VIF.StimResp];
    allGLM_pupil_Boxc_VIF = [allGLM_pupil_Boxc_VIF     ; DAT.st_GLMfit.VIF.Boxc];
    allGLM_pupil_Ramp_VIF = [allGLM_pupil_Ramp_VIF     ; DAT.st_GLMfit.VIF.Ramp];
    
    % also get pupil data that was used to fit the model, as well as fit
    % and residuals
    pupilLength     = cellfun(@(X) (length(X)), GLM_pupil, 'UniformOutput', true);
    allGLM_pupil{tmpsub}   = NaN(max(pupilLength), length(GLM_pupil));
    
    allGLM_Boxc_yhat{tmpsub}    = NaN(max(pupilLength), length(GLM_pupil));
    allGLM_Boxc_r{tmpsub}       = NaN(max(pupilLength), length(GLM_pupil));
    allGLM_Ramp_yhat{tmpsub}    = NaN(max(pupilLength), length(GLM_pupil));
    allGLM_Ramp_r{tmpsub}       = NaN(max(pupilLength), length(GLM_pupil));
    
    allGLM_conc_Boxc_yhat{tmpsub}    = NaN(max(pupilLength), length(GLM_pupil));
    allGLM_conc_Boxc_r{tmpsub}       = NaN(max(pupilLength), length(GLM_pupil));
    allGLM_conc_Ramp_yhat{tmpsub}    = NaN(max(pupilLength), length(GLM_pupil));
    allGLM_conc_Ramp_r{tmpsub}       = NaN(max(pupilLength), length(GLM_pupil));
    
    
    for itrial = 1:nTrials
        if ~isempty(GLM_pupil{itrial})
            allGLM_pupil{tmpsub}(1:length(GLM_pupil{itrial}), itrial) = GLM_pupil{itrial};
            
            allGLM_Boxc_yhat{tmpsub}(1:length(GLM_pupil{itrial}), itrial) = GLM_Boxc_yhat{itrial};
            allGLM_Boxc_r{tmpsub}(1:length(GLM_pupil{itrial}), itrial) = GLM_Boxc_r{itrial};
            allGLM_Ramp_yhat{tmpsub}(1:length(GLM_pupil{itrial}), itrial) = GLM_Ramp_yhat{itrial};
            allGLM_Ramp_r{tmpsub}(1:length(GLM_pupil{itrial}), itrial) = GLM_Ramp_r{itrial};
            
            allGLM_conc_Boxc_yhat{tmpsub}(1:length(GLM_pupil{itrial}), itrial) = GLM_conc_Boxc_yhat{itrial};
            allGLM_conc_Boxc_r{tmpsub}(1:length(GLM_pupil{itrial}), itrial) = GLM_conc_Boxc_r{itrial};
            allGLM_conc_Ramp_yhat{tmpsub}(1:length(GLM_pupil{itrial}), itrial) = GLM_conc_Ramp_yhat{itrial};
            allGLM_conc_Ramp_r{tmpsub}(1:length(GLM_pupil{itrial}), itrial) = GLM_conc_Ramp_r{itrial};
        end
    end

    clear DAT
end

%% save

filename_mat = [paths.pop 'allSub' fileExt_preprocess fileExt_GLM '.mat'];
disp(['saving ' filename_mat])
save([filename_mat],'allGLM*','-v7.3')

