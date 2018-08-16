function analyse_CDT(isub, single_participants, dataset, fileExt)
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
% analyse single trial data, for the continuous dots task

close all
tic
switch dataset
    case 'CD'
        getSubSpecs_CD
        getfilenames_CD
        channelConfig = 'biosemi96';
    case 'bigDots'
        getSubSpecs_bigDots
        getfilenames_bigDots
        channelConfig = 'biosemi64';
end

% setPlotSettings
setAnalysisSettings_bigDots

if ~isempty(single_participants)
    subject_folder  = subject_folder(single_participants);
    allsubj         = allsubj(single_participants);
end%

%% Define channels, having combined Brain Products and Biosemi data

switch channelConfig
    case 'BP'
        chStr.left_hemi = {...
            'FP1'    'AF7'    'AF3'    'F1'    'F3'    'F5'    'F7' ...
            'FT9'    'FT7'    'FC5'    'FC3'   'FC1'   'C1'    'C3'    'C5' ...
            'T7'     'TP9'    'TP7'    'CP5'   'CP3'   'CP1'   'P1'    'P3' ...
            'P5'     'P7'     'PO9'    'PO7'   'PO3'   'O1'};
        chStr.right_hemi = {...
            'FP2'    'AF8'    'AF4'    'F2'    'F4'    'F6'    'F8' ...
            'FT10'   'FT8'    'FC6'    'FC4'   'FC2'   'C2'    'C4'    'C6' ...
            'T8'     'TP10'   'TP8'    'CP6'   'CP4'   'CP2'   'P2'    'P4' ...
            'P6'     'P8'     'PO10'   'PO8'   'PO4'   'O2'};
        chStr.centre_chans = {...
            'Oz'    'POz'    'Pz'    'CPz'    'Fz'    'FCz'    'Cz'};
        chStr.elec_pairs = {...
            'FP1', 'FP2' ; 'AF7' ,'AF8' ; 'AF3' ,'AF4' ; 'F1'  ,'F2'  ;...
            'F3' , 'F4'  ; 'F5'  ,'F6'  ; 'F7'  ,'F8'  ; 'FT9' ,'FT10'; 'FT7' ,'FT8' ;...
            'FC5', 'FC6' ; 'FC3' ,'FC4' ; 'FC1' ,'FC2' ; 'C1'  ,'C2'  ;...
            'C3' , 'C4'  ; 'C5'  ,'C6'  ; 'T7'  ,'T8'  ; 'TP9' ,'TP10'; 'TP7' ,'TP8' ;...
            'CP5', 'CP6' ; 'CP3' ,'CP4' ; 'CP1' ,'CP2' ; 'P1'  ,'P2'  ;...
            'P3' , 'P4'  ; 'P5'  ,'P6'  ; 'P7'  ,'P8'  ; 'PO9' ,'PO10'; 'PO7' ,'PO8' ;...
            'PO3', 'PO4' ; 'O1'  ,'O2'};
    case 'biosemi64'
        chStr.left_hemi = {...
            'FP1'    'AF7'    'AF3'    'F1'    'F3'    'F5'    'F7' ...
            'FT7'    'FC5'    'FC3'    'FC1'    'C1'    'C3'    'C5' ...
            'T7'    'TP7'    'CP5'    'CP3'    'CP1'    'P1'    'P3' ...
            'P5'    'P7'    'PO7'    'PO3'    'O1'};
        chStr.right_hemi = {...
            'FP2'    'AF8'    'AF4'    'F2'    'F4'    'F6'    'F8' ...
            'FT8'    'FC6'    'FC4'    'FC2'    'C2'    'C4'    'C6' ...
            'T8'    'TP8'    'CP6'    'CP4'    'CP2'    'P2'    'P4' ...
            'P6'    'P8'    'PO8'    'PO4'    'O2'};
        chStr.centre_chans = {...
            'Oz'    'POz'    'Pz'    'CPz'    'Fz'    'FCz'    'Cz'};
        chStr.elec_pairs = {...
            'FP1', 'FP2' ; 'AF7' ,'AF8' ; 'AF3' ,'AF4' ; 'F1'  ,'F2'  ;...
            'F3' , 'F4'  ; 'F5'  ,'F6'  ; 'F7'  ,'F8'  ; 'FT7' ,'FT8' ;...
            'FC5', 'FC6' ; 'FC3' ,'FC4' ; 'FC1' ,'FC2' ; 'C1'  ,'C2'  ;...
            'C3' , 'C4'  ; 'C5'  ,'C6'  ; 'T7'  ,'T8'  ; 'TP7' ,'TP8' ;...
            'CP5', 'CP6' ; 'CP3' ,'CP4' ; 'CP1' ,'CP2' ; 'P1'  ,'P2'  ;...
            'P3' , 'P4'  ; 'P5'  ,'P6'  ; 'P7'  ,'P8'  ; 'PO7' ,'PO8' ;...
            'PO3', 'PO4' ; 'O1'  ,'O2'};
        
    case 'biosemi96'
        
        chStr.left_hemi = {...
            'FP1'    'AF7'    'AFP1'   'AF3'    'F1'    'F3'    'F5'    'F7'    'F9' ...
            'FFC1h'  'FFC5h'  'FT9'    'FT7'    'FC5'    'FC3'    'FC1'    'C1'    'C3'    'C5' ...
            'T7'    'TP9'     'TP7'    'CP5'    'CP3'    'CP1'    'TPP7h'  'CPP5h' 'CPP3h' 'CPP1h' ...
            'P1'    'P3'      'P5'     'P7'     'P9' ...
            'PPO9h' 'PPO5h'   'PPO1h'  'PO7'    'PO3'    ...
            'PO9'   'POO9h'   'O1'     'POO1'   'O9'  'OI1h'};
        
        chStr.right_hemi = {...
            'FP2'    'AF8'    'AFP2'   'AF4'    'F2'    'F4'    'F6'    'F8'    'F10' ...
            'FFC2h'  'FFC6h'  'FT10'   'FT8'    'FC6'    'FC4'    'FC2'    'C2'    'C4'    'C6' ...
            'T8'     'TP10'   'TP8'    'CP6'    'CP4'    'CP2'    'TPP8h'  'CPP6h' 'CPP4h' 'CPP2h' ...
            'P2'    'P4'      'P6'     'P8'     'P10' ...
            'PPO10h' 'PPO6h'   'PPO2h'  'PO8'    'PO4'    ...
            'PO10'   'POO10h'   'O2'     'POO2'   'O10'  'OI2h'};
        
        chStr.centre_chans = {...
            'Oz'    'POz'    'Pz'    'CPz'   'Cz'   'Fz' };
        
        chStr.elec_pairs = {...
            'FP1', 'FP2' ; 'AF7' ,'AF8' ; 'AFP1','AFP2' ;'AF3' ,'AF4' ; 'F1'  ,'F2'  ;...
            'F3' , 'F4'  ; 'F5'  ,'F6'  ; 'F7'  ,'F8'   ;'F9' , 'F10' ; 'FFC1h','FFC2h' ; 'FFC5h','FFC6h';...
            'FT9', 'FT10'; 'FT7' ,'FT8' ; 'FC5', 'FC6' ; 'FC3' ,'FC4' ; 'FC1' ,'FC2' ; ...
            'C1'  ,'C2'  ; 'C3' , 'C4'  ; 'C5'  ,'C6'  ; 'T7'  ,'T8'  ; 'TP9','TP10' ; 'TP7' ,'TP8' ;...
            'CP5', 'CP6' ; 'CP3' ,'CP4' ; 'CP1' ,'CP2' ; 'TPP7h','TPP8h' ; 'CPP5h', 'CPP6h'; 'CPP3h','CPP4h'; 'CPP1h','CPP2h';
            'P1'  ,'P2'  ; 'P3' , 'P4'  ; 'P5'  ,'P6'  ; 'P7'  ,'P8'  ; 'P9','P10' ; ...
            'PPO9h','PPO10h' ; 'PPO5h','PPO6h' ; 'PPO1h', 'PPO2h' ; 'PO7' ,'PO8' ;...
            'PO3', 'PO4' ; 'PO9','PO10'; 'POO9h','POO10h' ; 'O1'  ,'O2' ; 'POO1','POO2' ; 'O9','O10'; 'OI1h','OI2h' };
end

%
ch.left_hemi    = getChannelName(chStr.left_hemi,   channelConfig);
ch.right_hemi   = getChannelName(chStr.right_hemi,  channelConfig);
ch.centre_chans = getChannelName(chStr.centre_chans,channelConfig);
ch.elec_pairs   = getChannelName(chStr.elec_pairs,  channelConfig);

switch dataset
    case 'bigDots'
        [~,plot_chans, exclude_chans] = getChannelName;
    case 'CD'
        plot_chans = 1:96;
end
nChan=max(plot_chans);

if 0
    % plot channel locations on empty topoplot
    [~,plot_chans, exclude_chans] = getChannelName;
    tester = zeros(max(plot_chans),1);
    figure
    chanlocs = readlocs('cap64.loc'); %biosemi
    topoplot(tester,chanlocs,'maplimits', ...
    [min(tester)  max(tester)],'electrodes','labels','plotchans',plot_chans);
end
%% Triggers

switch dataset
    case 'bigDots'
        % Trigger 1: coherence 50, motion dir 270, ITI 3.06, patch 1
        % Trigger 2: coherence 50, motion dir 270, ITI 3.06, patch 2
        % Trigger 3: coherence 50, motion dir 270, ITI 5.17, patch 1
        % Trigger 4: coherence 50, motion dir 270, ITI 5.17, patch 2
        % Trigger 5: coherence 50, motion dir 270, ITI 7.29, patch 1
        % Trigger 6: coherence 50, motion dir 270, ITI 7.29, patch 2
        
        % ITI,left/right
        targcodes = zeros(3,2);
        targcodes(1,:) = [101,102];
        targcodes(2,:) = [103,104];
        targcodes(3,:) = [105,106];
        
    case 'CD'
        targcodes = 8;
end

%% get subject index
subIdx.sub  = single_participants(isub); %get matching subject
% paths.fig = [paths.s(subIdx.sub).fig 'singleTrial' filesep];
% if ~exist(paths.fig,'dir'),mkdir(paths.fig),end

%% get filenames and drug information, if applicable
switch dataset
    case 'CD'
        for iside = 1:2 %up down
            loadfilenames{iside} = [paths.data subject_folder{isub} filesep allsubj{isub} '_' dataset '_' side_tags{iside}(1) '_final.mat'];
        end
        savefilename = [paths.s(subIdx.sub).base allsubj{isub} fileExt '_ST.mat'];
    case 'bigDots'
        loadfilenames{1} = [paths.data subject_folder{isub} filesep allsubj{isub} 'big_dots_erp_final.mat'];
        savefilename = [paths.s(subIdx.sub).base allsubj{isub} fileExt '_ST.mat'];
end

%% start loop over separate data files
for ifile = 1:length(loadfilenames)
    %% load data
    load(loadfilenames{ifile});
    
    disp(['Subject ' num2str(isub) ': ' allsubj{isub} ' number of trials = ' num2str(length(find(allTrig)))])
    
    if strcmp(subject_folder{isub},'331M_CL') % really odd tiny artifact meant this trial was messing with CSD!
        allRT(53) = 0; allrespLR(53) = 0; allTrig(53) = 999;
    end
    
    %% initialise per subject
    initialise_variables
    
    %% ERP
    tmp_erp_csd=double(erp_LPF_35Hz_CSD);%used for CPP
    %     tmp_erp_csd=double(erp_LPF_8Hz_CSD);%used for CPP
    
    if doCSD
        tmp_erp=double(erp_LPF_35Hz_CSD);
        %         tmp_erp=double(erp_LPF_8Hz_CSD);
    else
        tmp_erp=double(erp_LPF_35Hz);
        %         tmp_erp=double(erp_LPF_8Hz);
    end
    
    % Baseline erp
    baseline_erp        = mean(tmp_erp(:,find(t>=BL_time(1)     & t<=BL_time(2)),:),2);
    baseline_erp_csd    = mean(tmp_erp_csd(:,find(t>=BL_time(1) & t<=BL_time(2)),:),2);
    
    tmp_erp                 = tmp_erp           - repmat(baseline_erp,      [1,size(tmp_erp,2),1]); % baseline full erp
    tmp_erp_csd             = tmp_erp_csd       - repmat(baseline_erp_csd,  [1,size(tmp_erp_csd,2),1]); % baseline full erp
    
    % get response locked erp and pupil
    tmp_erpr    = zeros(size(tmp_erp,1),length(tr),size(tmp_erp,3));
    tmp_erpr_csd= zeros(size(tmp_erp_csd,1),length(tr),size(tmp_erp_csd,3));
    for n=1:length(allRT)
        [~,RTsamp] = min(abs(t*fs/1000-allRT(n))); % get the sample point of the RT.
        RTsamp_idx(n) = RTsamp;
        if RTsamp+trs(1) >0 && RTsamp+trs(length(trs))<=length(t) && allRT(n)>0 % isub the RT larger than 1st stim RT point, smaller than last RT point.
            tmp_erpr(:,:,n)     = tmp_erp(:,RTsamp+trs,n);
            tmp_erpr_csd(:,:,n) = tmp_erp_csd(:,RTsamp+trs,n);
            validrlock(idx_trFile(n),1)=1;
        end
    end
    
    %% ERP
    
    % if not necessary, don't create these. They will be the biggest files
    if rawERP
        ERP(:,:,idx_trFile)       = tmp_erp;
        ERP_csd(:,:,idx_trFile)   = tmp_erp_csd;
        ERPr(:,:,idx_trFile)      = tmp_erpr;
        ERPr_csd(:,:,idx_trFile)  = tmp_erpr_csd;
    end
    
    %% trial indices & RT
    subject(idx_trFile,1)      = isub;
    hits(idx_trFile,1)         = allrespLR'==1;
    misses(idx_trFile,1)       = allrespLR'==3;
    
    subRT(idx_trFile,1)        = allRT*1000/fs;
    validRT(idx_trFile(allRT>(rtlim(1)*fs) & allRT<(rtlim(2)*fs)), 1) = true;
    
    validtr.validRT = false(nTrFile,1);
    validtr.validRT (allrespLR'==1 & validRT(idx_trFile,1) & validrlock(idx_trFile,1) ,1) = true; %
    
    trIdx = idx_trFile(validRT(idx_trFile));
    subRT_log(trIdx,1)    = log(subRT(trIdx,1));
    subRT_zscore(trIdx,1) = zscore(subRT_log(trIdx));
    
    % get condition indices side of target and iti
    switch dataset
        case {'bigDots'}
            for isideStim = 1:2 %target stimulus side
                for iti = 1:3
                    % calcs the indices of the triggers for each appropriate trial type.
                    sideStimtr  (idx_trFile(allTrig==targcodes(iti,isideStim)),1) = isideStim;
                    ititr       (idx_trFile(allTrig==targcodes(iti,isideStim)),1) = iti;
                    
                    % these subjects did not have iti 3
                    if (strcmp(subject_folder{isub},'AR_08_04_14') || strcmp(subject_folder{isub},'MH_14_04_14')) && iti==3
                        validtr.validRT (allTrig==targcodes(iti,isideStim),1)  = 0;
                        ititr(idx_trFile(allTrig==targcodes(iti,isideStim)),1) = 0;
                    end
                end
            end
        case 'CD'
            sideStimtr(idx_trFile,1) = allUD';
            ititr(idx_trFile,1)      = allITI';
    end
    
    % extract info from trials preceding and following current trial
    allBlocks = unique(blockIdx);
    for iblock = 1:length(allBlocks)
        trIdx = find(blockIdx==allBlocks(iblock));
        
        for itrial = 1:length(find(trIdx))
            saveIdx = 1:RTwinsize*2+1;
            trwin = trIdx(itrial) + [-RTwinsize:RTwinsize];
            saveIdx = saveIdx(trwin >= trIdx(1) & trwin <= trIdx(end));
            trwin = trwin(trwin >= trIdx(1) & trwin <= trIdx(end));
            
            subRTz_window(idx_trFile(trIdx(itrial)), saveIdx) = subRT_zscore(trwin);
            subRT_window(idx_trFile(trIdx(itrial)), saveIdx) = subRT(trwin);
            validRT_window(idx_trFile(trIdx(itrial)), saveIdx) = validRT(trwin);
        end
    end
    
    switch dataset
        case {'bigDots','CD'}
            
            for ifield = 1:length(artifact_times)
                eval(['validtr.' artifact_times{ifield} '(idx_trFile,1) = false(length(idx_trFile),1);'])
                eval(['validtr_eye.' artifact_times{ifield} '(idx_trFile,1) = false(length(idx_trFile),1);'])
                
                eval(['validtr.' artifact_times{ifield} '(idx_trFile( validtr.validRT & ' ...
                    ' transpose(artrej.' artifact_times{ifield} ') & ' ...
                    'transpose(artrej_ET.' artifact_times{ifield} ') & artrej_pupil.' artifact_times{ifield} '),1) = true;'])
                
                eval(['validtr_eye.' artifact_times{ifield} '(idx_trFile( validtr.validRT & ' ...
                    'transpose(artrej_ET.' artifact_times{ifield} ') & artrej_pupil.' artifact_times{ifield} '),1) = true;'])
            end
    end
    
    
    validtrials(ifile) = round(100*(length(find(validtr.neg100_RT_200(idx_trFile,1))))/(nTrFile));
    disp(['Subject ',allsubj{isub},' Total Valid Trials: ',num2str(length(find(validtr.neg100_RT_200(idx_trFile,1)))), ...
        ' = ',num2str(validtrials(ifile)),'%'])
    
    %% N2c & N2i
    
    clear ch_N2
    ch_N2 = getChannelName(ch_N2_str, channelConfig);
    
    for isideStim = 1:2
        ch_N2c      = ch_N2(isideStim);
        ch_N2_flip  = flip(ch_N2);
        ch_N2i      = ch_N2_flip(isideStim);
        
        trIdx = (sideStimtr(idx_trFile,1) == isideStim);
        N2c(:,idx_trFile(trIdx),1)     = squeeze(mean(tmp_erp(ch_N2c,:,trIdx),1));
        N2i(:,idx_trFile(trIdx),1)     = squeeze(mean(tmp_erp(ch_N2i,:,trIdx),1));
    end
    
    tIdx = find(t>=200 & t<=300);
    
    try
        N2c_topo(:,idx_trFile)    = squeeze(mean(tmp_erp(plot_chans,tIdx,:),2)); %get CPP response locked
    catch
        N2c_topo(:,idx_trFile)    = squeeze(mean(tmp_erp(:,tIdx,:),2)); %get CPP response locked
        
    end
    %         figure,clf
    %         subplot(1,2,1)
    %         topoplot(squeeze(mean(N2c_topo(:,sideStimtr==1),2)),chanlocs,'electrodes','off','plotchans',plot_chans,'numcontour',0);
    %         subplot(1,2,2)
    %         topoplot(squeeze(mean(N2c_topo(:,sideStimtr==2),2)),chanlocs,'electrodes','off','plotchans',plot_chans,'numcontour',0);
    %
    
    
    %% CPP
    
    % select the channel used for CPP and get ERPs from this/these channel(s)
    clear ch_CPP
    ch_CPP = getChannelName(ch_CPP_str, channelConfig);
    
    CPP(:,idx_trFile)     = squeeze(mean(tmp_erp(ch_CPP,:,:),1)); % get CPP, mean across channels if there are more than 1
    CPPr(:,idx_trFile)    = squeeze(mean(tmp_erpr(ch_CPP,:,:),1)); % get CPP, response-locked
    
    tIdx = find(trs>=-50 & trs<=50);
    try
        CPP_topo(:,idx_trFile)    = squeeze(mean(tmp_erpr(plot_chans,tIdx,:),2)); 
    catch
        CPP_topo(:,idx_trFile)    = squeeze(mean(tmp_erpr(:,tIdx,:),2)); 
    end
    %     topoplot(squeeze(mean(CPP_topo,2)),chanlocs,'electrodes','off','plotchans',plot_chans,'numcontour',0);
    
    CPP_csd(:,idx_trFile)     = squeeze(mean(tmp_erp_csd(ch_CPP,:,:),1)); % get CPP, mean across channels if there are more than 1
    CPPr_csd(:,idx_trFile)    = squeeze(mean(tmp_erpr_csd(ch_CPP,:,:),1)); %get CPP response locked
    
    
    %% Alpha
    
    % get spectral tranformation
    [spectral_t,  ~, alpha_TSE] = compute_SpectrotemporalEvolution(tmp_erp,bandlimits_alpha,ch,t,tSpectral,fs,BL_spectrum);
    
    %select the channel used for occipital alpha
    for isideCh = 1:2 % hemisphere, left and right
        ch_alpha(isub,isideCh,:) = getChannelName(mean_ch_alpha(abs(isideCh-3),:),channelConfig);
    end
    
    % get alpha
    tIdx = spectral_t>t_preTarget_alpha(1) & spectral_t<t_preTarget_alpha(2);
    alpha_preTarget_topo(:,idx_trFile) = squeeze(mean(alpha_TSE(1:nChan,tIdx,:),2));
    
    % do this for all trials, also the non-valid ones, this will be rejected after
    for isideCh = 1:2 % hemisphere, left and right, average over predetermined channels for each hemisphere
        if sum(ismember(ch_alpha(isub,isideCh,:), ch.right_hemi)) == nch_alpha % check that the correct channels will be put in the correct variables (matching right with right)
            alphaRh_preTarget(idx_trFile,1)     = squeeze(mean(mean(alpha_TSE(ch_alpha(isub,isideCh,:),tIdx,:),1),2));
        elseif sum(ismember(ch_alpha(isub,isideCh,:), ch.left_hemi)) == nch_alpha
            alphaLh_preTarget(idx_trFile,1)     = squeeze(mean(mean(alpha_TSE(ch_alpha(isub,isideCh,:),tIdx,:),1),2));
        end
    end
    alpha_preTarget(idx_trFile,1)   = squeeze(mean(mean(alpha_TSE(ch_alpha(isub,:,:),tIdx,:),1),2));
    
    %pretarget alpha asymmetry, asymmetry: right minus left. more positive = more right hemi alpha (= less dysynchronisation)
    alphaAsym_preTarget(idx_trFile,1) = ...
        (alphaRh_preTarget(idx_trFile,1) - alphaLh_preTarget(idx_trFile,1)) ./ ...
        ((alphaRh_preTarget(idx_trFile,1) + alphaLh_preTarget(idx_trFile,1))/2);
    
    %pretarget alpha asym, keeping all channel info, for topoplots
    alphaAsym_preTarget_topo(ch.elec_pairs(:,2),idx_trFile) = ...
        (squeeze(mean(alpha_TSE(ch.elec_pairs(:,2),tIdx,:),2)) - squeeze(mean(alpha_TSE(ch.elec_pairs(:,1),tIdx,:),2))) ./...
        ((squeeze(mean(alpha_TSE(ch.elec_pairs(:,2),tIdx,:),2)) + squeeze(mean(alpha_TSE(ch.elec_pairs(:,1),tIdx,:),2)))/2);
    
    %% Beta
    if beta_power == 1
        % get spectral tranformation
        [spectral_t, ~, beta_TSE]  = compute_SpectrotemporalEvolution(tmp_erp,bandlimits_beta,ch,t,tSpectral,fs,BL_spectrum);
        [spectral_tr,~, betar_TSE] = compute_SpectrotemporalEvolution(tmp_erpr,bandlimits_beta,ch,tr,tSpectral,fs,BL_spectrum);%lock beta to response
        
    elseif beta_power == 2
        % STFT calculation
        disp('Calculating STFT...')
        tmp_erp = double(tmp_erp); % chan x time x trial
        STFT = [];
        for itrial = 1:size(tmp_erp,3)
            cc=1;
            for tt_beta = 1:skip_step:size(tmp_erp,2)-(stftlen_STFT)
                tf = tt_beta:tt_beta+stftlen_STFT-1; % define time window
                ep = squeeze(tmp_erp(:,tf,itrial)); % chop out chan x time window
                nfft = size(ep,2);
                ep = detrend(ep')'; % detrend
                fftx = abs(fft(ep,[],2))./(stftlen_STFT/2);
                fftx = fftx(:,1:ceil((nfft+1)/2));
                ind = find(freq_temp_STFT>fs_STFT(1) & freq_temp_STFT<fs_STFT(end));
                %         [~,ind] = min(abs(freq_temp_STFT-fs_STFT)); % if you want SSVEP at one particular frequency
                STFT(:,cc,itrial) = mean(fftx(:,ind),2);
                cc=cc+1;
            end
        end
        % rename, to keep consistent whatever method is used
        beta_TSE = STFT;
%         beta_TSE_base = STFT-repmat(baseline_beta,[1,size(STFT,2),1]); % baseline full erp
        
        STFT_timer= -600:unique(diff(STFT_time)):80;
        %Response locked STFT time in samples
        STFT_timers = (1:length(STFT_timer))-length(find(STFT_timer<=0));
        STFTr = zeros(size(STFT,1),length(STFT_timer),size(STFT,3));
%         STFTr_base = zeros(size(STFT,1),length(STFT_timer),size(STFT,3));
        
        for itrial = 1:size(tmp_erp,3)
            if ~validrlock(idx_trFile(itrial))
                continue
            end
            [~,RTsamp] = min(abs(STFT_time*fs/1000-allRT(itrial))); % get the sample point of the RT.
            STFTr(:,:,itrial) = STFT(:,RTsamp+STFT_timers,itrial);
%             STFTr_base(:,:,itrial) = beta_TSE_base(:,RTsamp+STFT_timers,itrial);
        end
        
        betar_TSE = STFTr;
%         betar_TSE_base = STFTr_base;
        
        spectral_t = STFT_time;
        spectral_tr = STFT_timer;
    end
    
    %select the channel used for beta (motor response)
    ch_beta(isub) = getChannelName(ch_beta_str, channelConfig);
    
    % stim locked beta
    beta_postTarget(:,idx_trFile) = squeeze(mean(beta_TSE(ch_beta(isub), :, :),1));
    
    % response locked beta
    beta_preResponse(:,idx_trFile) = squeeze(mean(betar_TSE(ch_beta(isub), :, :),1));
    
    tIdx = spectral_tr>t_preResponse_beta(1) & spectral_tr<t_preResponse_beta(2);
    beta_preResponse_mean(idx_trFile,1)      = squeeze(mean(mean(betar_TSE(ch_beta(isub), tIdx, :),1),2));
        
    %post target beta, keeping all channel info, for topoplots
    tIdx = spectral_t>t_postTarget_beta(1) & spectral_t<t_postTarget_beta(2);
    beta_postTarget_topo(1:nChan,idx_trFile) = squeeze(mean(beta_TSE(1:nChan, tIdx, :),2));
    
    %post target beta, keeping all channel info, for topoplots
    tIdx = spectral_tr>t_preResponse_beta(1) & spectral_tr<t_preResponse_beta(2);
    beta_preResponse_topo(1:nChan,idx_trFile) = squeeze(mean(betar_TSE(1:nChan, tIdx, :),2));
    
    
    
    %% SPG chronux
    
    % N2
    [SPG.N2c_power,  SPG.tt,   SPG.ff]    = mtspecgramc(N2c, movingwin,   SPG.params);
    SPG.N2c_power = log10(SPG.N2c_power);
    [SPG.N2c_phase,  SPG.tt,   SPG.ff]    = mtspecgramc_phase(N2c,movingwin,   SPG.params);
    
    [SPG.N2i_power,  SPG.tt,   SPG.ff]    = mtspecgramc(N2i, movingwin,   SPG.params);
    SPG.N2i_power = log10(SPG.N2i_power);
    [SPG.N2i_phase,  SPG.tt,   SPG.ff]    = mtspecgramc_phase(N2i,movingwin,   SPG.params);        %
    
    % shorter time window for N2, control
    [SPG.N2c_256_power,  SPG.tt_256,   SPG.ff_256]    = mtspecgramc(N2c, movingwin_N2,   SPG.params);
    SPG.N2c_256_power = log10(SPG.N2c_256_power);
    [SPG.N2c_256_phase,  SPG.tt_256,   SPG.ff_256]    = mtspecgramc_phase(N2c, movingwin_N2,   SPG.params);
    
    [SPG.N2i_256_power,  SPG.tt_256,   SPG.ff_256]    = mtspecgramc(N2i, movingwin_N2,   SPG.params);
    SPG.N2i_256_power = log10(SPG.N2i_256_power);
    [SPG.N2i_256_phase,  SPG.tt_256,   SPG.ff_256]    = mtspecgramc_phase(N2i, movingwin_N2,   SPG.params);        %
    
    %CPP
    [SPG.CPP_power,  SPG.tt,   SPG.ff]    = mtspecgramc(CPP, movingwin,   SPG.params);
    SPG.CPP_power = log10(SPG.CPP_power);
    [SPG.CPP_phase,  SPG.tt,   SPG.ff]    = mtspecgramc_phase(CPP,movingwin,   SPG.params);
    
    %% pupil
    allET_trials(:,:,idx_trFile) = ET_trials;
end

%% get pupil, only process after loading of all files (to normalise to the average between the files)

[pupil, pupilr] = processPupilDiameter(allET_trials, pupilSet, validtr.neg100_RT_200, subRT);

%%% get linear projection
%%% response locked
tIdx = find(tr>=-500 & tr<=200);
validtr2use = validtr.neg100_RT_200;
lpPupil.lp.resp_locked_neg500_200(validtr2use) = linearProjection(pupilr.lp(tIdx,validtr2use),tr(tIdx));
lpPupil.bp.resp_locked_neg500_200(validtr2use) = linearProjection(pupilr.bp(tIdx,validtr2use),tr(tIdx));



%% save
save(savefilename, ...
    'trialIdx','blockIdx','blockTrialIdx','sideStimtr','motiontr','ititr','hits','misses','validrlock','subject',...
    'validtr','validtr_eye',...
    'subRT','subRT_log','subRT_zscore','validRT',...
    'subRT_window','subRTz_window','validRT_window',...
    'N2c','N2i','N2c_topo','N2i_topo',...
    'CPP','CPP_csd','CPPr','CPPr_csd','CPP_topo',...
    'spectral_t','spectral_tr','alpha_preTarget','alphaRh_preTarget','alphaLh_preTarget','alpha_preTarget_topo','alphaAsym_preTarget','alphaAsym_preTarget_topo',...
    'beta_postTarget','beta_preResponse','beta_preResponse_topo','beta_postTarget_topo', 'beta_preResponse_mean', ...
    'SPG',...
    'pupil','pupilr','lpPupil',...
    '-v7.3')

toc
