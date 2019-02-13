function analysis_PsPM_pupil(subject, RT, blockIdx, blockTrialIdx, validtr2use, tt, fs)
% analyses are based on Korn & Bach 2016 Journal of Vision
% PsPM toolbox: http://pspm.sourceforge.net/

addpath(genpath('C:\Jochem\repositories\toolboxes\PsPM_4.0.2\v4.0.2'))

getSubSpecs_bigDots
getfilenames_bigDots

%%% here we will model the pupil timecourse

%%% note on computing trial-by-trial estimates of coefficients (manual,
%%% p11) 
% 3.4.2 Multi-level modelling 
% Multi-level models, also termed linear mixed eects model (LME) or
% random-effect models, are becoming increasingly popular in psychology
% [11]. The idea is to take variance between trials into account: a summary
% statistic that is computed from individual data points with high variance
% is less precise that a statistic from low-variance data, and such
% information is lost when summarising within each condition. PsPM oers
% the possibility to compute trial-by-trial estimates (by modelling each
% trial as a separate condition in GLM, or by default in DCM for SCR). This
% approach has been used for SCR, PSR and SEBR [12], 3 GENERAL MODEL
% STRUCTURE 12 see also Homan et al. 2017 (Learning & Memory in press) and
% Tzovara et al. 2017 (in revision).




% MODEL with required fields
% model.modelfile:  a file name for the model output
% model.datafile:   a file name (single session) OR
%                   a cell array of file names
% model.timing:     a multiple condition file name (single session) OR
%                   a cell array of multiple condition file names OR
%                   a struct (single session) with fields .names, .onsets,
%                       and (optional) .durations and .pmod  OR
%                   a cell array of struct
% model.timeunits:  one of 'seconds', 'samples', 'markers'
% model.window:     a scalar in seconds that specifies over which time 
%                   window (starting with the events specified in
%                   model.timing) the model should be evaluated. Is only
%                   required if model.latency equals 'free'. Is ignored
%                   otherwise.

path_PsPM = [paths.s(subject.sub).base 'PsPM' filesep];
if ~exist(path_PsPM, 'dir')
    mkdir(path_PsPM)
end
model.modelfile = [path_PsPM allsubj{subject.sub} '_PsPM_GLM'];

% import settings

model.datafile = cell(1, length(allblocks{subject.sub}));
model.timing = cell(1, length(allblocks{subject.sub}));
icounter = 0;
for iblock = 1:length(allblocks{subject.sub})
    icounter = icounter+1;
    model.datafile{iblock} = files(subject.sub, allblocks{subject.sub}(iblock)).ET_files; % load ascii file
    % model.datafile{iblock} = fullfile([paths.data 'Samples_and_Events' filesep], files(1, allblocks{subject.sub}(iblock)).edf);
    
%     import{iblock}.eyelink_trackdist = 56; % 56 cm away from a 21 inch CRT (85 Hz, 1024 x 768 resolution)    
%     import{iblock}.type = 'pupil_l';
    

model.timing


end
    
% %     [data] = import_eyelink(model.datafile{iblock})
%     
%     [sts, import, sourceinfo] = pspm_get_eyelink(model.datafile, import);
%     
%     if strcmpi(sourceinfo.eyesObserved, 'r')
%         import{iblock}.type = 'pupil_r';
% %         [sts, import, sourceinfo] = pspm_get_eyelink(model.datafile{iblock}, import);
%     end
%     
%     options.overwrite = 1;
%     options.extrapolate = 0;
%     
%     % make sure we don't have to extrapolate, this gives weird results..
%     if isnan(import{iblock}.data(1))
%         import{iblock}.data(1) = nanmean(import{iblock}.data);
%     end
%     if isnan(import{iblock}.data(end))
%         import{iblock}.data(end) = nanmean(import{iblock}.data);
%     end
%     
    %     newdatafile = pspm_trim(import{iblock}, 'none', 'none', 'file', options)
%     [sts, import{iblock}.data] = pspm_interpolate(import{iblock}.data, options);
    
%     model.timing{iblock}.names = (blockTrialIdx(blockIdx==iblock));
% end
% 
% model.timing:     a multiple condition file name (single session) OR
%                   a cell array of multiple condition file names OR
%                   a struct (single session) with fields .names, .onsets,
%                       and (optional) .durations and .pmod  OR
%                   a cell array of struct
model.modelspec = 'psr';


pspm_glm(model)

% import.sr = 
%                       .sr
%                       .data


pspm_get_eyelink



[ntime, ntrial] = size(pupil);

% Analysis settings
prestim_range = [-200 0];  % pre-stimulus range (in ms) in which to measure average baseline pupil diameter

% Time settings
RTp_time = 2000;  % length of time (ms) after RT to pull on each trial

stim_times = [-200 3000];
resp_times = [-1500 2000];

% Constructing pupil impulse response function
w = 10.1;  % width of IRF (canonical = 10.1)
tmax = 930;  % time-to-peak of IRF (in ms) (canonical = 930)

pupil_IRF = [];
for t = 1:(fs*3)  % looping through 3 seconds of time
    t_ms = t/fs*1000;  % converting t to ms
    pupil_IRF(t,1) = (t_ms.^w).*(exp(-t_ms.*(w./tmax)));  % contructing IRF
end
%pupil_IRF = pupil_IRF./max(pupil_IRF).*0.4;
pupil_IRF = pupil_IRF.*1e-26;  % scaling IRF so that betas are on a manageable scale


% baseline pupil diameter
baseline    = mean(pupil(find(tt>=prestim_range(1) & tt<=prestim_range(2)),:),1);
pupil       = pupil  - repmat(baseline,[size(pupil,1)],1); % baseline full erp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PULLING PUPIL SEGMENTS AND CONSTRUCTING REGRESSORS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% z-scoring RTs
zRT = zscore(RT);

Pupil_conc=[];
Stim_reg=[]; Boxc_reg=[]; Boxc2thresh_reg=[]; Ramp_reg=[]; Ramp2thresh_reg=[]; Resp_reg=[];
Stim_reg_P=[]; Boxc_reg_P=[]; Boxc2thresh_reg_P=[]; Ramp_reg_P=[]; Ramp2thresh_reg_P=[]; Resp_reg_P=[];

trialLength = NaN(ntrial, 1);
for itrial = 1:ntrial
    % Pulling pupil epoch for this trial
    current_pupil = pupil(find( tt>=prestim_range(1) & tt<=(RT(itrial) + RTp_time) ),itrial);
    Pupil_conc(end+1:end+length(current_pupil),1) = current_pupil;
    
    trialLength(itrial) = length(current_pupil);
    
    % Creating STIMULUS regressor segments for this trial
    temp_reg = zeros(size(current_pupil));
    temp_reg(((-prestim_range(1)/1000)*fs)+1) = 1;
    Stim_reg(end+1:end+length(temp_reg),1) = temp_reg;
    Stim_reg_P(end+1:end+length(temp_reg),1) = temp_reg.*zRT(itrial);
    
    % Creating RESPONSE regressor segments for this trial
    temp_reg = zeros(size(current_pupil));
    temp_reg(ceil(((RT(itrial))-prestim_range(1))/1000*fs)) = 1;
    Resp_reg(end+1:end+length(temp_reg),1) = temp_reg;
    Resp_reg_P(end+1:end+length(temp_reg),1) = temp_reg.*zRT(itrial);
    
    % Creating DECISION regressor segments for this trial
    int_leng = length((((-prestim_range(1)/1000)*fs)+2):(floor(((RT(itrial))-prestim_range(1))/1000*fs)));  % getting length of decision interval
    boxc = ones(1,int_leng);                % plain boxcar regressor
    ramp = 1:int_leng;                      % pure ramp
    box2thresh = ones(1,int_leng)/int_leng; % plain boxcar regressor divided by length of trial (RT)
    ramp2thresh = (1:int_leng)./int_leng;   % ramp-to-threshold
    
    temp_box = zeros(size(current_pupil)); temp_box((((-prestim_range(1)/1000)*fs)+2):(floor(((RT(itrial))-prestim_range(1))/1000*fs))) = boxc;
    Boxc_reg(end+1:end+length(temp_box),1) = temp_box;
    Boxc_reg_P(end+1:end+length(temp_box),1) = temp_box.*zRT(itrial);
    
    temp_box = zeros(size(current_pupil)); temp_box((((-prestim_range(1)/1000)*fs)+2):(floor((RT(itrial)-prestim_range(1))/1000*fs))) = box2thresh;
    Boxc2thresh_reg(end+1:end+length(temp_box),1) = temp_box;
    Boxc2thresh_reg_P(end+1:end+length(temp_box),1) = temp_box.*zRT(itrial);

    temp_ramp = zeros(size(current_pupil)); temp_ramp((((-prestim_range(1)/1000)*fs)+2):(floor(((RT(itrial))-prestim_range(1))/1000*fs))) = ramp;
    Ramp_reg(end+1:end+length(temp_ramp),1) = temp_ramp;
    Ramp_reg_P(end+1:end+length(temp_ramp),1) = temp_ramp.*zRT(itrial);
    
    temp_ramp2thresh = zeros(size(current_pupil)); temp_ramp2thresh((((-prestim_range(1)/1000)*fs)+2):(floor(((RT(itrial))-prestim_range(1))/1000*fs))) = ramp2thresh;
    Ramp2thresh_reg(end+1:end+length(temp_ramp2thresh),1) = temp_ramp2thresh;
    Ramp2thresh_reg_P(end+1:end+length(temp_ramp2thresh),1) = temp_ramp2thresh.*zRT(itrial);
end


temp_conv=conv(Stim_reg,pupil_IRF); Stim_reg=temp_conv(1:length(Stim_reg));
temp_conv=conv(Resp_reg,pupil_IRF); Resp_reg=temp_conv(1:length(Resp_reg));
temp_conv=conv(Boxc_reg,pupil_IRF); Boxc_reg=temp_conv(1:length(Boxc_reg));
temp_conv=conv(Boxc2thresh_reg,pupil_IRF); Boxc2thresh_reg=temp_conv(1:length(Boxc2thresh_reg));
temp_conv=conv(Ramp_reg,pupil_IRF); Ramp_reg=temp_conv(1:length(Ramp_reg));
temp_conv=conv(Ramp2thresh_reg,pupil_IRF); Ramp2thresh_reg=temp_conv(1:length(Ramp2thresh_reg));

temp_conv=conv(Stim_reg_P,pupil_IRF); Stim_reg_P=temp_conv(1:length(Stim_reg_P));
temp_conv=conv(Resp_reg_P,pupil_IRF); Resp_reg_P=temp_conv(1:length(Resp_reg_P));
temp_conv=conv(Boxc_reg_P,pupil_IRF); Boxc_reg_P=temp_conv(1:length(Boxc_reg_P));
temp_conv=conv(Ramp_reg_P,pupil_IRF); Ramp_reg_P=temp_conv(1:length(Ramp_reg_P));
temp_conv=conv(Ramp2thresh_reg_P,pupil_IRF); Ramp2thresh_reg_P=temp_conv(1:length(Ramp2thresh_reg_P));


% FR three-param models, no PM
n = length(Pupil_conc); % number of observations
% 

%%% from old script
GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg],'linear',{'yhat','r','beta','tstat'});
GLM2 = regstats(Pupil_conc,[Stim_reg Ramp_reg Resp_reg],'linear',{'yhat','r','beta','tstat'});
GLM3 = regstats(Pupil_conc,[Stim_reg Ramp2thresh_reg Resp_reg],'linear',{'yhat','r','beta','tstat'});
GLM4 = regstats(Pupil_conc,[Stim_reg Boxc2thresh_reg Resp_reg],'linear',{'yhat','r','beta','tstat'});


% GLM1 = regstats(Pupil_conc_FR,[Stim_reg_FR Boxc_reg_FR Resp_reg_FR],'linear',{'yhat','beta','r'});
% GLM2 = regstats(Pupil_conc_FR,[Stim_reg_FR Ramp_reg_FR Resp_reg_FR],'linear',{'yhat','beta','r'});
% GLM3 = regstats(Pupil_conc_FR,[Stim_reg_FR Ramp2thresh_reg_FR Resp_reg_FR],'linear',{'yhat','beta','r'});

% GA_BIC_FR(1) = n+(n.*log(2*pi))+(n.*log(sum(GLM1.r.^2)./n))+(log(n).*(4+1));
% GA_BIC_FR(2) = n+(n.*log(2*pi))+(n.*log(sum(GLM2.r.^2)./n))+(log(n).*(4+1));
% GA_BIC_FR(3) = n+(n.*log(2*pi))+(n.*log(sum(GLM3.r.^2)./n))+(log(n).*(4+1));
% 
idx = 0;
pupil = NaN(max(trialLength), ntrial);
pupilFit = NaN(max(trialLength), ntrial);
for itrial = 1:ntrial
    current_pupil = Pupil_conc((idx+1):sum(trialLength(1:itrial)));
    current_pupil_fit = GLM4.yhat((idx+1):sum(trialLength(1:itrial)));
    
    pupil(1:trialLength(itrial), itrial) = current_pupil;
    pupilFit(1:trialLength(itrial), itrial) = current_pupil_fit;
    idx = sum(trialLength(1:itrial));
    
end


    
figure(2),clf
hold on
plot(nanmean(pupil, 2))
plot(nanmean(pupilFit, 2))



%%% from old script
GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg],'linear',{'yhat','r','beta','tstat'});
GLM2 = regstats(Pupil_conc,[Stim_reg Ramp_reg Resp_reg],'linear',{'yhat','r','beta','tstat'});
GLM3 = regstats(Pupil_conc,[Stim_reg Ramp2thresh_reg Resp_reg],'linear',{'yhat','r','beta','tstat'});
GLM4 = regstats(Pupil_conc,[Stim_reg Boxc2thresh_reg Resp_reg],'linear',{'yhat','r','beta','tstat'});

figure(1),clf
subplot(5,1,1)
hold on
plot(Stim_reg)
plot(Stim_reg_P)
subplot(5,1,2)
hold on
plot(Resp_reg)
plot(Resp_reg_P)

subplot(5,1,3)
hold on
plot(Boxc_reg)
plot(Boxc_reg_P)

subplot(5,1,4)
hold on
plot(Ramp_reg)
plot(Ramp_reg_P)

subplot(5,1,5)
hold on
plot(Ramp2thresh_reg)
plot(Ramp2thresh_reg_P)

subplot(3,1,1)
plot([Pupil_conc Stim_reg])
subplot(3,1,2)
plot([Pupil_conc Boxc_reg])
subplot(3,1,3)
plot([Pupil_conc Resp_reg])


figure(1),clf
hold on
plot(Pupil_conc, 'linewidth',2)
plot([Stim_reg, Boxc_reg, Resp_reg], 'linewidth',1.5)
plot([GLM1.yhat], '--', 'linewidth',1.5)
plot([GLM1.r], '-.', 'linewidth',1.5)
plot([GLM1.r.^2], '-.', 'linewidth',1.5)
h = legend({'data','stim_reg','box_reg','resp_reg','box_yhat','box_resid','box_resid^2'});
set(h,'Interpreter','none')




% FR three-param models w/ Stim PM
GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg Stim_reg_P],'linear',{'r'});
GLM2 = regstats(Pupil_conc,[Stim_reg Ramp_reg Resp_reg Stim_reg_P],'linear',{'r'});
GLM3 = regstats(Pupil_conc,[Stim_reg Ramp2thresh_reg Resp_reg Stim_reg_P],'linear',{'r'});

GA_BIC_FR_StimPM(1) = n+(n.*log(2*pi))+(n.*log(sum(GLM1.r.^2)./n))+(log(n).*(5+1));
GA_BIC_FR_StimPM(2) = n+(n.*log(2*pi))+(n.*log(sum(GLM2.r.^2)./n))+(log(n).*(5+1));
GA_BIC_FR_StimPM(3) = n+(n.*log(2*pi))+(n.*log(sum(GLM3.r.^2)./n))+(log(n).*(5+1));

% FR three-param models w/ Dec PM
GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg Boxc_reg_P],'linear',{'r'});
GLM2 = regstats(Pupil_conc,[Stim_reg Ramp_reg Resp_reg Ramp_reg_P],'linear',{'r'});
GLM3 = regstats(Pupil_conc,[Stim_reg Ramp2thresh_reg Resp_reg Ramp2thresh_reg_P],'linear',{'r'});

GA_BIC_FR_DecPM(1) = n+(n.*log(2*pi))+(n.*log(sum(GLM1.r.^2)./n))+(log(n).*(5+1));
GA_BIC_FR_DecPM(2) = n+(n.*log(2*pi))+(n.*log(sum(GLM2.r.^2)./n))+(log(n).*(5+1));
GA_BIC_FR_DecPM(3) = n+(n.*log(2*pi))+(n.*log(sum(GLM3.r.^2)./n))+(log(n).*(5+1));

% FR three-param models w/ Resp PM
GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg Resp_reg_P],'linear',{'r'});
GLM2 = regstats(Pupil_conc,[Stim_reg Ramp_reg Resp_reg Resp_reg_P],'linear',{'r'});
GLM3 = regstats(Pupil_conc,[Stim_reg Ramp2thresh_reg Resp_reg Resp_reg_P],'linear',{'r'});

GA_BIC_FR_RespPM(1) = n+(n.*log(2*pi))+(n.*log(sum(GLM1.r.^2)./n))+(log(n).*(5+1));
GA_BIC_FR_RespPM(2) = n+(n.*log(2*pi))+(n.*log(sum(GLM2.r.^2)./n))+(log(n).*(5+1));
GA_BIC_FR_RespPM(3) = n+(n.*log(2*pi))+(n.*log(sum(GLM3.r.^2)./n))+(log(n).*(5+1));

% FR three-param models w/ Stim & Dec PM
GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg Stim_reg_P Boxc_reg_P],'linear',{'r'});
GLM2 = regstats(Pupil_conc,[Stim_reg Ramp_reg Resp_reg Stim_reg_P Ramp_reg_P],'linear',{'r'});
GLM3 = regstats(Pupil_conc,[Stim_reg Ramp2thresh_reg Resp_reg Stim_reg_P Ramp2thresh_reg_P],'linear',{'r'});

GA_BIC_FR_StimDecPM(1) = n+(n.*log(2*pi))+(n.*log(sum(GLM1.r.^2)./n))+(log(n).*(6+1));
GA_BIC_FR_StimDecPM(2) = n+(n.*log(2*pi))+(n.*log(sum(GLM2.r.^2)./n))+(log(n).*(6+1));
GA_BIC_FR_StimDecPM(3) = n+(n.*log(2*pi))+(n.*log(sum(GLM3.r.^2)./n))+(log(n).*(6+1));

% FR three-param models w/ Dec & Resp PM
GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg Boxc_reg_P Resp_reg_P],'linear',{'r'});
GLM2 = regstats(Pupil_conc,[Stim_reg Ramp_reg Resp_reg Ramp_reg_P Resp_reg_P],'linear',{'r'});
GLM3 = regstats(Pupil_conc,[Stim_reg Ramp2thresh_reg Resp_reg Ramp2thresh_reg_P Resp_reg_P],'linear',{'r'});

GA_BIC_FR_DecRespPM(1) = n+(n.*log(2*pi))+(n.*log(sum(GLM1.r.^2)./n))+(log(n).*(6+1));
GA_BIC_FR_DecRespPM(2) = n+(n.*log(2*pi))+(n.*log(sum(GLM2.r.^2)./n))+(log(n).*(6+1));
GA_BIC_FR_DecRespPM(3) = n+(n.*log(2*pi))+(n.*log(sum(GLM3.r.^2)./n))+(log(n).*(6+1));

% FR three-param models w/ Stim & Resp PM
GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg Stim_reg_P Resp_reg_P],'linear',{'r'});
GLM2 = regstats(Pupil_conc,[Stim_reg Ramp_reg Resp_reg Stim_reg_P Resp_reg_P],'linear',{'r'});
GLM3 = regstats(Pupil_conc,[Stim_reg Ramp2thresh_reg Resp_reg Stim_reg_P Resp_reg_P],'linear',{'r'});

GA_BIC_FR_StimRespPM(1) = n+(n.*log(2*pi))+(n.*log(sum(GLM1.r.^2)./n))+(log(n).*(6+1));
GA_BIC_FR_StimRespPM(2) = n+(n.*log(2*pi))+(n.*log(sum(GLM2.r.^2)./n))+(log(n).*(6+1));
GA_BIC_FR_StimRespPM(3) = n+(n.*log(2*pi))+(n.*log(sum(GLM3.r.^2)./n))+(log(n).*(6+1));

% FR three-param models w/ ALL PM
GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg Stim_reg_P Boxc_reg_P Resp_reg_P],'linear',{'r','rsquare'});
GLM2 = regstats(Pupil_conc,[Stim_reg Ramp_reg Resp_reg Stim_reg_P Ramp_reg_P Resp_reg_P],'linear',{'r'});
GLM3 = regstats(Pupil_conc,[Stim_reg Ramp2thresh_reg Resp_reg Stim_reg_P Ramp2thresh_reg_P Resp_reg_P],'linear',{'r'});

GA_BIC_FR_AllPM(1) = n+(n.*log(2*pi))+(n.*log(sum(GLM1.r.^2)./n))+(log(n).*(7+1));
GA_BIC_FR_AllPM(2) = n+(n.*log(2*pi))+(n.*log(sum(GLM2.r.^2)./n))+(log(n).*(7+1));
GA_BIC_FR_AllPM(3) = n+(n.*log(2*pi))+(n.*log(sum(GLM3.r.^2)./n))+(log(n).*(7+1));

GA_allPM_rsquares(2) = GLM1.rsquare;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUNNING BEST-FITTING 3-PARAM GLMs W/ AND W/OUT PM AND PULLING BETAS%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_fullPM = [1 0 0 0 0 0;...  % Specifying model - include [0 0 0 0 0 0] as first line if want to include constant
    0 1 0 0 0 0;...
    0 0 1 0 0 0;...
    0 0 0 1 0 0;...
    0 0 0 0 1 0;...
    0 0 0 0 0 1];

GLM1 = regstats(Pupil_conc_DL,[Stim_reg_DL Ramp_reg_DL Resp_reg_DL Stim_reg_DL_P Ramp_reg_DL_P Resp_reg_DL_P],model_fullPM,{'beta','tstat'});
if sum(model_fullPM(1,:))==0
    GA_DL_main_ts(1:6) = GLM1.tstat.t(2:7);  % storing t-scores
else GA_DL_main_ts(1:6) = GLM1.tstat.t(1:6);
end

GLM1 = regstats(Pupil_conc,[Stim_reg Boxc_reg Resp_reg Stim_reg_P Boxc_reg_P Resp_reg_P],model_fullPM,{'beta','tstat'});
if sum(model_fullPM(1,:))==0
    GA_FR_main_ts(1:6) = GLM1.tstat.t(2:7);  % storing t-scores
else GA_FR_main_ts(1:6) = GLM1.tstat.t(1:6);
end








figure(1),clf
hold on
plot(Pupil_conc)
plot(Stim_reg)
plot(Resp_reg)
plot(Boxc_reg)
% plot(Ramp_reg_FR)
% plot(Ramp2thresh_reg_FR)


keyboard



