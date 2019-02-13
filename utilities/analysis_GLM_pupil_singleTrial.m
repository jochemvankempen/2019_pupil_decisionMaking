function [st_GLMfit, conc_GLMfit, bin_GLMfit] = analysis_GLM_pupil_singleTrial(pupil, cfg, RT, validtr2use, tt, fs)
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
[~, ntrial] = size(pupil);

% Analysis settings
prestim_range = cfg.prestim_range;  % pre-stimulus range (in ms) in which to measure average baseline pupil diameter

% Time settings
RTp_time = 2500;  % length of time (ms) after RT to pull on each trial

% Constructing pupil impulse response function
w = 10.1;  % width of IRF (canonical = 10.1)
tmax = 930;  % time-to-peak of IRF (in ms) (canonical = 930)

pupil_IRF = [];
for t = 1:(fs*3)  % looping through 3 seconds of time
    t_ms(t) = t/fs*1000;  % converting t to ms
    pupil_IRF(t,1) = (t_ms(t).^w).*(exp(-t_ms(t).*(w./tmax)));  % contructing IRF
end
%pupil_IRF = pupil_IRF./max(pupil_IRF).*0.4;
pupil_IRF = pupil_IRF.*1e-26;  % scaling IRF so that betas are on a manageable scale


% baseline pupil diameter
baseline    = mean(pupil(find(tt>=prestim_range(1) & tt<=prestim_range(2)),:),1);
pupil       = pupil  - repmat(baseline,[size(pupil,1)],1); % baseline full erp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PULLING PUPIL SEGMENTS AND CONSTRUCTING REGRESSORS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st_pupil        = cell(ntrial, 1);
trialTimeIdx    = cell(ntrial, 1);
trialTime       = cell(ntrial, 1);
    
Stim_reg            = cell(ntrial, 1);
Resp_reg            = cell(ntrial, 1);
Boxc_reg            = cell(ntrial, 1);
Ramp_reg            = cell(ntrial, 1);
Ramp2thresh_reg     = cell(ntrial, 1);
Rampflip_reg        = cell(ntrial, 1);
Ramp2threshflip_reg = cell(ntrial, 1);
BoxcScaledDown_reg  = cell(ntrial, 1);
RampScaled_reg      = cell(ntrial, 1);
RampflipScaled_reg  = cell(ntrial, 1);

trialLength = NaN(ntrial, 1);

% initialise
% murphy et al., 2016 components
conc_Pupil=[];
conc_Stim_reg=[]; 
conc_Boxc_reg=[]; conc_Ramp_reg=[]; conc_Ramp2thresh_reg=[]; conc_Rampflip_reg=[]; conc_Ramp2threshflip_reg=[]; 
conc_BoxcScaledDown_reg=[]; conc_RampScaled_reg=[]; conc_RampflipScaled_reg=[];
conc_Resp_reg=[];

% single trial structures
st_GLMfit.VIF.StimResp = NaN(ntrial, 2);
st_GLMfit.VIF.Boxc = NaN(ntrial, 3);
st_GLMfit.VIF.Ramp = NaN(ntrial, 3);

for itrial = 1:ntrial
    
    st_GLMfit.data(itrial).Pupil = [];
    st_GLMfit.data(itrial).TimeIdx = [];
    
    st_GLMfit.fit(itrial).StimResp = struct([]);
    st_GLMfit.fit(itrial).Boxc = struct([]);
    st_GLMfit.fit(itrial).Ramp = struct([]);
    if ~validtr2use(itrial)
        continue
    end
    
    % Pulling pupil epoch for this trial
    trialTimeIdx{itrial} = tt>=prestim_range(1) & tt<=(RT(itrial) + RTp_time);
    current_pupil       = pupil(trialTimeIdx{itrial},itrial);
    trialTime{itrial}   = tt(trialTimeIdx{itrial});
    trialLength(itrial) = length(current_pupil);
    
    st_pupil{itrial} = current_pupil;
    conc_Pupil(end+1:end+length(current_pupil),1) = current_pupil;
    
    % Creating STIMULUS regressor segments for this trial
    temp_reg = zeros(size(current_pupil));
    temp_reg(((-prestim_range(1)/1000)*fs)+1) = 1;
    Stim_reg{itrial}   = temp_reg;
    conc_Stim_reg(end+1:end+length(temp_reg),1) = temp_reg;
    
    % Creating RESPONSE regressor segments for this trial
    temp_reg = zeros(size(current_pupil));
    temp_reg(ceil((RT(itrial)-prestim_range(1))/1000*fs)) = 1;
    Resp_reg{itrial} = temp_reg;
    conc_Resp_reg(end+1:end+length(temp_reg),1) = temp_reg;
    
    % Creating DECISION regressor segments for this trial
    int_idx =  (((-prestim_range(1)/1000)*fs)+2):(floor(((RT(itrial))-prestim_range(1))/1000*fs));
    int_leng = length(int_idx);  % getting length of decision interval

    boxc = ones(1,int_leng);                        % plain boxcar regressor. Model 1 in Murphy et al 2016 Nat Comm
    ramp = 1:int_leng;                              % pure ramp. Model 2
    ramp2thresh = (1:int_leng)./int_leng;           % ramp-to-threshold. Model 3
    ramp_flip = flip(1:int_leng);                   % pure ramp, reversed. Model 4
    ramp2thresh_flip = flip((1:int_leng)./int_leng);% ramp-to-threshold. reversed. Model 5
    boxc_scaledDown = ones(1,int_leng)./int_leng;   % boxcar regressor, scaled down with duration trial. Model 6
    rampScaled = (1:int_leng)./int_leng^2;          % pure ramp. scaled down. Model 7
    ramp_flip_Scaled = flip(1:int_leng)./int_leng^2;% pure ramp. reversed. scaled down. Model 8
    
    % 1. Box
    temp_box = zeros(size(current_pupil)); temp_box(int_idx) = boxc;
    Boxc_reg{itrial} = temp_box;

    conc_Boxc_reg(end+1:end+length(temp_box),1) = temp_box;
        
    % 2. Ramp
    temp_ramp = zeros(size(current_pupil)); temp_ramp(int_idx) = ramp;
    Ramp_reg{itrial} = temp_ramp;

    conc_Ramp_reg(end+1:end+length(temp_ramp),1) = temp_ramp;
        
    % 3. Ramp2thresh
    temp_ramp2thresh = zeros(size(current_pupil)); temp_ramp2thresh(int_idx) = ramp2thresh;
    Ramp2thresh_reg{itrial} = temp_ramp2thresh;
    
    conc_Ramp2thresh_reg(end+1:end+length(temp_ramp2thresh),1) = temp_ramp2thresh;
    
    % 4. Ramp reversed
    temp_rampflip = zeros(size(current_pupil)); temp_rampflip(int_idx) = ramp_flip;
    Rampflip_reg{itrial} = temp_rampflip;
    
    conc_Rampflip_reg(end+1:end+length(temp_rampflip),1) = temp_rampflip;

     % 5. Ramp2thresh reversed
    temp_ramp2threshflip = zeros(size(current_pupil)); temp_ramp2threshflip(int_idx) = ramp2thresh_flip;
    Ramp2threshflip_reg{itrial} = temp_ramp2threshflip;
    
    conc_Ramp2threshflip_reg(end+1:end+length(temp_ramp2threshflip),1) = temp_ramp2threshflip;
   
    % 6. Box. scaled down
    temp_box_scaledDown = zeros(size(current_pupil)); temp_box_scaledDown(int_idx) = boxc_scaledDown;
    BoxcScaledDown_reg{itrial} = temp_box_scaledDown;
    
    conc_BoxcScaledDown_reg(end+1:end+length(temp_box_scaledDown),1) = temp_box_scaledDown;

    % 7. Ramp. scaled down
    temp_ramp = zeros(size(current_pupil)); temp_ramp(int_idx) = rampScaled;
    RampScaled_reg{itrial} = temp_ramp;
    
    conc_RampScaled_reg(end+1:end+length(temp_ramp),1) = temp_ramp;

    % 8. Ramp. reversed. scaled down
    temp_ramp = zeros(size(current_pupil)); temp_ramp(int_idx) = ramp_flip_Scaled;
    RampflipScaled_reg{itrial} = temp_ramp;
    
    conc_RampflipScaled_reg(end+1:end+length(temp_ramp),1) = temp_ramp;

    if 0
        % check timings
        figure(1), clf
        hold on
%         plot(tt(trialTimeIdx{itrial}), current_pupil)
        plot(tt(trialTimeIdx{itrial}), Stim_reg{itrial})
        plot(tt(trialTimeIdx{itrial}), Resp_reg{itrial})
        plot(tt(trialTimeIdx{itrial}), Boxc_reg{itrial})
        plot(tt(trialTimeIdx{itrial}), Ramp_reg{itrial})
        plot(tt(trialTimeIdx{itrial}), Ramp2thresh_reg{itrial})
        plot(tt(trialTimeIdx{itrial}), Rampflip_reg{itrial})
        plot(tt(trialTimeIdx{itrial}), Ramp2threshflip_reg{itrial})
        plot(tt(trialTimeIdx{itrial}), BoxcScaledDown_reg{itrial})
        
    end
    
    % convolution
    temp_conv=conv(Stim_reg{itrial}, pupil_IRF); Stim_reg{itrial}=temp_conv(1:length(Stim_reg{itrial}));
    temp_conv=conv(Resp_reg{itrial}, pupil_IRF); Resp_reg{itrial}=temp_conv(1:length(Resp_reg{itrial}));
    temp_conv=conv(Boxc_reg{itrial}, pupil_IRF); Boxc_reg{itrial}=temp_conv(1:length(Boxc_reg{itrial}));
    temp_conv=conv(Ramp_reg{itrial}, pupil_IRF); Ramp_reg{itrial}=temp_conv(1:length(Ramp_reg{itrial}));
    temp_conv=conv(Ramp2thresh_reg{itrial}, pupil_IRF); Ramp2thresh_reg{itrial}=temp_conv(1:length(Ramp2thresh_reg{itrial}));
    temp_conv=conv(Rampflip_reg{itrial}, pupil_IRF); Rampflip_reg{itrial}=temp_conv(1:length(Rampflip_reg{itrial}));
    temp_conv=conv(Ramp2threshflip_reg{itrial}, pupil_IRF); Ramp2threshflip_reg{itrial}=temp_conv(1:length(Ramp2threshflip_reg{itrial}));
    temp_conv=conv(BoxcScaledDown_reg{itrial}, pupil_IRF); BoxcScaledDown_reg{itrial}=temp_conv(1:length(BoxcScaledDown_reg{itrial}));
    temp_conv=conv(RampScaled_reg{itrial}, pupil_IRF); RampScaled_reg{itrial}=temp_conv(1:length(RampScaled_reg{itrial}));
    temp_conv=conv(RampflipScaled_reg{itrial}, pupil_IRF); RampflipScaled_reg{itrial}=temp_conv(1:length(RampflipScaled_reg{itrial}));
    
    if cfg.meanCenter_predictors
        % mean center
        Stim_reg{itrial} = Stim_reg{itrial} - mean(Stim_reg{itrial});
        Resp_reg{itrial} = Resp_reg{itrial} - mean(Resp_reg{itrial});
        Boxc_reg{itrial} = Boxc_reg{itrial} - mean(Boxc_reg{itrial});
        Ramp_reg{itrial} = Ramp_reg{itrial} - mean(Ramp_reg{itrial});
        Ramp2thresh_reg{itrial} = Ramp2thresh_reg{itrial} - mean(Ramp2thresh_reg{itrial});
        Rampflip_reg{itrial} = Rampflip_reg{itrial} - mean(Rampflip_reg{itrial});
        Ramp2threshflip_reg{itrial} = Ramp2threshflip_reg{itrial} - mean(Ramp2threshflip_reg{itrial});
        BoxcScaledDown_reg{itrial} = BoxcScaledDown_reg{itrial} - mean(BoxcScaledDown_reg{itrial});
        RampScaled_reg{itrial} = RampScaled_reg{itrial} - mean(RampScaled_reg{itrial});
        RampflipScaled_reg{itrial} = RampflipScaled_reg{itrial} - mean(RampflipScaled_reg{itrial});
        
    end
    
    if cfg.normalise_pupil == 1
        % normalise
        tmp_pupil = (st_pupil{itrial} - mean(st_pupil{itrial}))/std(st_pupil{itrial});
    elseif cfg.normalise_pupil == 2
        tmp_pupil = (st_pupil{itrial} - mean(st_pupil{itrial}));
    else
        tmp_pupil = st_pupil{itrial};
    end

    % single trial regression
    st_GLMfit.data(itrial).Pupil = tmp_pupil;
    st_GLMfit.data(itrial).TimeIdx = trialTimeIdx{itrial};

    if cfg.orthogonalise_predictors
        st_GLMfit.fit(itrial).StimResp = regstats(tmp_pupil,spm_orth([Stim_reg{itrial} Resp_reg{itrial}]),'linear',{'yhat','r','beta','tstat','rsquare'});
        st_GLMfit.fit(itrial).Boxc = regstats(tmp_pupil,spm_orth([Stim_reg{itrial} Boxc_reg{itrial} Resp_reg{itrial}]),'linear',{'yhat','r','beta','tstat','rsquare'});
        st_GLMfit.fit(itrial).Ramp = regstats(tmp_pupil,spm_orth([Stim_reg{itrial} Ramp_reg{itrial} Resp_reg{itrial}]),'linear',{'yhat','r','beta','tstat','rsquare'});
    else
        st_GLMfit.fit(itrial).StimResp = regstats(tmp_pupil,[Stim_reg{itrial} Resp_reg{itrial}],'linear',{'yhat','r','beta','tstat','rsquare'});
        st_GLMfit.fit(itrial).Boxc = regstats(tmp_pupil,[Stim_reg{itrial} Boxc_reg{itrial} Resp_reg{itrial}],'linear',{'yhat','r','beta','tstat','rsquare'});
        st_GLMfit.fit(itrial).Ramp = regstats(tmp_pupil,[Stim_reg{itrial} Ramp_reg{itrial} Resp_reg{itrial}],'linear',{'yhat','r','beta','tstat','rsquare'});
    end
    
    % Variance inflation factor
    if cfg.orthogonalise_predictors
        st_GLMfit.VIF.StimResp(itrial,:) = varianceInflationCoefficients(spm_orth([Stim_reg{itrial} Resp_reg{itrial}]));
        st_GLMfit.VIF.Boxc(itrial,:) = varianceInflationCoefficients(spm_orth([Stim_reg{itrial} Boxc_reg{itrial} Resp_reg{itrial}]));
        st_GLMfit.VIF.Ramp(itrial,:) = varianceInflationCoefficients(spm_orth([Stim_reg{itrial} Ramp_reg{itrial} Resp_reg{itrial}]));
    else
        st_GLMfit.VIF.StimResp(itrial,:) = varianceInflationCoefficients([Stim_reg{itrial} Resp_reg{itrial}]);
        st_GLMfit.VIF.Boxc(itrial,:) = varianceInflationCoefficients([Stim_reg{itrial} Boxc_reg{itrial} Resp_reg{itrial}]);
        st_GLMfit.VIF.Ramp(itrial,:) = varianceInflationCoefficients([Stim_reg{itrial} Ramp_reg{itrial} Resp_reg{itrial}]);
    end
end

if 0
    %% check regression shape
    figure(2), clf
    hold on
    
    RT2test = [-1.5 0 2]
    
    for itrial = 1:length(RTest)
        
        [~, idx] = min(abs(RT2test(itrial) - zRT));
%         plot(tt(trialTimeIdx{idx}), Boxc_reg{idx})
%         plot(tt(trialTimeIdx{idx}), Ramp_reg{idx})
%         plot(tt(trialTimeIdx{idx}), Ramp2thresh_reg{idx})
%         plot(tt(trialTimeIdx{idx}), Rampflip_reg{idx})
%         plot(tt(trialTimeIdx{idx}), Ramp2threshflip_reg{idx})
%         plot(tt(trialTimeIdx{idx}), BoxcScaledDown_reg{idx})
        plot(tt(trialTimeIdx{idx}), RampScaled_reg{idx})
%         plot(tt(trialTimeIdx{idx}), RampflipScaled_reg{idx})
    end
end

% across trial convolution
temp_conv=conv(conc_Stim_reg,pupil_IRF); conc_Stim_reg=temp_conv(1:length(conc_Stim_reg));
temp_conv=conv(conc_Resp_reg,pupil_IRF); conc_Resp_reg=temp_conv(1:length(conc_Resp_reg));
temp_conv=conv(conc_Boxc_reg,pupil_IRF); conc_Boxc_reg=temp_conv(1:length(conc_Boxc_reg));
temp_conv=conv(conc_Ramp_reg,pupil_IRF); conc_Ramp_reg=temp_conv(1:length(conc_Ramp_reg));
temp_conv=conv(conc_Ramp2thresh_reg,pupil_IRF); conc_Ramp2thresh_reg=temp_conv(1:length(conc_Ramp2thresh_reg));
temp_conv=conv(conc_Rampflip_reg,pupil_IRF); conc_Rampflip_reg=temp_conv(1:length(conc_Rampflip_reg));
temp_conv=conv(conc_Ramp2threshflip_reg,pupil_IRF); conc_Ramp2threshflip_reg=temp_conv(1:length(conc_Ramp2threshflip_reg));
temp_conv=conv(conc_BoxcScaledDown_reg,pupil_IRF); conc_BoxcScaledDown_reg=temp_conv(1:length(conc_BoxcScaledDown_reg));
temp_conv=conv(conc_RampScaled_reg,pupil_IRF); conc_RampScaled_reg=temp_conv(1:length(conc_RampScaled_reg));
temp_conv=conv(conc_RampflipScaled_reg,pupil_IRF); conc_RampflipScaled_reg=temp_conv(1:length(conc_RampflipScaled_reg));

if cfg.meanCenter_predictors
    % mean center
    conc_Stim_reg = conc_Stim_reg - mean(conc_Stim_reg);
    conc_Resp_reg = conc_Resp_reg - mean(conc_Resp_reg);
    conc_Boxc_reg = conc_Boxc_reg - mean(conc_Boxc_reg);
    conc_Ramp_reg = conc_Ramp_reg - mean(conc_Ramp_reg);
    conc_Ramp2thresh_reg = conc_Ramp2thresh_reg - mean(conc_Ramp2thresh_reg);
    conc_Rampflip_reg = conc_Rampflip_reg - mean(conc_Rampflip_reg);
    conc_Ramp2threshflip_reg = conc_Ramp2threshflip_reg - mean(conc_Ramp2threshflip_reg);
    conc_BoxcScaledDown_reg = conc_BoxcScaledDown_reg - mean(conc_BoxcScaledDown_reg);
    conc_RampScaled_reg = conc_RampScaled_reg - mean(conc_RampScaled_reg);
    conc_RampflipScaled_reg = conc_RampflipScaled_reg - mean(conc_RampflipScaled_reg);
    
end

if cfg.normalise_pupil == 1
    % normalise
    conc_Pupil = (conc_Pupil - mean(conc_Pupil))/std(conc_Pupil);
elseif cfg.normalise_pupil == 2
    conc_Pupil = (conc_Pupil - mean(conc_Pupil));
end

%%% across trial regression
n = length(conc_Pupil); % number of observations

% three-param models
if cfg.orthogonalise_predictors
    conc_GLMfit.stimResp        = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Boxc            = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_Boxc_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Ramp            = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_Ramp_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Ramp2thresh     = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_Ramp2thresh_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Rampflip        = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_Rampflip_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Ramp2threshflip = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_Ramp2threshflip_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.BoxcScaledDown  = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_BoxcScaledDown_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.RampScaled      = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_RampScaled_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.RampflipScaled  = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_RampflipScaled_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
else
    conc_GLMfit.stimResp        = regstats(conc_Pupil,[conc_Stim_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Boxc            = regstats(conc_Pupil,[conc_Stim_reg conc_Boxc_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Ramp            = regstats(conc_Pupil,[conc_Stim_reg conc_Ramp_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Ramp2thresh     = regstats(conc_Pupil,[conc_Stim_reg conc_Ramp2thresh_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Rampflip        = regstats(conc_Pupil,[conc_Stim_reg conc_Rampflip_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.Ramp2threshflip = regstats(conc_Pupil,[conc_Stim_reg conc_Ramp2threshflip_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.BoxcScaledDown  = regstats(conc_Pupil,[conc_Stim_reg conc_BoxcScaledDown_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.RampScaled      = regstats(conc_Pupil,[conc_Stim_reg conc_RampScaled_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
    conc_GLMfit.RampflipScaled  = regstats(conc_Pupil,[conc_Stim_reg conc_RampflipScaled_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
end
% BIC scores
conc_GLMfit.BIC(1) = n+(n.*log(2*pi))+(n.*log(sum(conc_GLMfit.stimResp.r.^2)./n))+(log(n).*(3+1)); % model 0. Just to see whether solely stim and resp model would explain sufficient variance
conc_GLMfit.BIC(2) = n+(n.*log(2*pi))+(n.*log(sum(conc_GLMfit.Boxc.r.^2)./n))+(log(n).*(4+1)); % model 1
conc_GLMfit.BIC(3) = n+(n.*log(2*pi))+(n.*log(sum(conc_GLMfit.Ramp.r.^2)./n))+(log(n).*(4+1)); % model 2
conc_GLMfit.BIC(4) = n+(n.*log(2*pi))+(n.*log(sum(conc_GLMfit.Ramp2thresh.r.^2)./n))+(log(n).*(4+1)); % model 3
conc_GLMfit.BIC(5) = n+(n.*log(2*pi))+(n.*log(sum(conc_GLMfit.Rampflip.r.^2)./n))+(log(n).*(4+1)); % model 4
conc_GLMfit.BIC(6) = n+(n.*log(2*pi))+(n.*log(sum(conc_GLMfit.Ramp2threshflip.r.^2)./n))+(log(n).*(4+1)); % model 5
conc_GLMfit.BIC(7) = n+(n.*log(2*pi))+(n.*log(sum(conc_GLMfit.BoxcScaledDown.r.^2)./n))+(log(n).*(4+1)); % model 6
conc_GLMfit.BIC(8) = n+(n.*log(2*pi))+(n.*log(sum(conc_GLMfit.RampScaled.r.^2)./n))+(log(n).*(4+1)); % model 7
conc_GLMfit.BIC(9) = n+(n.*log(2*pi))+(n.*log(sum(conc_GLMfit.RampflipScaled.r.^2)./n))+(log(n).*(4+1)); % model 8

% Variance inflation factor
if cfg.orthogonalise_predictors
    conc_GLMfit.VIF{1} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_Resp_reg]));
    conc_GLMfit.VIF{2} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_Boxc_reg conc_Resp_reg]));
    conc_GLMfit.VIF{3} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_Ramp_reg conc_Resp_reg]));
    conc_GLMfit.VIF{4} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_Ramp2thresh_reg conc_Resp_reg]));
    conc_GLMfit.VIF{5} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_Rampflip_reg conc_Resp_reg]));
    conc_GLMfit.VIF{6} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_Ramp2threshflip_reg conc_Resp_reg]));
    conc_GLMfit.VIF{7} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_BoxcScaledDown_reg conc_Resp_reg]));
    conc_GLMfit.VIF{8} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_RampScaled_reg conc_Resp_reg]));
    conc_GLMfit.VIF{9} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_RampflipScaled_reg conc_Resp_reg]));
else
    conc_GLMfit.VIF{1} = varianceInflationCoefficients([conc_Stim_reg conc_Resp_reg]);
    conc_GLMfit.VIF{2} = varianceInflationCoefficients([conc_Stim_reg conc_Boxc_reg conc_Resp_reg]);
    conc_GLMfit.VIF{3} = varianceInflationCoefficients([conc_Stim_reg conc_Ramp_reg conc_Resp_reg]);
    conc_GLMfit.VIF{4} = varianceInflationCoefficients([conc_Stim_reg conc_Ramp2thresh_reg conc_Resp_reg]);
    conc_GLMfit.VIF{5} = varianceInflationCoefficients([conc_Stim_reg conc_Rampflip_reg conc_Resp_reg]);
    conc_GLMfit.VIF{6} = varianceInflationCoefficients([conc_Stim_reg conc_Ramp2threshflip_reg conc_Resp_reg]);
    conc_GLMfit.VIF{7} = varianceInflationCoefficients([conc_Stim_reg conc_BoxcScaledDown_reg conc_Resp_reg]);
    conc_GLMfit.VIF{8} = varianceInflationCoefficients([conc_Stim_reg conc_RampScaled_reg conc_Resp_reg]);
    conc_GLMfit.VIF{9} = varianceInflationCoefficients([conc_Stim_reg conc_RampflipScaled_reg conc_Resp_reg]);
end

%%% -----------------------------------------------------------------------
% perform the same analysis but then binned by RT, bins are considered
% separate conditions
for ibin = 1:length(cfg.nbin)
    
    binning2use = nan(length(RT),1);
    switch cfg.bin2use
        case 'RT'
            [~, binning2use(validtr2use), ~] = binVal(RT(validtr2use), cfg.nbin(ibin), 'equal');
        case 'baselinePupil'
            [~, binning2use(validtr2use), ~] = binVal(baseline(validtr2use), cfg.nbin(ibin), 'equal');
    end
    binning2use(isnan(binning2use)) = 0;
    
    bin_GLMfit(ibin).nbin = cfg.nbin(ibin);
    bin_GLMfit(ibin).bin2use = cfg.bin2use;
    
    for ibinsize = 1:cfg.nbin(ibin)
        conc_Pupil=[];
        conc_Stim_reg=[]; conc_Boxc_reg=[]; conc_Ramp_reg=[]; conc_Resp_reg=[];
        
        trIdx =  find(binning2use == ibinsize);
        
        % concatenate regressors
        for itrial = 1:length(trIdx)
            conc_Pupil(end+1:end+length(st_pupil{trIdx(itrial)}),1) = st_pupil{trIdx(itrial)};
            conc_Stim_reg(end+1:end+length(st_pupil{trIdx(itrial)}),1) = Stim_reg{trIdx(itrial)};
            conc_Resp_reg(end+1:end+length(st_pupil{trIdx(itrial)}),1) = Resp_reg{trIdx(itrial)};
            conc_Boxc_reg(end+1:end+length(st_pupil{trIdx(itrial)}),1) = Boxc_reg{trIdx(itrial)};
            conc_Ramp_reg(end+1:end+length(st_pupil{trIdx(itrial)}),1) = Ramp_reg{trIdx(itrial)};
        end
        
        if cfg.meanCenter_predictors
            % mean center
            conc_Pupil = conc_Pupil - mean(conc_Pupil);
            conc_Stim_reg = conc_Stim_reg - mean(conc_Stim_reg);
            conc_Resp_reg = conc_Resp_reg - mean(conc_Resp_reg);
            conc_Boxc_reg = conc_Boxc_reg - mean(conc_Boxc_reg);
            conc_Ramp_reg = conc_Ramp_reg - mean(conc_Ramp_reg);
        end
        
        if cfg.normalise_pupil==1
            % normalise
            conc_Pupil = (conc_Pupil - mean(conc_Pupil))/std(conc_Pupil);
        elseif cfg.normalise_pupil == 2
            conc_Pupil = (conc_Pupil - mean(conc_Pupil));
        end
        
        
        %%% across trial regression
        n = length(conc_Pupil); % number of observations
        
        % parameters
        bin_GLMfit(ibin).GLM(ibinsize).binSize        =  length(find(trIdx));
        
        % three-param models
        if cfg.orthogonalise_predictors
            bin_GLMfit(ibin).GLM(ibinsize).stimResp        = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
            bin_GLMfit(ibin).GLM(ibinsize).Boxc            = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_Boxc_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
            bin_GLMfit(ibin).GLM(ibinsize).Ramp            = regstats(conc_Pupil,spm_orth([conc_Stim_reg conc_Ramp_reg conc_Resp_reg]),'linear',{'yhat','r','beta','tstat','rsquare'});
        else
            bin_GLMfit(ibin).GLM(ibinsize).stimResp        = regstats(conc_Pupil,[conc_Stim_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
            bin_GLMfit(ibin).GLM(ibinsize).Boxc            = regstats(conc_Pupil,[conc_Stim_reg conc_Boxc_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
            bin_GLMfit(ibin).GLM(ibinsize).Ramp            = regstats(conc_Pupil,[conc_Stim_reg conc_Ramp_reg conc_Resp_reg],'linear',{'yhat','r','beta','tstat','rsquare'});
        end
        
        % BIC scores
        bin_GLMfit(ibin).GLM(ibinsize).BIC(1) = n+(n.*log(2*pi))+(n.*log(sum(bin_GLMfit(ibin).GLM(ibinsize).stimResp.r.^2)./n))+(log(n).*(3+1)); % model 0. Just to see whether solely stim and resp model would explain sufficient variance
        bin_GLMfit(ibin).GLM(ibinsize).BIC(2) = n+(n.*log(2*pi))+(n.*log(sum(bin_GLMfit(ibin).GLM(ibinsize).Boxc.r.^2)./n))+(log(n).*(4+1)); % model 1
        bin_GLMfit(ibin).GLM(ibinsize).BIC(3) = n+(n.*log(2*pi))+(n.*log(sum(bin_GLMfit(ibin).GLM(ibinsize).Ramp.r.^2)./n))+(log(n).*(4+1)); % model 2
        
        if cfg.orthogonalise_predictors
            bin_GLMfit(ibin).GLM(ibinsize).VIF{1} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_Resp_reg]));
            bin_GLMfit(ibin).GLM(ibinsize).VIF{2} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_Boxc_reg conc_Resp_reg]));
            bin_GLMfit(ibin).GLM(ibinsize).VIF{3} = varianceInflationCoefficients(spm_orth([conc_Stim_reg conc_Ramp_reg conc_Resp_reg]));
        else
            bin_GLMfit(ibin).GLM(ibinsize).VIF{1} = varianceInflationCoefficients([conc_Stim_reg conc_Resp_reg]);
            bin_GLMfit(ibin).GLM(ibinsize).VIF{2} = varianceInflationCoefficients([conc_Stim_reg conc_Boxc_reg conc_Resp_reg]);
            bin_GLMfit(ibin).GLM(ibinsize).VIF{3} = varianceInflationCoefficients([conc_Stim_reg conc_Ramp_reg conc_Resp_reg]);
        end
    end
end
