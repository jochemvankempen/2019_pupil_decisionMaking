function batch_preprocess(allsubj, subIdx, paths, files, blocks, badchans, location, dataset, fileExt)
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
% preprocessing of CDT paradigm.

%
% Trigger 1: coherence 50, motion dir 270, ITI 3.06, patch 1
% Trigger 2: coherence 50, motion dir 270, ITI 3.06, patch 2
% Trigger 3: coherence 50, motion dir 270, ITI 5.17, patch 1
% Trigger 4: coherence 50, motion dir 270, ITI 5.17, patch 2
% Trigger 5: coherence 50, motion dir 270, ITI 7.29, patch 1
% Trigger 6: coherence 50, motion dir 270, ITI 7.29, patch 2
tic

plotVisible = 'off';
plotFigures = 1;

% setPlotSettings
setAnalysisSettings_bigDots

saveBigFile = 1; %if 0, doesn't save 8Hz erp. Smaller file sizes

savefilename = [paths.s(subIdx).savebase allsubj{subIdx} fileExt '.mat'];

if exist(savefilename,'file')
    fprintf('file %s already exists\n', savefilename)
%     return
end

%% CSD prep

try
    load([paths.readdata 'csd_montage.mat'])
catch
    E = textread('chans64.asc','%s');
    M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',E);  % reading in the montage for the CSD toolbox
    % MapMontage(M);
    [G_TCD,H_TCD] = GetGH(M);
    
    E = textread('chans65_monash.asc','%s');
    M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',E);  % reading in the montage for the CSD toolbox
    % MapMontage(M);
    [G_monash,H_monash] = GetGH(M);
    
    save([paths.savedata 'csd_montage.mat'],'G_TCD','H_TCD','G_monash','H_monash');
end

switch location
    case 'Monash'
        G_CSD = G_monash;
        H_CSD = H_monash;
    case 'TCD'
        G_CSD = G_TCD;
        H_CSD = H_TCD;
end
%% Artifact rejection

ARchans = [1:64];  % just the ones we care about, and RH opposite LH to be symmetric
artifth = 100;

switch location
    case 'Monash'
        chanlocs = readlocs('actiCAP65_ThetaPhi.elp','filetype','besa'); %DN for actiCAP
        nchan   = 65;
    case 'TCD'
        chanlocs = readlocs('cap64.loc');
        nchan   = 64;
end
chanlocs = chanlocs(1:nchan)';

%% Initialise
% # % Make sure all initialised
numtr=0;
allRT=[]; allrespLR=[]; allTrig=[]; allblock_count = [];
erp_LPF_8Hz = []; erp_LPF_35Hz = []; erp_LPF_8Hz_CSD = []; erp_LPF_35Hz_CSD = [];

ET_trials = [];
blockIdx=[];

% keep track of channels on which the threshold is exceeded, causing trial rejection
% initialise artifact rejection vectors
for ifield = 1:length(artifact_times)
    eval(['artrej.'         artifact_times{ifield} ' = [];'])
    eval(['artifchans.'     artifact_times{ifield} ' = [];'])
    eval(['artrej_ET.'      artifact_times{ifield} ' = [];'])
    eval(['artrej_pupil.'   artifact_times{ifield} ' = [];'])
end

%% Begin loop
for iblock= blocks 
    artifacts_this_block = []; artifactsET_this_block =[]; RT_this_block = []; trial_this_block = [];
    
    disp(files(subIdx,iblock).subID)
    
    % load EEG
    switch location
        case 'Monash'
            EEG = pop_loadbv(paths.s(subIdx).raw,files(subIdx,iblock).vhdr);
            loadbvSK_DN % this takes 65 channels
            EEG = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers
            while 1
                if EEG.event(1).type==1 | EEG.event(1).type==-88
                    EEG.event(1) = [];
                else
                    break;
                end
            end
            
        case 'TCD'
            filename=[paths.s(subIdx).raw files(subIdx,iblock).bdf];
            EEG = pop_biosig(filename, 'blockepoch', 'off','channels',[1:nchan]);
    end
    EEG.data = double(EEG.data);
    
    load([paths.s(subIdx).raw files(subIdx, iblock).mat],'trialCond','par');
    
    if 1%~isempty(files(subIdx, iblock).ET_files)
        [~, name] = fileparts(files(subIdx, iblock).ET_files);
        if strcmpi(name, 'SW50M1') || strcmpi(name, 'SW50M16')
            screen_res = [1024 768];%somehow there is no information in the first block. Screen res taken from second block
        elseif strcmpi(name, 'NB52M1') %|| strcmpi(name, 'SW50M16')
            screen_res = [1024 768];%somehow there is no information in the first block. Screen res taken from second block
        elseif strcmpi(allsubj{subIdx}, 'HN015') || strcmpi(allsubj{subIdx}, 'HN019') || strcmpi(allsubj{subIdx}, 'HN023')
            screen_res = [1024 768];%No information in any of this subjects' blocks. assumer the same as the other subjects
        end
        fid = fopen(files(subIdx, iblock).ET_files);
        ET_text = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s','Headerlines',22,'ReturnOnError',0);
        fclose(fid);
        for i = 1:size(ET_text{1,3},1)
            if strcmp('GAZE_COORDS',ET_text{1,3}(i))
                screen_res(1) = str2num(cell2mat(ET_text{1,6}(i)))+1;
                screen_res(2) = str2num(cell2mat(ET_text{1,7}(i)))+1;
                continue
            end
        end
        if screen_res(1)==1024, ranger = 76; elseif screen_res(1)==1280, ranger = 98; else disp(screen_res), keyboard, end
        middle = screen_res/2;
        
        if 1%~exist(files(subIdx, iblock).ET_matfiles, 'file') %DN: if ET matfile NOT been saved
            fixEyelinkMessages %then calculate and save it now
        end
    else
        screen_res = [1024 768];%
        ranger = 76;
        middle = screen_res/2;
    end
    
    trialCond = trialCond+100;
    
    switch location
        case 'Monash'
            numev = length(EEG.event);
            % Fish out the event triggers and times
            clear trigs stimes
            for i=1:numev
                trigs(i)=EEG.event(i).type;
                stimes(i)=round(EEG.event(i).latency);
            end
            
        case 'TCD'
            EEG2 = pop_resample(EEG,fs);
            resamp_data = [];
            for elec = 1:size(EEG.data,1)
                resamp_data(elec,:) = resample(EEG.data(elec,:),500,512,0);
            end
            EEG.data = resamp_data;
            EEG.pnts = EEG2.pnts;
            EEG.times = EEG2.times;
            EEG.srate = EEG2.srate;
            EEG.event = EEG2.event;
            
            if size(resamp_data,2)~=EEG.pnts
                disp('Resample size mismatch')
                keyboard
            end
            
            numev = length(EEG2.event);
            % Fish out the event triggers and times
            clear trigs stimes
            for i=1:numev
                trigs(i)=EEG2.event(i).type;
                stimes(i)=round(EEG2.event(i).latency);
            end
            clear EEG2 % unnecessary now
    end
    
    % interpolate bad channels
    if ~isempty(badchans)
        EEG.chanlocs = chanlocs;
        EEG=eeg_interp(EEG,[badchans],'spherical');
    end
    
    % ET stuff
    if ~isempty(files(subIdx, iblock).ET_files)
        load(files(subIdx, iblock).ET_matfiles) %DN: load the ET mat file
        
        % fix missing trig issue
        EEG_trigs=[]; ET_trigs=[];
        for i = 1:length(EEG.event), EEG_trigs(i) = EEG.event(i).type; end
        for i = 1:length(event), ET_trigs(i) = event(i,2); end
        ET_trigs = ET_trigs(find(ET_trigs>100 & ET_trigs<107)); EEG_trigs = EEG_trigs(find(EEG_trigs>100 & EEG_trigs<107));
        
        if length(ET_trigs)>length(EEG_trigs)
            % compare et_event with targtriggers
            clear x trigIdx
            for itrig = 1:length(EEG_trigs)
                if find(ET_trigs(itrig:end)==EEG_trigs(itrig),1,'first')
                    x(itrig) = find(ET_trigs(itrig:end)==EEG_trigs(itrig),1,'first');
                end
                trigIdx(itrig) = itrig + x(itrig)-1;
            end
            last_event = ET_trigs(trigIdx(end));
            if ET_trigs(end)==ET_trigs(end-1),
                event = event(1:end-2,:);
                save(files(subIdx, iblock).ET_matfiles,'event','-append');
            end
        end
        
        
        plotter = 0;
        %Add an extra 4 rows into the EEG struct - 'TIME'
        %'GAZE_X' 'GAZE_Y' 'AREA'. This will add these as extra channels onto %EEG.data. So the final channel is the pupil area (i.e. diameter):
        EEG_temp = pop_importeyetracker(EEG,files(subIdx, iblock).ET_matfiles,[first_event ...
            last_event],[1:4] ,{'TIME' 'GAZE_X' 'GAZE_Y' 'AREA'},1,1,0,plotter);
        if 0 % cannot be used with parfor, only necessary to check import of eye events
            [output_cell,~,~] = command_window_text();
            text = output_cell{length(output_cell)-1};
            numsamp = sscanf(text,'%*s%*s%*s%*s%*s%*f%*s%*s%f');
            if numsamp<30
                beep
                disp([allsubj{s},', block ',num2str(iblock),': ET sync issue'])
                figure, plot(ET_trigs), hold on, plot(EEG_trigs)
                keyboard
            end
        end
        
        ET_data     = EEG_temp.data([nchan+1:nchan+4],:);
        tmpETevent  = EEG_temp.event;
        clear EEG_temp
        pause(1)
    end
    % Check fft for mains frequencies etc
    %     ep = squeeze(EEG.data(31,:)); % time
    %     nfft = size(ep,2);
    %     fftx = abs(fft(ep,[],2))./(nfft/2);
    %     fftx = fftx(:,1:ceil((nfft+1)/2));
    %     freq_temp = (0:ceil((nfft+1)/2)-1)*fs/nfft;
    %     figure, plot(freq_temp,fftx)
    
    % Get rid of mains frequency
    EEG = pop_eegfiltnew(EEG,49,51,[],1,0,0,0);
    
    % First HP Filter
    if HPF
        EEG = pop_eegfiltnew(EEG,HPFcutoff,0,[]); % filter to 0.1Hz
        disp('HPF finished')
    end
    
    EEG_LPF_8Hz = EEG; EEG_LPF_35Hz = EEG; clear EEG;
    if LPF
        EEG_LPF_8Hz = pop_eegfiltnew(EEG_LPF_8Hz,0,LPFcutoff_8Hz,[]);
        EEG_LPF_35Hz = pop_eegfiltnew(EEG_LPF_35Hz,0,LPFcutoff_35Hz,[]);
    end
    
    % average-reference the whole continuous data (safe to do this now after interpolation and filtering):
    EEG_LPF_8Hz.data    = EEG_LPF_8Hz.data -    repmat(mean(EEG_LPF_8Hz.data(:,:),1),[nchan,1]);
    EEG_LPF_35Hz.data   = EEG_LPF_35Hz.data -   repmat(mean(EEG_LPF_35Hz.data(:,:),1),[nchan,1]);
    
    targtrigs = [];
    for n=1:length(trigs)
        if any(targcodes(:)==trigs(n))
            targtrigs = [targtrigs n];
        end
    end
    
    if trigs(targtrigs(end))==trialCond(1)
        motion_on = targtrigs(1:end-1); % GL: indices of trigs when motion on. get rid of last trig, it was a repeat
    else
        motion_on = targtrigs;
    end
    
    % this interpolation script is currently run twice, not very
    % efficient.. But doesn't take too long
    if ~isempty(files(subIdx, iblock).ET_files)
        [ET_data] = interpolateET(ET_data, [stimes(motion_on)' trigs(motion_on)'] , tmpETevent, eyelinkevent, fs, files(subIdx,iblock).ET_matfiles, dataset);
    else
        ET_data = zeros(10,size(EEG_LPF_35Hz.data,2));
    end
    
    for n=1:length(motion_on)
        clear ep_LPF_8Hz ep_LPF_35Hz ep_LPF_8Hz_CSD ep_LPF_35Hz_CSD
        locktime = stimes(motion_on(n));
        if motion_on(n)<length(trigs)
            if trigs(motion_on(n)+1)==12
                response_time = stimes(motion_on(n)+1)-locktime; % time in samples from beginning of motion to response.
                response_time = floor(response_time);
                if response_time>rtlim(2)*fs
                    response_time = rtlim(2)*fs;
                end
            else
                response_time = rtlim(2)*fs;
            end
        else
            response_time = rtlim(2)*fs;
        end
        RT_this_block(n) = response_time;
        try
            ep_LPF_8Hz  = EEG_LPF_8Hz.data(:,locktime+ts);   % chop out an epoch
            ep_LPF_35Hz = EEG_LPF_35Hz.data(:,locktime+ts);
            ep_ET       = ET_data(:,locktime+ts_pupil);
            %
        catch
            disp('EEG ended too soon2')
            numtr = numtr+1;
            allTrig(numtr) = 0;
            targMotion(numtr) = 0;
            allblock_count(numtr) = iblock;
            allrespLR(numtr) = 0;
            allRT(numtr) = 0;
            
            erp_LPF_8Hz(:,:,numtr)      = zeros(nchan,ERP_samps);
            erp_LPF_35Hz(:,:,numtr)     = zeros(nchan,ERP_samps);
            erp_LPF_8Hz_CSD(:,:,numtr)  = zeros(nchan,ERP_samps);
            erp_LPF_35Hz_CSD(:,:,numtr) = zeros(nchan,ERP_samps);
            ET_trials(:,:,numtr)        = zeros(10,ERP_samps);
            
            continue;
        end
        
        % new trial
        numtr = numtr+1;
        trial_this_block = [trial_this_block numtr];
        
        BLamp           = mean(ep_LPF_35Hz(:,find(t>=BL_time(1) & t<BL_time(2))),2); % record baseline amplitude for each channel, ONLY FOR ART REJECT
        ep_LPF_35Hz_BL  = ep_LPF_35Hz - repmat(BLamp,[1,length(t)]);
                
        ep_test = [find(ts<=(response_time+0.1*fs))];
        if isempty(ep_test)
            disp('Empty epoch for art rejection')
            keyboard
        end
        
        %loop over previously defined time windows (struct tt, defined in
        %setAnalysisSettings_bigDots)
        for ifield = 1:length(artifact_times)
            
            switch artifact_times{ifield}
                case {'neg500_0','neg100_0'}
                    %define time window for artifact rejection
                    eval(['ts_artif = (ts >= tt.' artifact_times{ifield} '(1) & ts <= tt.' artifact_times{ifield} '(2));'])
                    eval(['ts_pupil_artif = (ts_pupil >= tt.' artifact_times{ifield} '(1) & ts_pupil <= tt.' artifact_times{ifield} '(2));'])
                otherwise
                    % add 2nd element to RT
                    eval(['ts_artif = (ts >= tt.' artifact_times{ifield} '(1) & ts <= (response_time + tt.' artifact_times{ifield} '(2)));'])
                    eval(['ts_pupil_artif = (ts_pupil >= tt.' artifact_times{ifield} '(1) & ts_pupil <= (response_time + tt.' artifact_times{ifield} '(2)));'])
            end
            
            % find artifacts in specified time windows in EEG channels
            eval(['artifchans_thistrial.' artifact_times{ifield} ' = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts_artif)))  ,[],2)>artifth));'])
            
            % concatenate them across blocks
            eval(['artifchans.' artifact_times{ifield} ' = [artifchans.' artifact_times{ifield} ' artifchans_thistrial.' artifact_times{ifield} '];'])
            
            % 0 = reject, 1 = keep
            eval(['if length(artifchans_thistrial.' artifact_times{ifield} ') > 0, artrej.' artifact_times{ifield} '(numtr) = 0;  else artrej.' artifact_times{ifield} '(numtr) = 1; end'])
            
            % find trials where fixation was broken
            artif_ET = find(ep_ET(2,find( ts_pupil_artif )) < middle(1)-ranger  |  ep_ET(2,find( ts_pupil_artif )) > middle(1)+ranger  | ... % x
                ep_ET(3,find( ts_pupil_artif )) < middle(2)-ranger  |  ep_ET(3,find( ts_pupil_artif )) > middle(2)+ranger ); % y
            
            % 0 = reject, 1 = keep
            eval(['if length(artif_ET) > 0, artrej_ET.' artifact_times{ifield} '(numtr) = 0; else artrej_ET.' artifact_times{ifield} '(numtr) = 1;      end'])
        end
        
        % for plotting in interpolate ET data script
        if length(artifchans_thistrial.neg500_RT_200)    > 0, artifacts_this_block(n) = 0;       else artifacts_this_block(n) = 1;   end
        if artrej_ET.neg500_RT_200(numtr) == 0, artifactsET_this_block(n) = 0;     else artifactsET_this_block(n) = 1; end
        
        ep_LPF_8Hz  = double(ep_LPF_8Hz);
        ep_LPF_35Hz = double(ep_LPF_35Hz);
        
        ep_LPF_8Hz_CSD  = CSD(ep_LPF_8Hz,G_CSD,H_CSD);
        ep_LPF_35Hz_CSD = CSD(ep_LPF_35Hz,G_CSD,H_CSD);
        
        erp_LPF_8Hz(:,:,numtr)      = ep_LPF_8Hz;
        erp_LPF_35Hz(:,:,numtr)     = ep_LPF_35Hz;
        erp_LPF_8Hz_CSD(:,:,numtr)  = ep_LPF_8Hz_CSD;
        erp_LPF_35Hz_CSD(:,:,numtr) = ep_LPF_35Hz_CSD;
        ET_trials(:,:,numtr)        = ep_ET(:,:);
        
        allTrig(numtr) = trigs(motion_on(n));
        allblock_count(numtr) = iblock;
        
        try % change this
            if trigs(motion_on(n)+1)==12
                allrespLR(numtr) = 1;
                allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));
            elseif trigs(motion_on(n)+1)==13 % they pressed the wrong button
                allrespLR(numtr) = 2;
                allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));
            elseif trigs(motion_on(n)+1)==28 % fixation break (In CD paradigm)
                breakFix(numtr) = 1; 
            else
                allrespLR(numtr) = 3; % no response, to mark it out from artifact trials.
                allRT(numtr) = 0;
            end
        catch
            allrespLR(numtr) = 0;
            allRT(numtr) = 0;
        end
    end
    if ~isempty(files(subIdx, iblock).ET_files)
        [~,tmpET_artrej] = interpolateET(ET_data, [stimes(motion_on)' trigs(motion_on)'] , tmpETevent, eyelinkevent, fs, files(subIdx,iblock).ET_matfiles, dataset, artifacts_this_block, artifactsET_this_block, RT_this_block, tt);
        
        
        for ifield = 1:length(artifact_times)
            eval(['artrej_pupil.' artifact_times{ifield} ' = [artrej_pupil.' artifact_times{ifield} '  trial_this_block( tmpET_artrej.' artifact_times{ifield} ' )];'])
        end
    else
        for ifield = 1:length(artifact_times)
            eval(['artrej_pupil.' artifact_times{ifield} ' = [artrej_pupil.' artifact_times{ifield} '  trial_this_block];'])
        end
    end
    blockIdx = [blockIdx ; ones(length(trial_this_block),1)*iblock];
end

tmpartrej_pupil = artrej_pupil;
clear artrej_pupil

for ifield = 1:length(artifact_times)
    eval(['artrej_pupil.' artifact_times{ifield} ' = ones(numtr, 1);'])
    eval(['artrej_pupil.' artifact_times{ifield} '(tmpartrej_pupil.' artifact_times{ifield} ') = 0;'])
end


rejected_trials = length(find(artrej.neg100_RT_200==0));
% plot artifact trials per drug condition, interpolated channels are
% boldface
if plotFigures
    figure(1),clf;
    set(gcf,'visible',plotVisible)
    
    subplot(2,1,1)
    hist(artifchans.neg100_RT_200,[1:nchan]); title([allsubj{subIdx} ': ' num2str(rejected_trials) ' artifacts = ',num2str(round(100*(rejected_trials/length(artrej.neg100_RT_200)))),'%']) % s from runafew
    disp([allsubj{subIdx},' number of trials: ',num2str(length(artrej.neg100_RT_200))])
    
    [counts,centers] = hist(artifchans.neg100_RT_200,[1:nchan]);
    subplot(2,1,2)
    topoplot(counts,chanlocs,'plotchans',[1:nchan],'electrodes','numbers');
    if ~isempty(badchans)
        H = get(gca);
        clear tmpH
        tmpH = {H.Children(1:nchan).String};
        for iC = 1:length(badchans)
            cIdx = strcmpi(num2str(badchans(iC)), [tmpH(:)']);
            H.Children(cIdx).FontWeight = 'bold';
            H.Children(cIdx).FontSize = 14;
        end
    end
    colorbar;
    title([allsubj{subIdx}],'interpret','none')
    
    
    D.fig = [paths.s(subIdx).fig 'preprocess' filesep];
    if ~exist(D.fig,'dir'),mkdir(D.fig),end
    
    saveFigName = ['fig_artRej_' allsubj{subIdx} fileExt];
    fileExt = 'jpeg';
    print(gcf,['-d' fileExt],[paths.s(subIdx).fig saveFigName '.' fileExt])
end

erp_LPF_8Hz         = single(erp_LPF_8Hz);
erp_LPF_35Hz        = single(erp_LPF_35Hz);
erp_LPF_8Hz_CSD     = single(erp_LPF_8Hz_CSD);
erp_LPF_35Hz_CSD    = single(erp_LPF_35Hz_CSD);


%% reorganise chanlocs
switch dataset
    case {'bigDots'}
        switch location
            case 'Monash'
                chanlocs_TCD = readlocs('cap64.loc');
                chanlocs_Monash = readlocs('actiCAP65_ThetaPhi.elp','filetype','besa'); %DN for actiCAP
                
                counter = 1; counter2 = 1;
                order_TCD=[]; order_Monash=[]; unmatched_elecs_Monash=[];
                for elec_Monash = 1:65
                    yesser = 0;
                    for elec_TCD = 1:64
                        if strcmp(chanlocs_TCD(1,elec_TCD).labels,chanlocs_Monash(1,elec_Monash).labels)
                            order_TCD(counter) = elec_TCD;
                            order_Monash(counter) = elec_Monash;
                            counter = counter+1;
                            yesser = 1;
                        end
                    end
                    if yesser==0
                        unmatched_elecs_Monash(counter2) = elec_Monash;
                        counter2 = counter2+1;
                    end
                end
                % tester = zeros(65,1);
                % figure
                % topoplot(tester,chanlocs_Monash,'maplimits', ...
                %     [min(tester)  max(tester)],'electrodes','numbers');
                % figure
                % topoplot(tester,chanlocs_Monash,'maplimits', ...
                %     [min(tester)  max(tester)],'electrodes','labels');
                
                erp_TCD = single(zeros(size(erp_LPF_8Hz)));
                erp_TCD(order_TCD,:,:) = erp_LPF_8Hz(order_Monash,:,:);
                erp_LPF_8Hz = erp_TCD(1:64,:,:);
                erp_TCD = single(zeros(size(erp_LPF_35Hz)));
                erp_TCD(order_TCD,:,:) = erp_LPF_35Hz(order_Monash,:,:);
                erp_LPF_35Hz = erp_TCD(1:64,:,:);
                erp_TCD = single(zeros(size(erp_LPF_8Hz_CSD)));
                erp_TCD(order_TCD,:,:) = erp_LPF_8Hz_CSD(order_Monash,:,:);
                erp_LPF_8Hz_CSD = erp_TCD(1:64,:,:);
                erp_TCD = single(zeros(size(erp_LPF_35Hz_CSD)));
                erp_TCD(order_TCD,:,:) = erp_LPF_35Hz_CSD(order_Monash,:,:);
                erp_LPF_35Hz_CSD = erp_TCD(1:64,:,:);
                
                % % Baseline erp
                % baseline_erp = mean(erp_LPF_35Hz(:,find(t>=BL_time(1) & t<=BL_time(2)),:),2);
                % erp_temp = erp_LPF_35Hz-repmat(baseline_erp,[1,size(erp_LPF_35Hz,2),1]); % baseline full erp
                %
                % erp_temp = squeeze(mean(erp_temp(:,:,find(artrej_BL_resp==1)),3));
                % figure
                % plottopo(erp_temp(:,:),'chanlocs',chanlocs_TCD,'limits',[t(1) t(end) ...
                %     min(min(erp_temp(:,:)))  max(max(erp_temp(:,:)))], ...
                %     'title',['ERP'],'ydir',1)
        end
end

if saveBigFile
    save(savefilename,...
        'erp_LPF_8Hz','erp_LPF_35Hz','erp_LPF_8Hz_CSD','erp_LPF_35Hz_CSD', ...
        'allRT','allrespLR','allTrig','allblock_count','t','ET_trials', 'blockIdx', ...
        'artifchans',...
        'artrej',...
        'artrej_ET', ...
        'artrej_pupil')
else
    save(savefilename,...
        'erp_LPF_35Hz','erp_LPF_35Hz_CSD', ...
        'allRT','allrespLR','allTrig','allblock_count','t','ET_trials', 'blockIdx', ...
        'artifchans',...
        'artrej',...
        'artrej_ET', ...
        'artrej_pupil')
end
toc
return;