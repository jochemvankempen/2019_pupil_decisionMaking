function batch_preprocess_CD(allsubj, subIdx, paths, files, blocks, badchans, dataset, erpfileExt)
% These scripts reproduce the analysis in the paper: van Kempen et al.,
% (2018) 'Behavioural and neural signatures of perceptual evidence
% accumulation are modulated by pupil-linked arousal'. 
% 
% Many of these scripts are based on the original scripts for the papers
% Newman et al. (2017), Journal of Neuroscience.
% https://github.com/gerontium/big_dots, and Loughnane et al. (2018), The
% Journal of Neuroscience
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
% preprocessing of CD paradigm.

% Triggers
targcodes = [8]; % contrast change onset event
resp_string = {'12','13'};

setAnalysisSettings_bigDots % use same settings as bigDots. 
% Note that this is only valid because I'm analysing behavioural data. If
% you analyse EEG data, these settings need to be altered

saveBigFile = 1; %if 0, doesn't save 8Hz erp. Smaller file sizes

savefilename = [paths.s(subIdx).base allsubj{subIdx} '_' dataset '_' erpfileExt '.mat'];

if exist(savefilename,'file')
    fprintf('file %s already exists\n', savefilename)
%     return
end

%% CSD

try
    load([paths.data 'csd_montage.mat'])
catch
    E = textread([paths.base 'chans_JBhead_96_noground.asc'],'%s');
    M = ExtractMontage([paths.base 'JBhead96_sym.csd'],E);  % reading in the montage for the CSD toolbox
    MapMontage(M);
    [G_CSD,H_CSD] = GetGH(M);
    
    save([paths.data 'csd_montage.mat'],'G_CSD','H_CSD');
end

%% Artifact rejection

ARchans = [1:97];
artifth = 100;

% frontal_chans = [1,2,65,66,33,34,35,36,67,3,37,4,38,5,39,6,40,7,68,41,42,45,46];
% ch_CPP = [53];
% ch_lr = [26,59,87;24,56,84];
% ch_rl = [24,56,84;26,59,87];

chanlocs    = readlocs('JBhead96_sym.loc');
nchan       = 97;
chanlocs    = chanlocs(1:nchan)';

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
% disp(subject_folder{isub})
   
for iblock=blocks %6%
    artifacts_this_block = []; artifactsET_this_block =[]; RT_this_block = []; trialThisBlock = [];
    
    disp(files(subIdx,iblock).subID)
    
    EEG = pop_loadbv(paths.s(subIdx).raw,files(subIdx,iblock).vhdr);
    loadbvSK_DN_JBhead_97
    EEG = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers    
    EEG.data = double(EEG.data);
    
    load([paths.s(subIdx).raw files(subIdx, iblock).mat]);
    
%     if ~isempty(files(subIdx, iblock).ET_files)
        FixEyelinkMessages_CD %then calculate and save it now
%     end
    
    % Eyetracker stuff
    % scres = 800 x 600: 400,300 is middle. 75 pixels is 3 degrees.
    screen_res = [800 600];
    ranger = 75;
    middle = screen_res/2;

    numev = length(EEG.event);
    % Fish out the event triggers and times
    clear trigs stimes
    for i=1:numev
        trigs(i)=EEG.event(i).type;
        stimes(i)=round(EEG.event(i).latency);
    end
    stimes(find(ismember(trigs,[-88,1,5])))=[];
    trigs(find(ismember(trigs,[-88,1,5])))=[];       
    
    % ET stuff
    if ~isempty(files(subIdx, iblock).ET_files)
        load(files(subIdx, iblock).ET_matfiles) %DN: load the ET mat file
        
        % fix missing trig issue
        EEG_trigs=[]; ET_trigs=[];
        for i = 1:length(EEG.event), EEG_trigs(i) = EEG.event(i).type; end
        for i = 1:length(event), ET_trigs(i) = event(i,2); end
        ET_trigs = ET_trigs(find(ET_trigs==8)); EEG_trigs = EEG_trigs(find(EEG_trigs==8));
        
        if length(ET_trigs)>length(EEG_trigs),
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
%         [output_cell,~,~] = command_window_text();
%         text = output_cell{length(output_cell)-1};
%         numsamp = sscanf(text,'%*s%*s%*s%*s%*s%*f%*s%*s%f');
%         if numsamp<24
%             beep
%             disp([allsubj{subIdx},', block ',num2str(iblock),': ET sync issue'])
%             figure, plot(ET_trigs), hold on, plot(EEG_trigs)
%             keyboard
%         end
        
        ET_data     = EEG_temp.data([nchan+1:nchan+4],:);
        tmpETevent  = EEG_temp.event;
%         clear EEG_temp
        pause(1)
    end
    
    % Baseline position?
    temp = ET_data([2,3],stimes(1):stimes(end));
    ET_data_mean(1) = mean(temp(1,find(temp(1,:)>1)),2); ET_data_mean(2) = mean(temp(2,find(temp(2,:)>1)),2);
    ET_data(2, find(ET_data(2,:)>0)) = ET_data(2, find(ET_data(2,:)>0))-ET_data_mean(1)+400;
    ET_data(3, find(ET_data(3,:)>0)) = ET_data(3, find(ET_data(3,:)>0))-ET_data_mean(2)+300;

  
    % interpolate bad channels
    if ~isempty(badchans)
        EEG.chanlocs = chanlocs;
        EEG=eeg_interp(EEG,[badchans],'spherical');
    end
    
    %filtering
    % There are peaks at 25, 50, 75, 100, 180, 200. Oddly a suppression at 60.
    % Get rid of mains frequency
    EEG = pop_eegfiltnew(EEG,59,61,[],1,0,0,0);
    EEG = pop_eegfiltnew(EEG,119,121,[],1,0,0,0);
    EEG = pop_eegfiltnew(EEG,179,181,[],1,0,0,0);
    % HP Filter
    if HPF
        EEG = pop_eegfiltnew(EEG,HPFcutoff,0,[]); % filter to 0.1Hz
    end
      
    % Get rid of mains frequency
    EEG_LPF_8Hz = EEG; EEG_LPF_35Hz = EEG; clear EEG;
    
    % First LP Filter
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
    
    %%% the pupil data is only recorded from the startpoint of the first
    %%% trial, until the last trial. This causes problems with filtering of
    %%% eye data etc. Therefore, we delete those.. 
    targtrigs([1 end]) = [];
        
    if ~isempty(files(subIdx, iblock).ET_files)
        [ET_data] = interpolateET(ET_data, [stimes(targtrigs)' trigs(targtrigs)'] , tmpETevent, blinkevent, fs, files(subIdx,iblock).ET_matfiles, dataset);
    else
        ET_data = zeros(8,size(EEG_LPF_35Hz.data,2));
    end
    
    for n=1:length(targtrigs)
        clear ep_LPF_8Hz ep_LPF_35Hz ep_LPF_8Hz_CSD ep_LPF_35Hz_CSD
        locktime = stimes(targtrigs(n));
        if (targtrigs(n)+1 <= length(trigs)) && trigs(targtrigs(n)+1)==12
            response_time = stimes(targtrigs(n)+1)-locktime; % time in samples from beginning of motion to response.
            response_time = floor(response_time);
            if response_time>rtlim(2)*fs
                response_time = rtlim(2)*fs;
            end
        else
            response_time = rtlim(2)*fs;
        end
        RT_this_block(n) = response_time;
        try
            ep_LPF_8Hz  = EEG_LPF_8Hz.data(:,locktime+ts);   % chop out an epoch
            ep_LPF_35Hz = EEG_LPF_35Hz.data(:,locktime+ts);
            ep_ET       = ET_data(:,locktime+ts);
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
            ET_trials(:,:,numtr)        = zeros(8,ERP_samps);

            continue;
        end
        
        % new trial
        numtr = numtr+1;
        trialThisBlock = [trialThisBlock numtr];
        
%         BLamp           = mean(ep_LPF_8Hz(:,find(t>=BL_time(1) & t<BL_time(2))),2); % record baseline amplitude for each channel, ONLY FOR ART REJECT
%         ep_LPF_8Hz_BL   = ep_LPF_8Hz - repmat(BLamp,[1,length(t)]);
        
        BLamp           = mean(ep_LPF_35Hz(:,find(t>=BL_time(1) & t<BL_time(2))),2); % record baseline amplitude for each channel, ONLY FOR ART REJECT
        ep_LPF_35Hz_BL  = ep_LPF_35Hz - repmat(BLamp,[1,length(t)]);
        
        ep_test = [find(ts<=(response_time+0.1*fs))];
        if isempty(ep_test)
            disp('Empty epoch for art rejection')
            keyboard
        end
        
        for ifield = 1:length(artifact_times)
            
            switch artifact_times{ifield}
                case {'neg500_0','neg100_0'}
                    %define time window for artifact rejection
                    eval(['ts_artif = (ts >= tt.' artifact_times{ifield} '(1) & ts <= tt.' artifact_times{ifield} '(2));'])
                otherwise
                    eval(['ts_artif = (ts >= tt.' artifact_times{ifield} '(1) & ts <= (response_time + tt.' artifact_times{ifield} '(2)));'])
            end
                    
            % find artifacts in specified time windows in EEG channels
            eval(['artifchans_thistrial.' artifact_times{ifield} ' = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts_artif)))  ,[],2)>artifth));'])
            
            % concatenate them across blocks
            eval(['artifchans.' artifact_times{ifield} ' = [artifchans.' artifact_times{ifield} ' artifchans_thistrial.' artifact_times{ifield} '];'])
            
            % 0 = reject, 1 = keep
            eval(['if length(artifchans_thistrial.' artifact_times{ifield} ') > 0, artrej.' artifact_times{ifield} '(numtr) = 0;  else artrej.' artifact_times{ifield} '(numtr) = 1; end'])

            % find trials where fixation was broken
            artif_ET = find(ep_ET(2,find( ts_artif )) < middle(1)-ranger  |  ep_ET(2,find( ts_artif )) > middle(1)+ranger  | ... % x
                ep_ET(3,find( ts_artif )) < middle(2)-ranger  |  ep_ET(3,find( ts_artif )) > middle(2)+ranger ); % y
            
            % 0 = reject, 1 = keep
            eval(['if length(artif_ET) > 0, artrej_ET.' artifact_times{ifield} '(numtr) = 0; else artrej_ET.' artifact_times{ifield} '(numtr) = 1;      end'])
        end
        
        % for plotting in interpolate ET data script
        if length(artifchans_thistrial.neg500_RT_200)   > 0,    artifacts_this_block(n) = 0;       else artifacts_this_block(n) = 1;   end
        if artrej_ET.neg500_RT_200(numtr)               == 0,   artifactsET_this_block(n) = 0;     else artifactsET_this_block(n) = 1; end
        
        ep_LPF_8Hz  = double(ep_LPF_8Hz);
        ep_LPF_35Hz = double(ep_LPF_35Hz);

        ep_LPF_8Hz_CSD  = CSD(ep_LPF_8Hz,G_CSD,H_CSD);
        ep_LPF_35Hz_CSD = CSD(ep_LPF_35Hz,G_CSD,H_CSD);

        erp_LPF_8Hz(:,:,numtr)      = ep_LPF_8Hz;
        erp_LPF_35Hz(:,:,numtr)     = ep_LPF_35Hz;
        erp_LPF_8Hz_CSD(:,:,numtr)  = ep_LPF_8Hz_CSD;
        erp_LPF_35Hz_CSD(:,:,numtr) = ep_LPF_35Hz_CSD;
        ET_trials(:,:,numtr)        = ep_ET(:,:);
        
        allTrig(numtr) = trigs(targtrigs(n));
        allblock_count(numtr) = iblock;
        allITI(numtr) = trialITI(n);
        allUD(numtr) = par.UPorDOWN;
        
        try % change this
            if trigs(targtrigs(n)+1)==12
                allrespLR(numtr) = 1;
                allRT(numtr) = stimes(targtrigs(n)+1)-stimes(targtrigs(n));
            elseif trigs(targtrigs(n)+1)==13 % they pressed the wrong button
                allrespLR(numtr) = 2;
                allRT(numtr) = stimes(targtrigs(n)+1)-stimes(targtrigs(n));
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
        [~,tmpET_artrej] = interpolateET(ET_data, [stimes(targtrigs)' trigs(targtrigs)'] , tmpETevent, blinkevent, fs, files(subIdx,iblock).ET_matfiles, dataset, artifacts_this_block, artifactsET_this_block, RT_this_block, tt);
            
    
        for ifield = 1:length(artifact_times)
            eval(['artrej_pupil.' artifact_times{ifield} ' = [artrej_pupil.' artifact_times{ifield} '  trialThisBlock( tmpET_artrej.' artifact_times{ifield} ' )];'])
            
            switch dataset
                case 'CD'
                    % for some reason, pupil diameter only started to be
                    % recorded during the first trial and was stopped
                    % during the last trial of the block. Mark those trials
                    % as errors
                    eval(['artrej_pupil.' artifact_times{ifield} ' = [artrej_pupil.' artifact_times{ifield} '  trialThisBlock( [1 end] )];'])
            end
        end
        
    else
        
        for ifield = 1:length(artifact_times)
            eval(['artrej_pupil.' artifact_times{ifield} ' = [artrej_pupil.' artifact_times{ifield} '  trialThisBlock];'])
        end

    end
    blockIdx = [blockIdx ; ones(length(trialThisBlock),1)*iblock];
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
figure(1),clf;
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

saveFigName = ['fig_artRej_' allsubj{subIdx} erpfileExt];
fileExt = 'jpeg';
print(gcf,['-d' fileExt],[paths.s(subIdx).fig saveFigName '.' fileExt])

erp_LPF_8Hz         = single(erp_LPF_8Hz);
erp_LPF_35Hz        = single(erp_LPF_35Hz);
erp_LPF_8Hz_CSD     = single(erp_LPF_8Hz_CSD);
erp_LPF_35Hz_CSD    = single(erp_LPF_35Hz_CSD);

if saveBigFile
    save(savefilename,...
        'erp_LPF_8Hz','erp_LPF_35Hz','erp_LPF_8Hz_CSD','erp_LPF_35Hz_CSD', ...
        'allRT','allrespLR','allTrig','allblock_count','t','ET_trials', 'blockIdx', 'allITI', 'allUD',  ...
        'artifchans',...
        'artrej',...
        'artrej_ET', ...
        'artrej_pupil')
else
    save(savefilename,...
        'erp_LPF_35Hz','erp_LPF_35Hz_CSD', ...
        'allRT','allrespLR','allTrig','allblock_count','t','ET_trials', 'blockIdx', 'allITI', 'allUD',  ...
        'artifchans',...
        'artrej',...
        'artrej_ET', ...
        'artrej_pupil')
end
return;