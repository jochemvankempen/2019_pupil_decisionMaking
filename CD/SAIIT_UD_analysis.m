% SAIIT_UD analysis

% par.CD_RESP  = 1;
% par.CD_FIXON = 4;
% par.CD_TGOFF = 5;   % target off
% par.CD_TG = 8;   % target   % one for each target type
% par.CD_BUTTONS = [12 13];
% par.CD_BEEP = 29; % shouldn't be in this task.
% par.secs_btw_targs = [4 7 10];
% par.targrampdur = 1.2;  % in sec. Return ramp will be at double rate. Choose multiple of 0.2!

%% Triggers
targcodes = [8]; % Continuous Dots
resp_string = {'12','13'};
%% EEG parameters
fs = 500; % sample rate
nchan = 97;
LPFcutoff_35Hz=35;
[H_LP_35Hz,G_LP_35Hz]=butter(4,(LPFcutoff_35Hz*2)/fs,'low');
LPFcutoff_8Hz=8;
[H_LP_8Hz,G_LP_8Hz]=butter(4,(LPFcutoff_8Hz*2)/fs,'low');

%% Time
ts = -0.3*fs:1.7*fs;
t = ts*1000/fs;

BLint = [-100 0];   % baseline interval in ms
default_response_time = 1.7-0.2;

%% Artifact rejection
ARchans = [1:97];

artifth = 100;
artifchans=[]; % keep track of channels on which the threshold is exceeded, causing trial rejection
frontal_chans = [1,2,65,66,33,34,35,36,67,3,37,4,38,5,39,6,40,7,68,41,42,45,46];
ch_CPP = [53];
ch_lr = [26,59,87;24,56,84];
ch_rl = [24,56,84;26,59,87];

artif_trials=[];

erp = []; erp_CSD = []; erp = single(erp); erp_CSD = single(erp_CSD); ET_trials = [];

allRT=[]; allrespLR=[]; allITI=[]; allUD=[]; saccadeT={}; % note allRT will be in time
numtr = 0; eye_artif_trials=0;

chanlocs = readlocs('JBhead96_sym.loc');
for f=1:length(files)    
    disp(f)
    filename=[files{f}];
    EEG = pop_loadbv(paths{f},files{f});
    loadbvSK_DN_JBhead_97
    EEG = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers. -88s are boundaries
    load([matfiles{f}])

    if ~exist(ET_matfiles{f}, 'file') %DN: if ET matfile NOT been saved
        FixEyelinkMessages_UD %then calculate and save it now
    end
    
%     for i = 1:97
%         ep = squeeze(EEG.data(i,:)); % time
%         nfft = size(ep,2);
%         fftx = abs(fft(ep,[],2))./(nfft/2);
%         fftx = fftx(:,1:ceil((nfft+1)/2));
%         freq_temp = (0:ceil((nfft+1)/2)-1)*fs/nfft;
%         fftx_all(i,:) = fftx;
%     end
%     figure, plot(freq_temp,mean(fftx_all,1))

    % There are peaks at 25, 50, 75, 100, 180, 200. Oddly a suppression at 60.
    % Get rid of mains frequency
    EEG = pop_eegfiltnew(EEG,59,61,[],1,0,0,0);
    EEG = pop_eegfiltnew(EEG,119,121,[],1,0,0,0);
    EEG = pop_eegfiltnew(EEG,179,181,[],1,0,0,0);
    % HP Filter
    EEG = pop_eegfiltnew(EEG,0.1,0,[]); % filter to 0.1Hz 
    
    % interpolate bad channels
    if ~isempty(badchans)
        EEG=eeg_interp(EEG,[badchans],'spherical');
    end    

    % LP Filter
    EEG = pop_eegfiltnew(EEG,0,8,[]);
    
    % average-reference the whole continuous data (safe to do this now after interpolation):
    EEG = pop_reref(EEG,[]);
    
    numev = length(EEG.event);
    % Fish out the event triggers and times
    clear trigs stimes
    for i=1:numev
        trigs(i)=EEG.event(i).type;
        stimes(i)=round(EEG.event(i).latency);
    end
    stimes(find(ismember(trigs,[-88,1,5])))=[];
    trigs(find(ismember(trigs,[-88,1,5])))=[];

    % Eyetracker stuff
    load(ET_matfiles{f}) %DN: load the ET mat file
    EEG = pop_importeyetracker(EEG,ET_matfiles{f},[first_event ...
        last_event],[1:4] ,{'TIME' 'GAZE_X' 'GAZE_Y' 'AREA'},0,1,0,0);
    [output_cell,~,~] = command_window_text();
    text = output_cell{length(output_cell)-1};
    numsamp = sscanf(text,'%*s%*s%*s%*s%*s%*f%*s%*s%f');
    if numsamp<20
        disp([allsubj{s},', block ',num2str(f),': ET sync issue'])
        keyboard
    end
    % Baseline position?
    temp = EEG.data([nchan+2,nchan+3],stimes(1):stimes(end));
    ET_data_mean(1) = mean(temp(1,find(temp(1,:)>1)),2); ET_data_mean(2) = mean(temp(2,find(temp(2,:)>1)),2);
    EEG.data(nchan+2,find(EEG.data(nchan+2,:)>0)) = EEG.data(nchan+2,find(EEG.data(nchan+2,:)>0))-ET_data_mean(1)+400;
    EEG.data(nchan+3,find(EEG.data(nchan+3,:)>0)) = EEG.data(nchan+3,find(EEG.data(nchan+3,:)>0))-ET_data_mean(2)+300;
    
    % get epochs
    EEG = pop_epoch(EEG,{'8'},[-0.3 1.702]);
    
    % Eyetracker stuff
    % scres = 800 x 600: 400,300 is middle. 75 pixels is 3 degrees.
    [EEG,Indexes] = pop_eegthresh(EEG,1,[nchan+2 nchan+3],[400-75 300-75],[400+75 300+75],EEG.times(1),EEG.times(end),0,1);
    eye_artif_trials = eye_artif_trials+length(Indexes);
    
    if length(Indexes)>18 % if loads of trials are missing, then it's a crap block and the eyetracker may have been badly set up.
        continue
    end
    plotter = 0;
    EEG = pop_detecteyemovements(EEG,[],[nchan+2 nchan+3],6,4,0.0405,1,0,25,3,plotter,1,0);
%     pause(0.5)
%     EEG = pop_detecteyemovements(EEG,left_eye_xy,right_eye_xy,vfac,mindur,...
%         degperpixel,smooth,globalthresh,clusterdist,clustermode,
%     plotfig,writesac,writefix)

    for n = 1:EEG.trials
        ep = squeeze(EEG.data(1:nchan,:,n));
        ep_ET = squeeze(EEG.data([nchan+1:nchan+4],:,n));
        ind2 = find(ismember(EEG.epoch(n).eventtype,{'12'}));
        if length(ind2)==1, response_time = str2num(EEG.epoch(n).eventtype{1,ind2}); else response_time = length(EEG.times)-0.1*fs; end
            
%         artifchans_thistrial = ARchans(find(max(abs(ep(ARchans,find(t>=-100 & ts<=(response_time+0.1*fs)))),[],2)>artifth));
        artifchans_thistrial = ARchans(find(max(abs(ep(ARchans,find(ts<=(response_time+0.1*fs)))),[],2)>artifth));
        artifchans = [artifchans artifchans_thistrial];
        
        numtr=numtr+1;
        if length(artifchans_thistrial) > 0, artif_trials(numtr) = 0; else artif_trials(numtr) = 1; end
        erp(:,:,numtr) = ep;
        ep_CSD = CSD(ep,G_CSD,H_CSD);
        erp_CSD(:,:,numtr) = ep_CSD;
        ET_trials(:,:,numtr) = ep_ET(:,:);
         
        allITI(numtr) = trialITI(n);
        allUD(numtr) = par.UPorDOWN;
        
        ind1 = find(ismember(EEG.epoch(n).eventtype,{'8'}));
%         if length(ind1)>1, ind1 = ind1(end); end
        
        ind2 = find(ismember(EEG.epoch(n).eventtype,{'12'}));
        if length(ind2)~=1 % | length(ind1)>1
            allrespLR(numtr) = 3;
            allRT(numtr) = 0;
        else
            resp = str2num(EEG.epoch(n).eventtype{1,ind2});
            allrespLR(numtr) = 1;
            allRT(numtr) = EEG.epoch(n).eventlatency{1,ind2}-EEG.epoch(n).eventlatency{1,ind1};
        end
        ind3 = find(ismember(EEG.epoch(n).eventtype,'saccade'));
        if isempty(ind3)
            % need 9 parts
            saccadeT{numtr} = zeros(1,9);
        else
            for i = 1:length(ind3)
                saccadeT{numtr}(i,:) = [EEG.epoch(n).eventlatency{1,ind3(i)},EEG.epoch(n).eventduration{1,ind3(i)},EEG.epoch(n).eventsac_amplitude{1,ind3(i)}, ...
                    EEG.epoch(n).eventsac_angle{1,ind3(i)},EEG.epoch(n).eventsac_endpos_x{1,ind3(i)},EEG.epoch(n).eventsac_endpos_y{1,ind3(i)}, ...
                    EEG.epoch(n).eventsac_startpos_x{1,ind3(i)},EEG.epoch(n).eventsac_startpos_y{1,ind3(i)},EEG.epoch(n).eventsac_vmax{1,ind3(i)}];
            end
        end
    end
end
    
%% Plots
rejected_trials = length(find(artif_trials==0));
figure;
hist(artifchans,[1:nchan]); title([allsubj{s} ': ' num2str(rejected_trials) ' artifacts = ',num2str(round(100*(rejected_trials/length(allRT)))),'%']) % s from runafew
disp([allsubj{s},' number of trials: ',num2str(length(find(allRT)))])
disp([allsubj{s},' number of rejected eye trials: ',num2str(eye_artif_trials)])

[counts,centers] = hist(artifchans,[1:nchan]);
figure;
topoplot(counts,chanlocs,'plotchans',[1:nchan],'electrodes','numbers');
title(subject_folder{s})

% baseline_erp = mean(erp_CSD(:,find(t>=BLint(1) & t<=BLint(2)),:),2);
% erp_CSD = erp_CSD-repmat(baseline_erp,[1,size(erp_CSD,2),1]); % baseline full erp
% erp_mean = squeeze(mean(erp_CSD(:,:,find(artif_trials==1 & allRT>0)),3));
% figure
% plottopo(erp_mean(:,:),'chanlocs',chanlocs,'limits',[t(1) t(end) ...
%     min(min(min(erp_mean(:,:))))  max(max(max(erp_mean(:,:))))], ...
%     'title',[allsubj{s}],'ydir',1)
% t1 = 800; t2 = 1100;
% plot_mean = squeeze(mean(erp_mean(:,find(t>t1 & t<t2)),2));
% figure
% topoplot(plot_mean,chanlocs,'maplimits', ...
%     [min(plot_mean) max(plot_mean)], 'electrodes','numbers','plotchans',[1:97]);

% return
pause(1)
erp_8Hz = erp;
erp_CSD_8Hz = erp_CSD;
artif_trials_8Hz = artif_trials;
%%
save([path_temp subject_folder{s} '\' allsubj{s} '_epochs_UD'], ...
    'erp_8Hz','erp_CSD_8Hz','artif_trials_8Hz','-append')
return
%% Save
save([path_temp subject_folder{s} '\' allsubj{s} '_epochs_UD'], ...
    'erp','erp_CSD','ET_trials','saccadeT',...
    'allRT','allrespLR','allITI','allUD','artifchans','artif_trials','eye_artif_trials','t')
return;


%% Code graveyard
%     clear temp; c=1;
%     for i=1:length(EEG.event)
%         if ~ismember(EEG.event(i).type,[-88,1,5])
%             temp(c) = EEG.event(i);
%             c=c+1;
%         end
%     end
%     EEG.event = temp;
%     
%     numev = length(EEG.event);
%     % Fish out the event triggers and times
%     clear trigs stimes
%     for i=1:numev
%         trigs(i)=EEG.event(i).type;
%         stimes(i)=round(EEG.event(i).latency);
%     end
%     stimes(find(ismember(trigs,[-88,1,5])))=[];
%     trigs(find(ismember(trigs,[-88,1,5])))=[];


    
    
    
%     % get triggers
%     targ_on = find(trigs==par.CD_TG);
% %     targ_off = find(trigs==par.CD_TGOFF);
%     if ~isequal(length(targ_on),length(trialITI))
%         disp('Trigger discrepancy')
%         length(targ_on),length(trialITI)
%         keyboard
%     end
%     
%     % get epochs
%     for n=1:length(targ_on)
%         numtr=numtr+1;
%         locktime = stimes(targ_on(n));
%         try
%             if trigs(targ_on(n)+1)==12 || trigs(targ_on(n)+1)==13
%                 response_time = stimes(targ_on(n)+1)-locktime; % time in samples from beginning of motion to response.
%                 response_time = floor(response_time);
%                 if response_time>default_response_time*fs
%                     response_time = default_response_time*fs;
%                 end
%             else
%                 response_time = default_response_time*fs;
%             end
%         catch
%             disp('EEG ended too soon')
%             response_time = default_response_time*fs;
%         end
%         
%         ep_LPF_35Hz = EEG_LPF_35Hz.data(1:nchan,locktime+ts);   % chop out an epoch
%         ep_LPF_8Hz = EEG_LPF_8Hz.data(1:nchan,locktime+ts);   % chop out an epoch
%         
%         % baselining to whole epoch for art reject
%         BLamp = mean(ep_LPF_35Hz,2);
%         ep_LPF_35Hz_base = ep_LPF_35Hz - repmat(BLamp,[1,length(t)]); % baseline correction
%         
%         % baselining to whole epoch for art reject
%         BLamp = mean(ep_LPF_8Hz,2);
%         ep_LPF_8Hz_base = ep_LPF_8Hz - repmat(BLamp,[1,length(t)]); % baseline correction
%         
%         ep_test = [find(t>=BLint(1) & t<BLint(2))];
%         if isempty(ep_test)
%             disp('Empty epoch for art rejection')
%             keyboard
%         end
%         
%         % artifact rejection
%         ep_test = [find(ts>=-0.3*fs & ts<=response_time+0.2*fs)];
%         if isempty(ep_test)
%             disp('Empty epoch for art rejection2')
%             keyboard
%         end
%         
%         artifchans_thistrial = ARchans(find(max(abs(ep_LPF_35Hz_base(ARchans,find(ts>=-0.3*fs & ts<=response_time+0.2*fs))),[],2)>artifth));
%         artifchans = [artifchans artifchans_thistrial];
%         
%         if length(artifchans_thistrial) > 0
%             rejected_trials=rejected_trials+1;
%             art_reject_trials(numtr) = 0;
%         else
%             art_reject_trials(numtr) = 1;
%         end
%             
%         ep_LPF_35Hz_CSD = CSD(ep_LPF_35Hz,G_CSD,H_CSD);
%         ep_LPF_8Hz_CSD = CSD(ep_LPF_8Hz,G_CSD,H_CSD);
%         
%         % assign epochs
%         erp_LPF_35Hz(:,:,numtr) = single(ep_LPF_35Hz);        
%         erp_LPF_8Hz(:,:,numtr) = single(ep_LPF_8Hz);
%         
%         % CSD epochs
%         erp_LPF_35Hz_CSD(:,:,numtr) = single(ep_LPF_35Hz_CSD);
%         erp_LPF_8Hz_CSD(:,:,numtr) = single(ep_LPF_8Hz_CSD);
%         
%         allITI(numtr) = trialITI(n);
%         allUD(numtr,:) = par.UPorDOWN;
%         
%         try % change this
%             if trigs(targ_on(n)+1)==12
%                 allrespLR(numtr) = 1;
%                 allRT(numtr) = stimes(targ_on(n)+1)-stimes(targ_on(n));
%             elseif trigs(targ_on(n)+1)==13
%                 allrespLR(numtr) = 2;
%                 allRT(numtr) = stimes(targ_on(n)+1)-stimes(targ_on(n));
%             else
%                 allrespLR(numtr) = 3;
%                 allRT(numtr) = 0;
%             end
%         catch
%             allrespLR(numtr) = 0;
%             allRT(numtr) = 0;
%         end
%     end