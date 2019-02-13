% SAIIT UP/DOWN - detect gradual increases OR decreases in contrast (separate blocks)
% Contents:
% 1) load monitor params, set task parameters
% 2) Get filename
% 3) set up triggers and open PTB window
% 4) make stimuli
% 5) make stimulus sequences, save as vectors of textures
% 6) designate trigger codes
% 7) make randomized sequence of trial types
% 8) present instructions to subject
% 9) Start trials

clear;
close all
SITE = 'T';     % T = TCD, C = City college, E = EGI in City College
commandwindow;

load parUD                  % load last parameters (in structure called "par")
if SITE=='C'|SITE=='E'
    TheUsualParamsCRT_Dell_lores      % this script defines some useful parameters of the monitor, trigger port, etc common to all experiments
    load gammafnDell_lores100   % load the gamma function parameters for this monitor - other options are gammafn (LCD?) gammafnCRT (flames)
    par.BGcolor=midgray;
elseif SITE == 'T'
    TheUsualParamsCRT_Laptop
    load gammafnCRT
    % Gamma function for TCD???
    par.BGcolor=midgray;
end

%%%%%%%%% IMPORTANT SETTINGS
par.videoFrate = 60;   % Monitor refresh rate. Original was 100Hz
par.FlickF = 20;      % Flicker frequencies of two stimuli in Hz

par.numtargets = 24;
par.secs_btw_targs = [4 7 10];
par.spatfreq = 1;       % Spatial frequency of gratings
par.outerrad_deg = 8;   % in DEGREES
par.innerrad_deg = 3;   % in DEGREES
par.targrampdur = 1.2;  % in sec. Return ramp will be at double rate. Choose multiple of 0.2!
par.BLcontrast = 0.7;    % contrast
par.targChange = 0.3; % Max 0.5 for complete disappearance of other grating
par.useEL = 0;  % use the eye tracker?
par.eyeFBK = 0; %sound feedback, use only for trial block

% Other Settings
par.leadintime = 1000; % how long to pause before experiment starts

% ASK EXPERIMENTER TO INPUT SESSION/BLOCK (LOGFILE) NAME
dlg_title = 'SAIIT UD';
while 1
    prompt = {'Enter SUBJECT/RUN/TASK IDENTIFIER:','EEG? (1=yes, 0=no)','UP (1) OR DOWN (2)?'};
    def = {par.runID,num2str(par.recordEEG),num2str(par.UPorDOWN)};
    answer = inputdlg(prompt,dlg_title,1,def);
    par.runID = answer{1};
    par.recordEEG = str2num(answer{2});
    par.UPorDOWN = str2num(answer{3});
    if exist([par.runID '.mat'],'file'), 
        dlg_title = [par.runID '.mat EXISTS ALREADY - CHOOSE ANOTHER, OR DELETE THAT ONE IF IT IS RUBBISH']
    else
        break;
    end
end
tic

% SOUND STUFF
Fs = 22050; % Hz
High = 0.3*sin(2*pi*500*[0:1/Fs:0.1]);
si = hanning(Fs/100)';
env = [si(1:round(length(si)/2)) ones(1,length(High)-2*round(length(si)/2)) fliplr(si(1:round(length(si)/2)))];
hHigh = audioplayer(High.*env, Fs);

Low = 0.4*sin(2*pi*200*[0:1/Fs:0.3]);
si = hanning(Fs/100)';
env = [si(1:round(length(si)/2)) ones(1,length(Low)-2*round(length(si)/2)) fliplr(si(1:round(length(si)/2)))];
hLow = audioplayer(Low.*env, Fs);

% Set up for triggers
if par.recordEEG
    if SITE=='C'
        % USB port (posing as Serial Port) for triggers
        [port, errmsg] = IOPort('OpenSerialPort', 'COM4','BaudRate=115200');
        IOPort('Write', port, uint8([setpulsedur 2 0 0 0]))   % pulse width given by last 4 numbers (each a byte, little-endian)
    elseif SITE=='E'
        port = hex2dec('1030');
        lptwrite(port,0);
    elseif SITE=='T'
        % Parallel Port for triggers - set to zero to start
        port = 888;
        lptwrite(port,0);
    end
end

if par.useEL, ELCalibrateDialog, end

if par.useEL
    %%%%%%%%% EYETRACKING PARAMETERS
    par.FixWinSize = 3;    % RADIUS of fixation (circular) window in degrees
    par.TgWinSize = 3;    % RADIUS of fixation (circular) window in degrees
    ELsetupCalib
    Eyelink('Command', 'clear_screen 0')
    Eyelink('command', 'draw_box %d %d %d %d 15', center(1)-deg2px*par.FixWinSize, center(2)-deg2px*par.FixWinSize, center(1)+deg2px*par.FixWinSize, center(2)+deg2px*par.FixWinSize);
end

%  **********************  MAKE STIMULI
Rout = round(deg2px*par.outerrad_deg);  % radii in pix
Rin = round(deg2px*par.innerrad_deg);
D=Rout*2+1;                             % full stimulus size "D"
% First make the complete stimulus:
numchecks = [64 6];
[x,y] = meshgrid([1:D]-(D+1)/2,[1:D]-(D+1)/2);
[th,r]=cart2pol(x,-y);
th(find(th<0)) = th(find(th<0))+2*pi;
th1 = 0; th2 = 2*pi;
A = sin(numchecks(1)*pi.*(th-th1)./(th2-th1)) .* sin(numchecks(2)*pi.*(r-Rin)./(Rout-Rin));

% figure
% imshow(A(:,:,:))

midLum = GrayLevel2Lum(par.BGcolor,Cg,gam,b0);   % The very middle luminance on the monitor in cd/m^2
lumAmpl = floor(midLum);   % luminance amplitude (divergence from midLum) in cd/m^2
A(find(A>0))=lumAmpl; A(find(A<0))=-lumAmpl;    % convert sinusoidal luminance modulation of the spatial pattern to square wave

% figure
% imshow(A(:,:,:))

A(find(r>Rout | r<Rin)) = 0;
% A(Rout:Rout+2,Rout:Rout+2)=lumAmpl;
figure
imshow(mat2gray(A(:,:,:)))

% Now we'll make frame sequences, which comprise just a vector of multipliers for the pattern stimuli we've generated
% above (mostly 0 and 1 for off and on, and -1 for reversed pattern)
framesperflickercycle = round(par.videoFrate./par.FlickF);
BLframeseq = []; TGframeseq = [];

% baseline frame sequence length (this will be repeated again and again in the ITI)
BLframeseqlen = LCM_SK(framesperflickercycle);

% A standard baseline frame sequence:
% for f=1:length(par.FlickF)
%     ONframes = floor(framesperflickercycle(f)/2);
%     BLframeseq(:,f) = repmat([ones(1,ONframes) zeros(1,framesperflickercycle(f)-ONframes)],1,BLframeseqlen/framesperflickercycle(f))';
% end
BLframeseq=[1:3]';

% make the baseline frame sequence
B=A;
B(find((th>pi/2 & th<13*2*pi/32)))=0;
B(find((th>17*2*pi/32 & th<24*2*pi/32) | (th>27*2*pi/32 & th<31*2*pi/32))) = 0;
S(:,:,1) = B;
B=A;
B(find((th>3*2*pi/32 & th<8*2*pi/32) | (th>17*2*pi/32 & th<21*2*pi/32) | (th>24*2*pi/32 & th<31*2*pi/32))) = 0;
S(:,:,2) = B;
B=A;
B(find((th>0 & th<8*2*pi/32) | (th>13*2*pi/32 & th<17*2*pi/32) | (th>24*2*pi/32 & th<27*2*pi/32) | (th>31*2*pi/32 & th<32*2*pi/32))) = 0;
S(:,:,3) = B;
% B=A;
% B(find((th>0 & th<3*2*pi/32) | (th>8*2*pi/32 & th<17*2*pi/32) | (th>21*2*pi/32 & th<24*2*pi/32) | (th>31*2*pi/32 & th<32*2*pi/32))) = 0;
% S(:,:,4) = B;

figure
imshow(mat2gray(S(:,:,:)))

% Opens a graphics window on the main monitor
window = Screen('OpenWindow', whichScreen, par.BGcolor);

if abs(hz-par.videoFrate)>1
    error(['The monitor is NOT SET to the desired frame rate of ' num2str(par.videoFrate) ' Hz. Change it.'])
end

for n=1:size(BLframeseq,1)
    stim = midLum + S(:,:,BLframeseq(n))*par.BLcontrast;
    % Fixation point
    stim(Rout:Rout+2,Rout:Rout+2)=GrayLevel2Lum(255,Cg,gam,b0);
    figure
    imshow(mat2gray(stim(:,:,:)))
    clear Screen
    return
    BLstim(n) = Screen('MakeTexture', window, round(Lum2GrayLevel(stim,Cg,gam,b0)));
end
% I = mat2gray(stim(:,:,:));
% image(tempo(:,:,:))

stimrect = round([-1 -1 1 1]*D/2);

% make targets:
TGframeseqlen = floor(par.targrampdur*1.5*par.videoFrate/BLframeseqlen)*BLframeseqlen;  % Target frame sequence length, just the ramp down
rampdown = round(TGframeseqlen*2/3); rampup = TGframeseqlen - rampdown;
if par.UPorDOWN==1
    CIF=[par.BLcontrast+[1:rampdown]*par.targChange/rampdown fliplr(par.BLcontrast+[1:rampup]*par.targChange/rampup)];
else
    CIF=[par.BLcontrast-[1:rampdown]*par.targChange/rampdown fliplr(par.BLcontrast-[1:rampup]*par.targChange/rampup)];
end
% for f=1:length(par.FlickF)
%     ONframes = floor(framesperflickercycle(f)/2);
%     TGframeseq(:,f) = repmat([ones(1,ONframes) zeros(1,framesperflickercycle(f)-ONframes)],1,TGframeseqlen/framesperflickercycle(f))';
% end
TGframeseq = repmat([1:3]',[floor(TGframeseqlen/framesperflickercycle),1]);

for n=1:size(TGframeseq,1)
    stim = midLum + (CIF(n))*S(:,:,TGframeseq(n,1));
     % Fixation point
    stim(Rout:Rout+2,Rout:Rout+2)=GrayLevel2Lum(255,Cg,gam,b0);
    targstim(n,1) = Screen('MakeTexture', window, Lum2GrayLevel(stim,Cg,gam,b0));
end

% for n=1:size(TGframeseq,1)
%     Screen('DrawTexture', window, targstim(n), [], [center center] + stimrect);
%     Screen('Flip', window);
% end
% return

%  ************************************************* CODES AND TRIAL SEQUENCE
% trigger codes - can only use these 15: [1 4 5 8 9 12 13 16 17 20 21 24 25 28 29]
par.CD_RESP  = 1;
par.CD_FIXON = 4;
par.CD_TGOFF = 5;   % target off
par.CD_TG = 8;   % target   % one for each target type
par.CD_BUTTONS = [12 13];
par.CD_BEEP = 29;

% TRIAL SEQUENCE RANDOMIZATION
% Factors varying trial to trial: ITI (3) x Tilt (2, left/right)
% first make the smallest block of trial types that cover all possibilities:
block = [1:3];%[ones(1,length(par.secs_btw_targs)) 2*ones(1,length(par.secs_btw_targs)) ; ...
%         repmat(1:length(par.secs_btw_targs),1,2)];
% Then repeat that smallest block enough times to get the desired number of trials:
temp = repmat(block,[1,ceil(par.numtargets/size(block,2))]);
temp = temp(:,randperm(size(temp,2)));  % shuffle
trialITI = par.secs_btw_targs(temp(1,:));   % in seconds
% trialLR = temp(1,:);

% *********************************************************************************** START TASK
% Instructions:
Screen('DrawText', window, 'Fixate on the central dot.', 0.05*scres(1), 0.25*scres(2), 255);
if par.UPorDOWN==1
    Screen('DrawText', window, 'Press left button with right hand as soon as you are', 0.05*scres(1), 0.35*scres(2), 255);
    Screen('DrawText', window, 'sure that the stimulus has INCREASED in contrast.', 0.05*scres(1), 0.45*scres(2), 255);
else
    Screen('DrawText', window, 'Press left button with right hand as soon as you are', 0.05*scres(1), 0.35*scres(2), 255);
    Screen('DrawText', window, 'sure that the stimulus has DECREASED in contrast.', 0.05*scres(1), 0.45*scres(2), 255);
end



% Things 
Screen('DrawText', window, 'Press to begin', 0.05*scres(1), 0.65*scres(2), 255);
Screen('Flip', window); %that we'll save on a trial by trial basis
clear ITIstartT TargOnT RespLR RespT
numResp=1;

% Waits for the user to press a button before starting
[clicks,x,y,whichButton] = GetClicks(whichScreen,0);
if par.recordEEG, sendtrigger(par.CD_RESP,port,SITE,0), end
if par.useEL, Eyelink('Message', ['TASK_START']); end
RespT(1) = GetSecs;
RespLR(1) = whichButton;  if RespLR(numResp)==3, RespLR(numResp)=2; end  % The first response will be the one that sets the task going, after subject reads instructions

%%%%%%%%%%%%%%%%%%%% START TRIALS

% initial lead-in:
Screen('FillRect',window, 255, fixRect);
Screen('Flip', window);
WaitSecs(par.leadintime/1000);

imageArray = {};
% Start Task:
portUP=0; lastTTL=0; ButtonDown=0;
for n=1:par.numtargets
    % First show standard during ITI - the baseline frame sequence
    for m=1:round(trialITI(n)*par.videoFrate/BLframeseqlen)
        for s=1:BLframeseqlen
            if par.recordEEG, if SITE=='T'|SITE=='E', if portUP & GetSecs-lastTTL>0.01, lptwrite(port,0); portUP=0; end, end, end
            Screen('DrawTexture', window, BLstim(s), [], [center center] + stimrect);
            if m==1 & s==1
                % Whole bunch of triggers
                if par.recordEEG, sendtrigger(par.CD_TGOFF,port,SITE,1); portUP=1; lastTTL=GetSecs; end
                if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) 'TGOFF' num2str(par.CD_TGOFF)]); end
    %             Screen('FillRect',window, 255, syncRect);
                [VBLTimestamp ITIstartT(n)] = Screen('Flip', window);
            else
                
                
                Screen('Flip', window);
%                 imageArray = [imageArray; {Screen('GetImage', window)}];
            end
            checkButton
            if par.eyeFBK
                checkeyeSK
                if isnan(x) & isnan(y) % blink
                    play(hHigh)
                    if par.recordEEG, sendtrigger(par.CD_BEEP,port,SITE,1); portUP=1; lastTTL=GetSecs; end
                elseif sqrt(x^2+y^2)>deg2px*par.FixWinSize
                    play(hLow)
                    if par.recordEEG, sendtrigger(par.CD_BEEP,port,SITE,1); portUP=1; lastTTL=GetSecs; end
                end
            end
        end
    end
    % present target
    for s=1:TGframeseqlen
        if par.recordEEG, if SITE=='T'|SITE=='E', if portUP & GetSecs-lastTTL>0.01, lptwrite(port,0); portUP=0; end, end, end
        Screen('DrawTexture', window, targstim(s,1), [], [center center] + stimrect);
        if s==1
            if par.recordEEG, sendtrigger(par.CD_TG,port,SITE,1); portUP=1; lastTTL=GetSecs; end
            if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) 'TG' num2str(par.CD_TG)]); end
            [VBLTimestamp TargOnT(n)] = Screen('Flip', window);
        else
            Screen('Flip', window);
        end
        checkButton
        if par.eyeFBK
            checkeyeSK
            if isnan(x) & isnan(y) % blink
                play(hHigh)
                if par.recordEEG, sendtrigger(par.CD_BEEP,port,SITE,1); portUP=1; lastTTL=GetSecs; end
            elseif sqrt(x^2+y^2)>deg2px*par.FixWinSize
                play(hLow)
                if par.recordEEG, sendtrigger(par.CD_BEEP,port,SITE,1); portUP=1; lastTTL=GetSecs; end
            end
        end
    end
end
% Lead-out
for m=1:round(par.secs_btw_targs(1)*par.videoFrate/BLframeseqlen) % shortest ITI...
    for s=1:BLframeseqlen
        if par.recordEEG, if SITE=='T'|SITE=='E', if portUP & GetSecs-lastTTL>0.01, lptwrite(port,0); portUP=0; end, end, end
        Screen('DrawTexture', window, BLstim(s), [], [center center] + stimrect);
        if m==1 & s==1
            if par.recordEEG, sendtrigger(par.CD_TGOFF,port,SITE,1); portUP=1; lastTTL=GetSecs; end
            if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) 'TGOFF' num2str(par.CD_TGOFF)]); end
            %             Screen('FillRect',window, 255, syncRect);
            [VBLTimestamp ITIstartT(par.numtargets+1)] = Screen('Flip', window);
        else
            Screen('Flip', window);
        end
        checkButton
        if par.eyeFBK
            checkeyeSK
            if isnan(x) & isnan(y) % blink
                play(hHigh)
                if par.recordEEG, sendtrigger(par.CD_BEEP,port,SITE,1); portUP=1; lastTTL=GetSecs; end
            elseif sqrt(x^2+y^2)>deg2px*par.FixWinSize
                play(hLow)
                if par.recordEEG, sendtrigger(par.CD_BEEP,port,SITE,1); portUP=1; lastTTL=GetSecs; end
            end
        end
    end
end

% I = mat2gray(imageArray{1});
% image(I)

if par.useEL, 
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    ELdownloadDataFile
end
cleanup
toc 
save([par.runID],'ITIstartT','TargOnT','RespT','RespLR','trialITI','par') 
save parUD par
