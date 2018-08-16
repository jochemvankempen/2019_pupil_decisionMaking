function [ET_data, artTrial] = interpolateET(ET_data, targtrigs, ET_event, blinkevent, fs, filename, dataset, artifChan, artifET, RT, t_artif)
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
% find blinks, interpolate them, and filter data
% also, checks whether blinks occured during/around the trial and outputs
% trial indices with artifacts.

% ET_event = events from pop_importeyetracker
% blinkevent = manually extracted events from fixEyelinkMessages.m

%
plotFigures = 1;
plotVisible = 'off';

blinkstring = {'L_blink','R_blink'};
padding = .200; %how long around a blink do we want to interpolate. In seconds


%% events from pop_importeyetracker
for i=1:length(ET_event)
    if (strcmpi(ET_event(i).type,blinkstring{1}) || strcmpi(ET_event(i).type,blinkstring{2}))
        ET_trigs(i) = 99;
    else
        try
            ET_trigs(i) = str2double(ET_event(i).type);
        catch
            ET_trigs(i) = NaN;
        end
    end
    ET_stimes(i)    = round(ET_event(i).latency);
end

ET_events = [ET_stimes(~isnan(ET_trigs) & (ET_trigs>1) & (ET_trigs<107))' ET_trigs(~isnan(ET_trigs) & (ET_trigs>1) & (ET_trigs<107))'];

switch dataset
    case {'bigDots'}
        ET_targettrigsIdx   = (ET_trigs>100 & ET_trigs<107);%trial start
        ET_targettrigs      = ET_trigs(ET_targettrigsIdx);
        ET_targettrigsTimes = ET_stimes(ET_targettrigsIdx);
        
        if ET_targettrigs(end)==ET_targettrigs(1)
            ET_targettrigs      = ET_targettrigs(1:end-1); % GL: indices of trigs when motion on. get rid of last trig, it was a repeat
            ET_targettrigsTimes = ET_targettrigsTimes(1:end-1); %
        end
    case 'CD'
        ET_targettrigsIdx   = (ET_trigs == 8);%trial start
        ET_targettrigs      = ET_trigs(ET_targettrigsIdx);
        ET_targettrigsTimes = ET_stimes(ET_targettrigsIdx);
end


ET_blinkTimes = ET_stimes(ET_trigs==99);%blinks

%% sync times ET_event and blinkevent (calculated in fixEyelinkMessages)
switch dataset
    case {'bigDots'}
        eyetrigsIdx = find(blinkevent(:,2)>100 & blinkevent(:,2)<107);
        eyetrigs    = blinkevent(eyetrigsIdx,2);
    case 'CD'
        eyetrigsIdx = find(blinkevent(:,2) == 8);
        eyetrigs    = blinkevent(eyetrigsIdx,2);
end

% first compare blinkevent with targtriggers
for itrig = 1:size(targtrigs,1)
    if find(eyetrigs(itrig:end)==targtrigs(itrig,2),1,'first')
        x(itrig) = find(eyetrigs(itrig:end)==targtrigs(itrig,2),1,'first');
    end
    trigIdx(itrig) = itrig + x(itrig)-1;
end

blinktimes2sync = blinkevent(eyetrigsIdx(trigIdx),1);
newblinktimes   = round(interp1(blinktimes2sync,targtrigs(:,1), blinkevent(:,1),'linear','extrap'));

newblinkevent = [newblinktimes blinkevent(:,2)];

%perform check
if sum(ET_targettrigs(1:length(targtrigs(:,2))) == targtrigs(:,2)') ~= length(trigIdx)
    keyboard
end
if sum(ET_targettrigsTimes(1:length(targtrigs(:,2))) == targtrigs(:,1)') ~= length(trigIdx)
    keyboard
end


%% find blinks and determine the times to interpolate
pupilDat = double(ET_data(4,:));

time2interp =[-padding padding]*fs; % data window being interpolated, from start blink - time2interp(1) to end blink + time2interp(2). 

%initialise
blinktimes      = zeros(length(ET_blinkTimes),2);
blinkvector     = zeros(1,length(pupilDat));
interpvector    = zeros(1,length(pupilDat));

% loop over blink times to find start and end of blinks
for iblink = 1:length(ET_blinkTimes)
    
    blinktimes(iblink,1) = ET_blinkTimes(iblink);
    [~,blinkIdx] = min(abs(newblinkevent(:,1) - blinktimes(iblink,1)));
    
    if newblinkevent(blinkIdx+1,2) == 2
        blinktimes(iblink,2) = newblinkevent(blinkIdx+1,1);
    else
        disp('no end of blink found')
        twin = [-0.1 1000]*fs;
        twin2plot   = ET_blinkTimes(iblink) + twin;
        if twin2plot(2) <= length(pupilDat)
            dat2plot    = pupilDat(twin2plot(1):twin2plot(2));
        else
            dat2plot    = pupilDat(twin2plot(1):end);
        end
        t2plot      = twin2plot(1):twin2plot(2);
        [~,endBlink] = find(diff(dat2plot==0) == -1, 1, 'first');
        blinktimes(iblink,2) = t2plot(endBlink);
    end
    
    blinkvector((blinktimes(iblink,1) ) : (blinktimes(iblink,2) )) = 1;
    interpvector((blinktimes(iblink,1) + time2interp(1)) : (blinktimes(iblink,2) + time2interp(2))) = 1;
end

interptimes = reshape(find(diff(interpvector)), 2, length(find(diff(interpvector)))/2)';

% remove data to be interpolated
tOriginal = 1:length(pupilDat);
pupilDat(logical(interpvector)) = [];
t2interp = tOriginal;
t2interp(logical(interpvector)) = [];

% interpolate data
ET_data(5,:) = interp1(t2interp,pupilDat, tOriginal,'linear');

%% cut out period where we actually recorded
% start and end of pupil data is zeros, this messes with the filters

tmp = min(diff(ET_data(5,:)==0));
[~,startIdx] = find(diff(ET_data(5,:)==0) == tmp, 1, 'first'); % get start of  recording
tmp = max(diff(ET_data(5,:)==0));
[~,endIdx] = find(diff(ET_data(5,:)==0) == tmp, 1, 'last'); % get end of  recording

cut(1) = floor((startIdx + (ET_targettrigsTimes(1) - startIdx) / 2)); % take a point between start of recording and first trial
cut(2) = ceil((endIdx - (endIdx - ET_targettrigsTimes(end)) / 2));

tmpDat = double((ET_data(5,cut(1):cut(2)))); % cut out data

%prepare variables to put the filtered data back into
for i = 6:7
    ET_data(i,:) = zeros(size(ET_data(5,:)));
end

%%% filters
[b,a]      = butter(2, 0.01 / (fs / 2), 'high');
tmpDat_hp = filtfilt(b,a, double(tmpDat)); % filter with zero lag

[b,a] = butter(2, 6 / (fs / 2), 'low');
tmpDat_filt     = filtfilt(b,a, double(tmpDat)); % filter with zero lag
tmpDat_hp_filt  = filtfilt(b,a, double(tmpDat_hp)); % filter with zero lag

ET_data(6,cut(1):cut(2)) = tmpDat_filt; % put back in original matrix
ET_data(7,cut(1):cut(2)) = tmpDat_hp_filt; % 

%%% check whether blinks occured during specified time windows around the trial
if nargin == 11
    % initialise struct for artefact trials
    art_times = fields(t_artif);
    
    for ifield = 1:length(art_times)
        times2check=[];
        art_trials = [];
        eval(['art_time = t_artif.' art_times{ifield} ';'])
        
        switch art_times{ifield}
            case {'pretarg', 'iti'} % don't add the 2nd time to the RT for these fields
                times2check = repmat(ET_targettrigsTimes(1:length(RT))',1,2) + ([repmat(art_time(1), length(RT),1)       repmat(art_time(2),length(RT),1) ]);
            otherwise
                times2check = repmat(ET_targettrigsTimes(1:length(RT))',1,2) + ([repmat(art_time(1), length(RT),1)       RT' + art_time(2) ]);
        end
        
        for itrial = 1:size(times2check,1)

            for iblink = 1:size(blinktimes,1)
                if (blinktimes(iblink) >= times2check(itrial,1)) && (blinktimes(iblink) <= times2check(itrial,2))
                    art_trials = [art_trials itrial];
                end
            end
        end
        
        eval(['artTrial.' art_times{ifield} ' = unique(art_trials);'])
    end
    
    
    if plotFigures
        % plot figure with information on alignment, trial rejection, and filtering
        figure(1),clf
        set(gcf,'visible',plotVisible,'position',[40 502 1853 554])
        hold on, plot(ET_data(4,:),'linewidth',2),plot(ET_data(5,:),'linewidth',2), plot(ET_data(6,:),'linewidth',2), plot(ET_data(7,:),'linewidth',2);%, plot(ET_data(7,:))
        hTrial = plot(repmat(ET_targettrigsTimes,2,1), repmat(ylim',1,length(ET_targettrigsTimes)));
        YLIM = get(gca,'ylim');
        hblink = plot(ET_blinkTimes, repmat(YLIM(1),1,length(ET_blinkTimes)),'.','markersize',20,'color','k');
        artifacts = [artifChan + artifET];
        for iart = find(artifacts==0 | artifacts==1)
            set(hTrial(iart),'linewidth',3)
        end
        startTimes = repmat(ET_targettrigsTimes,2,1);
        hRT = plot(startTimes(:,1:length(RT))+repmat(RT,2,1), repmat(ylim',1,length(RT)));
        for itrial = 1:length(RT)
            set(hRT(itrial), 'Color', hTrial(itrial).Color)
        end
        
        for iart = artTrial.neg500_RT_200
            set(hTrial(iart),'Marker','.','markersize',30)
        end
        
        saveFigName = [filename(1:end-4) '_filt'];
        fileExt = 'jpeg';
        print(gcf,['-d' fileExt],[saveFigName '.' fileExt])
        
        savefig(gcf,[saveFigName])
        
        slashloc = findstr(filename,'\');
        figname = filename(slashloc(end)+1:end-4);
        saveFigName = ['E:\Monash\' dataset '\data\eyecheck\' figname];
        fileExt = 'jpeg';
        print(gcf,['-d' fileExt],[saveFigName '.' fileExt])
        
    end
    

end



