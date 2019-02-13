clear
close all
clc
pause(0.5);

path_temp = 'C:\Users\loughnge\Documents\SAIIT_UD\All Subjects\';

subject_folder = {'AA','AB','AC','AO','AV','DB','IV','JK','KB','MB',...
    'MH','NG','NL','RS','SB','SG','SH','SK','ST'};
allsubj = {'AA','AB','AC','AO','AV','DB','IV','JK','KB','MB',...
    'MH','NG','NL','RS','SB','SG','SH','SK','ST'};

allblocks = {{[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]},...
    {[2:6] [1:6]},{[2:6] [1:6]},{[1:5] [1:5]},{[1:6] [1:4,6]},{[1:6] [1:6]},...
    {[1:6] [1:6]},{[1:6] [1:6]},{[2:6] [2:6]},{[1:6] [1:6]},{[1:6] [1:6]},...
    {[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]}}; % Up then down

duds = [3,18]; % AC,SK, JK?? 3,18. IV chan 54 is shite
single_participants = []; % JK. SB, SH with chan 92?

allbadchans = {[36]...
    []...
    [8,89]...
    []...
    []... % 5
    []...
    [54]...
    [16,68]...
    [30]...
    [51,73]... % 10
    []...
    []...
    []...
    []...
    []... % 15
    []...
    [5,8,13,17,33,70,93]...
    []...
    []... % 19
    };

if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    allbadchans([duds]) = [];
    allblocks([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
    allbadchans = allbadchans([single_participants]);
    allblocks = allblocks([single_participants]);
end

% convert from locs to csd, don't need this now, already converted to csd.
% ConvertLocations('E:\Cued dots SK\recsdfilesfor96channelcap\JBhead96_sym.locs');
E = textread('C:\Program Files\MATLAB\R2015b\toolbox\CSDtoolbox\chans_JBhead_96_noground.asc','%s');
M = ExtractMontage('C:\Program Files\MATLAB\R2015b\toolbox\CSDtoolbox\JBhead96_sym.csd',E);  % reading in the montage for the CSD toolbox
MapMontage(M);
[G_CSD,H_CSD] = GetGH(M);

updown = {'U','D'};

h = waitbar(0,'Please wait...');
steps = length(allsubj);
step = 0;
for s=1:length(allsubj)
    tic
    if s==1
        waitbar(step/steps,h)
    else
        min_time = round((end_time*(steps-step))/60);
        sec_time = round(rem(end_time*(steps-step),60));
        waitbar(step/steps,h,[num2str(min_time),' minutes remaining'])
    end
    step=step+1;
    disp(['Subject: ',num2str(s)])
    disp(['Subject: ',allsubj{s}])
    
    badchans = allbadchans{s};
    k=1;
    clear files paths matfiles;
    for ud = 1:2
        blocks = allblocks{s}{ud};
        for bb = 1:length(blocks)
            paths{k} = [path_temp subject_folder{s} '\'];
            files{k} = [allsubj{s} '_' updown{ud} '_' num2str(blocks(bb)) '.vhdr'];
            matfiles{k} = [path_temp subject_folder{s} '\' allsubj{s} '_' updown{ud} '_' num2str(blocks(bb)) '.mat'];
            ET_files{k}=[path_temp subject_folder{s} '\' allsubj{s} '_' updown{ud} '_' num2str(blocks(bb)) '.asc'];
            ET_matfiles{k} = [path_temp subject_folder{s} '\' allsubj{s} '_' updown{ud} '_' num2str(blocks(bb)) '_ET.mat'];
            k=k+1;
        end
    end
    SAIIT_UD_analysis
    SAIIT_UD_analysis_ITI
        
    end_time = toc;
end
close(h)