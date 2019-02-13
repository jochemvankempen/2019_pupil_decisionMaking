function [outp, chan_overlap, chan_exlude] = getChannelName(inp, type)
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
% Function description
% convert channel name to number or vice versa
%
% INPUT
% inp: channel name or number (or vector/matrix/cell)
% type: string matching 'BP' or 'biosemi'
%
% OUTPUT
% outp: channel name or number
% chan_overlap: channels that are both present in BP and biosemi (optional)
% chan_exlude: channels that are not present in both BP and biosemi (optional)
%

%% brain products
% Channel locations acticap 65ch standard-2
% (1)FP1 , (2)FP2,
% (33)AF7, (34)AF3,  (GND)AFz, (35)AF4, (36)AF8
% (3)F7  , (37)F5 ,  (4)F3,    (38)F1,  (5)Fz,  (39)F2,   (6)F4,   (40)F6   (7)F8
% (41)FT9, (42)FT7,  (8)FC5,   (43)FC3, (9)FC1, (65)FCz, (10)FC2, (44)FC4, (11)FC6, (45)FT8, (46)FT10
% (12)T7 , (47)C5,   (13)C3,   (48)C1,  (14)Cz, (49)C2,   (15)C4,  (50)C6,  (16)T8
% (17)TP9, (51)TP7,  (18)CP5,  (52)CP3, (19)CP1,(53)CPz,  (20)CP2, (54)CP4, (21)CP6, (55)TP8, (22)TP10
% (23)P7,  (56)P5,   (24)P3,   (57)P1,  (25)Pz, (58)P2,   (26)P4,  (59)P6,  (27)P8
% (28)PO9, (60)PO7,  (61)PO3,  (62)POz, (63)PO4,(64)PO8,  (32) PO10
% (29)O1,  (30)Oz,   (31)O2
chan_BP = {...
    'FP1' ,'FP2' ,'F7' ,'F3' ,'Fz' ,...%1-5
    'F4'  ,'F8'  ,'FC5','FC1','FC2',...%6-10
    'FC6' ,'T7'  ,'C3' ,'Cz' ,'C4' ,...%11-15
    'T8'  ,'TP9' ,'CP5','CP1','CP2',...%16-20
    'CP6' ,'TP10','P7' ,'P3' ,'Pz' ,...%21-25
    'P4'  ,'P8'  ,'PO9','O1' ,'Oz' ,...%26-30
    'O2'  ,'PO10','AF7','AF3','AF4',...%31-35
    'AF8' ,'F5'  ,'F1' ,'F2' ,'F6' ,...%36-40
    'FT9' ,'FT7' ,'FC3','FC4','FT8',...%41-45
    'FT10','C5'  ,'C1' ,'C2' ,'C6' ,...%46-50
    'TP7' ,'CP3' ,'CPz','CP4','TP8',...%51-55
    'P5'  ,'P1'  ,'P2' ,'P6' ,'PO7',...%56-60
    'PO3' ,'POz' ,'PO4','PO8','FCz'...%61-65
    };
% %brain products
% chStr.left_hemi   = {...
%     'Fp1','AF7','AF3','F7','F5','F3','F1','FT9','FT7','FC5','FC3','FC1',...
%     'T7','C5','C3','C1','TP9','TP7','CP5','CP3','CP1',...
%     'P7','P5','P3','P1','PO9','PO7','PO3','O1'};
% chStr.right_hemi   = {...
%     'Fp2','AF4','AF8','F2','F4','F6','F8','FC2','FC4','FC6','FT8','FT10',...
%     'C2','C4','C6','T8','CP2','CP4','CP6','TP8','TP10',...
%     'P2','P4','P6','P8','PO4','PO8','O2','PO10'};
% chStr.centre_chans= {...
%     'Fz','FCz','Cz','CPz','Pz','POz','Oz'};
% chStr.elec_pairs  = {...
%     'FP1','FP2'; 'AF7','AF8'; 'AF3','AF4' ; 'F7','F8' ; 'F5','F6' ; 'F3','F4' ; 'F1','F2' ; ...
%     'FT9','FT10' ; 'FT7','FT8' ; 'FC5','FC6' ; 'FC3','FC4' ; 'FC1','FC2' ; ...
%     'T7','T8' ; 'C5','C6' ; 'C3','C4' ; 'C1','C2'; ...
%     'TP9','TP10';'TP7','TP8';'CP5','CP6';'CP3','CP4';'CP1','CP2';...
%     'P7','P8';'P5','P6';'P3','P4';'P1','P2';...
%     'PO9','PO10';'PO7','PO8';'PO3','PO4';'O1','O2'};


%% biosemi 64 channel (Monash)
chan_biosemi64 = {...
    'FP1'   , 'AF7' , 'AF3' , 'F1'  , 'F3'  , ...
    'F5'    , 'F7'  , 'FT7' , 'FC5' , 'FC3' , ...
    'FC1'   , 'C1'  , 'C3'  , 'C5'  , 'T7'  , ...
    'TP7'   , 'CP5' , 'CP3' , 'CP1' , 'P1'  , ...
    'P3'    , 'P5'  , 'P7'  , 'P9'  , 'PO7' , ...
    'PO3'   , 'O1'  , 'Iz'  , 'Oz'  , 'POz' , ...
    'Pz'    , 'CPz' , 'Fpz' , 'FP2' , 'AF8' , ...
    'AF4'   , 'AFz' , 'Fz'  , 'F2'  , 'F4'  , ...
    'F6'    , 'F8'  , 'FT8' , 'FC6' , 'FC4' , ...
    'FC2'   , 'FCz' , 'Cz'  , 'C2'  , 'C4'  , ...
    'C6'    , 'T8'  , 'TP8' , 'CP6' , 'CP4' , ...
    'CP2'   , 'P2'  , 'P4'  , 'P6'  , 'P8'  , ...
    'P10'   , 'PO8' , 'PO4' , 'O2'  ...
    };
%chan order compared to BP
%[1,34,7,5,38,40,42,9,11,46,44,15,13,48,50,52,17,19,56,54,23,21,31,58,60,27,29,64,2,3,36,35,6,4,39,41,8,10,45,43,14,12,49,51,16,18,32,55,53,22,20,57,59,25,26,30,63,62,47]

% left_hemi = [1:23,25:27];
% right_hemi = [34:36,39:46,49:60,62:64];
% centre_chans = [29:32,38,47,48];
% elec_pairs = [1,34;2,35;3,36;4,39;5,40;6,41;7,42; ...
%     8,43;9,44;10,45;11,46;12,49;13,50;14,51;15,52; ...
%     16,53;17,54;18,55;19,56;20,57;21,58;22,59;23,60; ...
%     25,62;26,63;27,64];
% chStr.left_hemi = {
%     'FP1'    'AF7'    'AF3'    'F1'    'F3'    'F5'    'F7'
%     'FT7'    'FC5'    'FC3'    'FC1'    'C1'    'C3'    'C5'
%     'T7'    'TP7'    'CP5'    'CP3'    'CP1'    'P1'    'P3'
%     'P5'    'P7'    'PO7'    'PO3'    'O1'};
% chStr.right_hemi = {
%     'FP2'    'AF8'    'AF4'    'F2'    'F4'    'F6'    'F8'
%     'FT8'    'FC6'    'FC4'    'FC2'    'C2'    'C4'    'C6'
%     'T8'    'TP8'    'CP6'    'CP4'    'CP2'    'P2'    'P4'
%     'P6'    'P8'    'PO8'    'PO4'    'O2'};
% chStr.centre_chans = {
%     'Oz'    'POz'    'Pz'    'CPz'    'Fz'    'FCz'    'Cz'};
% chStr.elec_pairs = {
%     'FP1' 'FP2' ; 'AF7' 'AF8' ; 'AF3' 'AF4' ; 'F1' 'F2' 
%     'F3'  'F4'  ; 'F5'  'F6'  ; 'F7'  'F8'  ; 'FT7' 'FT8'
%     'FC5' 'FC6' ; 'FC3' 'FC4' ; 'FC1' 'FC2' ; 'C1'  'C2' 
%     'C3'  'C4'  ; 'C5'  'C6'  ; 'T7'  'T8'  ; 'TP7' 'TP8'
%     'CP5' 'CP6' ; 'CP3' 'CP4' ; 'CP1' 'CP2' ; 'P1'  'P2' 
%     'P3'  'P4'  ; 'P5'  'P6'  ; 'P7'  'P8'  ; 'PO7' 'PO8'
%     'PO3' 'PO4' ; 'O1'  'O2'};

%% biosemi 96 channels

chan_biosemi96 = {...
    'Fp1', 'Fp2', 'F7', 'F3', 'Fz', ...
    'F4', 'F8', 'FC5', 'FC1', 'FC2', ...
    'FC6', 'T7', 'C3', 'Cz', 'C4', ...
    'T8', 'TP9', 'CP5', 'CP1', 'CP2', ...
    'CP6', 'TP10', 'P7', 'P3', 'Pz', ...
    'P4', 'P8', 'PO9', 'O1', 'Oz', ...
    'O2', 'PO10', 'AF7', 'AF3', 'AF4', ...
    'AF8', 'F5', 'F1', 'F2', 'F6', ...
    'FT9', 'FT7', 'FC3', 'FC4', 'FT8', ...
    'FT10', 'C5', 'C1', 'C2', 'C6', ...
    'TP7', 'CP3', 'CPz', 'CP4', 'TP8', ...
    'P5', 'P1', 'P2', 'P6', 'PO7', ...
    'PO3', 'POz', 'PO4', 'PO8', 'AFp1', ...
    'AFp2', 'F9', 'F10', 'FFC5h', 'FFC1h', ...
    'FFC2h', 'FFC6h', 'TPP7h', 'CPP5h', 'CPP3h', ...
    'CPP1h', 'CPP2h', 'CPP4h', 'CPP6h', 'TPP8h', ...
    'P9', 'P10', 'PPO9h', 'PPO5h', 'PPO1h', ...
    'PPO2h', 'PPO6h', 'PPO10h', 'POO9h', 'POO1', ...
    'POO2', 'POO10h', 'O9', 'OI1h', 'OI2h', ...
    'O10', 'REF_FCz', 'GND_AFz'};




%%
chan_exlude=[];
if nargin<1
    chan_overlap    = find(ismember(chan_biosemi64,chan_BP));
    chan_exlude     = find(~ismember(chan_biosemi64,chan_BP));
    outp=[];
    return
end

switch type
    case 'BP' %brain products
        cN = chan_BP;
        chan_overlap    = 1:length(chan_BP);
    case 'biosemi64'
        cN = chan_biosemi64;
%         chan_exlude     = find(~ismember(chan_biosemi64,chan_BP));
    case 'biosemi96'
        cN = chan_biosemi96;
end



if isnumeric(inp)
    outp = cN(inp);
else
    if iscell(inp)
        for iN = 1:size(inp,1)
            for iM = 1:size(inp,2)
                outp(iN,iM) = find(strcmpi(cN, inp(iN,iM)));
            end
        end
    else
        outp = find(strcmpi(cN, inp));
    end
end

