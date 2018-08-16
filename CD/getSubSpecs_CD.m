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
% % getSubSpecs

subject_folder = {'AA','AB','AC','AO','AV','DB','IV','JK','KB','MB',...
    'MH','NG','NL','RS','SB','SG','SH','SK','ST'};
allsubj = {'AA','AB','AC','AO','AV','DB','IV','JK','KB','MB',...
    'MH','NG','NL','RS','SB','SG','SH','SK','ST'};

%%

TCD_bigdots = {''};
Monash_bigdots = subject_folder;

allParadigms    = {'CD'};
paradigms    = {'CD'};

% available data for {subject}
allblocks = {{[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]},...
    {[2:6] [1:6]},{[2:6] [1:6]},{[1:5] [1:5]},{[1:6] [1:4,6]},{[1:6] [1:6]},...
    {[1:6] [1:6]},{[1:6] [1:6]},{[2:6] [2:6]},{[1:6] [1:6]},{[1:6] [1:6]},...
    {[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]},{[1:6] [1:6]}}; % Up then down

% use function getChannelName to match number and name of channel
% allbadchans, organised in {paradigm}{[session]}
% {1} CDT
allbadchans = {[36]... % AA
    []...% AB
    [8,89]...% AC
    []...% AO
    []... % AV
    []...% DB
    [54]...% IV
    [16,68]...% JK
    [30]...%KB, ger only had 30
    [51,73]... % MB
    []... % MH
    []... % NG
    []... % NL
    []... % RS
    []... % SB
    []... % SG
    [5,8,13,17,33,70,93]... % SH
    []... % SK
    []... % ST
    };



