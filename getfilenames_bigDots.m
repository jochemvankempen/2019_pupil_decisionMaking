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
%
% get subject file names, and change them from standard where necessary
% Change filenames to NaN if file is missing (based on xls notes of data)


% set the directories where data should be loaded/saved
switch upper(getComputerName)
    case 'FMS-ION-511027'
        paths.base  = ['E:\Monash\bigDots\'];
    case 'FMS-ION-ALEX02'
        paths.base  = ['D:\Monash\bigDots\'];  
    case 'ALEX1'
        paths.base  = ['E:\Monash\bigDots\'];  
end

paths.server  = ['\\campus\rdw\ion13\13\thieletapes3\Jochem\Monash\bigDots\'];

if exist('readDataFromServer','var') && readDataFromServer
    paths.readdata  = [paths.server 'data' filesep]; %
    paths.savedata  = [paths.base 'data' filesep]; %

else
    paths.readdata  = [paths.base 'data' filesep]; %
    paths.savedata  = [paths.base 'data' filesep]; %
end

paths.pop   = [paths.savedata 'population' filesep]; %


for isub2=1:length(subject_folder)
    % path_temp = 'C:\Users\Ger Loughnane\Documents\Main Files\PhD\Projects\Evidence Acculumation Project\Dots Analysis\Study Participants\';
    % path_temp = 'C:\Users\Ger Loughnane\Documents\Main Files\PhD\Projects\Evidence Acculumation Project\Dots Analysis\Study Participants\';
    paths.s(isub2).readbase       = [paths.readdata subject_folder{isub2} filesep];
    paths.s(isub2).savebase       = [paths.savedata subject_folder{isub2} filesep];
    paths.s(isub2).raw            = [paths.s(isub2).readbase 'raw' filesep];
    paths.s(isub2).fig            = [paths.s(isub2).readbase 'fig_new' filesep];
    
    if ~exist(paths.s(isub2).savebase, 'dir')
        mkdir(paths.s(isub2).savebase)
    end
    
    for iblock = allblocks{isub2}
        
        files(isub2,iblock).subID       = subject_folder{isub2};
        
        if ismember(subject_folder{isub2},TCD_bigdots)
           
            files(isub2,iblock).bdf   = [allsubj{isub2} '_' num2str(iblock) '.bdf'];
            files(isub2,iblock).edf   = [allsubj{isub2} '_' num2str(iblock) '.edf'];
            files(isub2,iblock).mat   = [allsubj{isub2} '_' num2str(iblock) '.mat'];
                        
            files(isub2,iblock).ET_files     = [paths.readdata 'Samples_and_Events' filesep allsubj{isub2} '_' num2str(iblock) '.asc'];
            files(isub2,iblock).ET_matfiles  = [paths.s(isub2).readbase allsubj{isub2} '_' num2str(iblock) '_ET.mat'];

        elseif ismember(subject_folder{isub2},Monash_bigdots)
            files(isub2,iblock).edf   = [allsubj{isub2} num2str(iblock) '.edf'];
            files(isub2,iblock).eeg   = [allsubj{isub2} num2str(iblock) '.eeg'];
            files(isub2,iblock).mat   = [allsubj{isub2} num2str(iblock) '.mat'];
            files(isub2,iblock).vhdr  = [allsubj{isub2} num2str(iblock) '.vhdr'];
            files(isub2,iblock).vmrk  = [allsubj{isub2} num2str(iblock) '.vmrk'];
            
            files(isub2,iblock).ET_files     = [paths.readdata 'Samples_and_Events' filesep allsubj{isub2} num2str(iblock) '.asc'];
            files(isub2,iblock).ET_matfiles  = [paths.s(isub2).readbase allsubj{isub2} num2str(iblock) '_ET.mat'];
                        
        end
        
        switch subject_folder{isub2}
            case 'PR_20_04_14'
                if iblock == 6
                    files(isub2,iblock).ET_files     = [paths.data 'Samples_and_Events' filesep allsubj{isub2} '_' num2str(iblock) 'b.asc'];
                end
        end
    end
end


