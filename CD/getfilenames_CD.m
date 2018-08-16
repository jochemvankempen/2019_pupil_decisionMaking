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
% 
% get subject file names, and change them from standard where necessary
% Change filenames to NaN if file is missing (based on xls notes of data)

% set the directories where data should be loaded/saved
switch upper(getComputerName)
    case 'FMS-ION-511027'
        paths.base  = ['E:\Monash\CD\'];
end


paths.data  = [paths.base 'data' filesep]; %
paths.pop   = [paths.data 'population' filesep]; %
clear files
for isub2=1:length(subject_folder)
    paths.s(isub2).base           = [paths.data subject_folder{isub2} filesep];
    paths.s(isub2).raw            = [paths.s(isub2).base 'raw' filesep];
    paths.s(isub2).fig            = [paths.s(isub2).base 'fig_new' filesep];
    
    for iside = 1:2 % up then down
        for iblock = allblocks{isub2}{iside}
            
            switch iside
                case 1
                    fileExtDir = 'U';
                case 2
                    fileExtDir = 'D';
            end
            files(isub2,iside,iblock).subID       = subject_folder{isub2};
            
            
            files(isub2,iside,iblock).edf   = [allsubj{isub2} '_' fileExtDir '_' num2str(iblock) '.edf'];
            files(isub2,iside,iblock).eeg   = [allsubj{isub2} '_' fileExtDir '_' num2str(iblock) '.eeg'];
            files(isub2,iside,iblock).mat   = [allsubj{isub2} '_' fileExtDir '_' num2str(iblock) '.mat'];
            files(isub2,iside,iblock).vhdr  = [allsubj{isub2} '_' fileExtDir '_' num2str(iblock) '.vhdr'];
            files(isub2,iside,iblock).vmrk  = [allsubj{isub2} '_' fileExtDir '_' num2str(iblock) '.vmrk'];
            
            files(isub2,iside,iblock).ET_files     = [paths.s(isub2).raw allsubj{isub2} '_' fileExtDir '_' num2str(iblock) '.asc'];
            files(isub2,iside,iblock).ET_matfiles  = [paths.s(isub2).raw allsubj{isub2} '_' fileExtDir '_' num2str(iblock) '_ET.mat'];
            
           
        end
    end
end


