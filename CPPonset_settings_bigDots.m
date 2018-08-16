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
% Set exceptions. For some subjects/bins, no CPP onset can be found, or the
% wrong one is found

clear win_mean_change consecutive_windows

win_mean_change = 0;
consecutive_windows=15;%15 works well for everybody else

switch bin2use
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% baseline pupil %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case {'pupil_lp_baseline_regress_iti_side'}  %
        if (sideInfo == 1 && nbin2use == 5)
            
            if any(strcmp(subject_folder{itmpsub},{'NT_16_04_14'})) || any(strcmp(subject_folder{itmpsub},{'OS_09_05_14'})) || any(strcmp(subject_folder{itmpsub},{'091M_SW'})) || ...
                    any(strcmp(subject_folder{itmpsub},{'279M_FT'})) || any(strcmp(subject_folder{itmpsub},{'114M_CS'}))
                consecutive_windows=50;%had to make it longer for these participants otherwise it records a false CPP onset
            end
            
            if any(strcmp(subject_folder{itmpsub},{'AC_13_05_14'})) || any(strcmp(subject_folder{itmpsub},{'037M_JD'}))
                win_mean_change = -1;
            end            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% average response locked Pupil %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
    case 'pupil_lp_RT_neg200_200_regress_bl_iti_side'
        
        if (sideInfo == 1 && nbin2use == 5)
            if any(strcmp(subject_folder{itmpsub},{'RM_06_05_14'})) || any(strcmp(subject_folder{itmpsub},{'289M_AS'})) || any(strcmp(subject_folder{itmpsub},{'108M_CY'})) || ...
                    any(strcmp(subject_folder{itmpsub},{'400M_ED'}))
                consecutive_windows=50;%had to make it longer for these participants otherwise it records a false CPP onset
            end
            if (any(strcmp(subject_folder{itmpsub},{'AS_23_04_14'})) && ibin==2) || (any(strcmp(subject_folder{itmpsub},{'OM_07_05_14'})) && ibin == 5)
                consecutive_windows=50;%had to make it longer for these participants otherwise it records a false CPP onset
            end
            if  any(strcmp(subject_folder{itmpsub},{'118M_CS'}))
                win_mean_change = -1;
            end
        end
        
end
