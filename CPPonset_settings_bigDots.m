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
        %%% GLM estimated
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'GLM_pupil_Ramp_stim_regress_iti_side'
        if (sideInfo == 1 && nbin2use == 5)
            if any(strcmp(subject_folder{itmpsub},{'NT_16_04_14'})) || any(strcmp(subject_folder{itmpsub},{'SH_25_05_14'})) || any(strcmp(subject_folder{itmpsub},{'059M_HP'})) || ...
                    any(strcmp(subject_folder{itmpsub},{'103M_JK'})) || any(strcmp(subject_folder{itmpsub},{'091M_SW'})) || any(strcmp(subject_folder{itmpsub},{'191M_DM'})) || ...
                    any(strcmp(subject_folder{itmpsub},{'352M_MK'})) || any(strcmp(subject_folder{itmpsub},{'414M_LA'}))
                consecutive_windows=50;%had to make it longer for these participants otherwise it records a false CPP onset
            end
            if (any(strcmp(subject_folder{itmpsub},{'114M_CS'}))  && ibin==5)
                consecutive_windows=100;%had to make it longer for these participants otherwise it records a false CPP onset
            end
        end
    case 'GLM_pupil_Ramp_stim_regress_bl_iti_side'
        if (sideInfo == 1 && nbin2use == 5)
            if any(strcmp(subject_folder{itmpsub},{'NT_16_04_14'})) || any(strcmp(subject_folder{itmpsub},{'SH_25_05_14'})) || any(strcmp(subject_folder{itmpsub},{'059M_HP'})) || ...
                    any(strcmp(subject_folder{itmpsub},{'093M_BR'})) || any(strcmp(subject_folder{itmpsub},{'091M_SW'})) || any(strcmp(subject_folder{itmpsub},{'331M_CL'})) || ...
                    any(strcmp(subject_folder{itmpsub},{'352M_MK'})) || any(strcmp(subject_folder{itmpsub},{'114M_CS'}))
                consecutive_windows=50;%had to make it longer for these participants otherwise it records a false CPP onset
            end
        end     
    case 'GLM_pupil_Ramp_stim_regress_bl_blPhase_iti_side'
        if (sideInfo == 1 && nbin2use == 5)
            if any(strcmp(subject_folder{itmpsub},{'MH_14_04_14'})) || any(strcmp(subject_folder{itmpsub},{'SB_08_05_14'})) || any(strcmp(subject_folder{itmpsub},{'279M_FT'}))
                consecutive_windows=50;%had to make it longer for these participants otherwise it records a false CPP onset
            end
            if ( (any(strcmp(subject_folder{itmpsub},{'059M_HP'})) || any(strcmp(subject_folder{itmpsub},{'037M_JD'})))  && ibin==2)
                consecutive_windows=50;%had to make it longer for these participants otherwise it records a false CPP onset
            end
            if (any(strcmp(subject_folder{itmpsub},{'114M_CS'}))  && (ibin==5))
                consecutive_windows=100;%had to make it longer for these participants otherwise it records a false CPP onset
            end
            if any(strcmp(subject_folder{itmpsub},{'ND_16_05_14'}))
                win_mean_change = -0.5;
            end
            if any(strcmp(subject_folder{itmpsub},{'037M_JD'})) || any(strcmp(subject_folder{itmpsub},{'384M_PD'}))
                win_mean_change = -1;
            end
        end
        
        
        
        
        
end
