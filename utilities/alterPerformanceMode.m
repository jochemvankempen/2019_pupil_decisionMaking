function alterPerformanceMode(goToState)
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
% change performance mode of PC
%
% goToState: string with options
%     - 'highPerformance'
%     - 'balanced'
%     - 'powerSaver'
%
% prevent PC from falling asleep.
% The GUID used in each case is different on each computer. It can be
% queried by typing 'powercfg -list' on a command line window. There is
% a GUID for each energy mode, just choose the ones that allow or
% prevent the computer from going to sleep mode.
% https://uk.mathworks.com/matlabcentral/fileexchange/36194-insomnia-prevent-computer-sleep-mode
%
% Jochem van Kempen, 2018
% 


recognizedPC = false;

switch upper(getComputerName)
    case 'FMS-ION-511027'
        recognizedPC    = true;
        highPerformance = '8c5e7fda-e8bf-4a96-9a85-a6e23a8c635c';
        balanced        = '381b4222-f694-41f0-9685-ff5bb260df2e';
        powerSaver      = 'a1841308-3541-4fab-bc81-f71556f20b4a';
end

if recognizedPC
    switch goToState
        case 'highPerformance'
            fprintf('switching to high performance mode, preventing sleep\n');
            system(['powercfg -setactive ' highPerformance]);
        case 'balanced'
            fprintf('switching to balanced mode\n');
            system(['powercfg -setactive ' balanced]);
        case 'powerSaver'
            fprintf('switching to powerSaver\n');
            system(['powercfg -setactive ' powerSaver]);
    end
end