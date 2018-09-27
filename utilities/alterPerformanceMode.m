function alterPerformanceMode(goToState)
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