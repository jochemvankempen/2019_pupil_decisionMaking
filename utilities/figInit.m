function [figHandle, figSet] = figInit(type,figNum,varargin)
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
% initialises figure or axis. 
% Consistent initialisation is useful to keep e.g. fontsizes the same
% across figures (by keeping the width/height of figure consistent).
%
% INPUT
%     fignum: number the figure will have
%     varargin: optional arguments
%         - figureType: 'Manuscript' (default), 'Poster' (will lead to different font sizes and line width, etc.)
%         - width: decimal number, width of figure will be a fraction of A4
%         - height: see width
%         - margin: 
%         - fontsize: fontsize to use for this figure/ax
%
%
% OUTPUT
% fighandle
% figSet: settings with linewidth, fontsize, etc
%
% example:
% to initialise figure, call figInit('fig', figNum, varargin)
% to initialise axes (e.g. subplot), call figInit('ax')
%
% Jochem van Kempen, 2017
 

getVarargin

if ~exist('figureType','var') 
    figureType = 'Manuscript';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% default settings
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% figure
defSet.units        = 'centimeters';
defSet.paperType    = 'a4letter';
        
defSet.FontName     = 'Times New Roman';
% defSet.FontName     = 'Arial';
defSet.transarency  = 0.3;

switch figureType
    case 'Manuscript'
        
        
        %%% axes
        defSet.LineWidth_ax     = 1.5;
        
        %%% font
        defSet.Fontsize_ax      = 11; % ticks
        defSet.Fontsize_ax_in   = 6; % ticks inset
        defSet.Fontsize_text    = 12; % labels/legends
        defSet.Fontsize_text_in = 10; % inside plot text
        defSet.Fontsize_title   = 14; % labels/legends
        defSet.Fontsize_panel   = 16; % panel label, e.g. Figure 1 panel A
        
        %%% plots
        defSet.LineWidth      = 2.5;% main plots
        defSet.LineWidth_in   = 1.5;% insets
        defSet.MarkerSize       = 8;
        
    case 'Poster'

        %%% axes
        defSet.LineWidth_ax     = 2;
        
        %%% font
        defSet.Fontsize_ax      = 10; % ticks
        defSet.Fontsize_ax_in   = 6; % ticks inset
        defSet.Fontsize_text    = 10; % labels/legends
        defSet.Fontsize_text_in = 8; % inside plot text
        defSet.Fontsize_title   = 14; % labels/legends
        defSet.Fontsize_panel   = 16; % panel label, e.g. Figure 1 panel A
        
        %%% plots
        defSet.LineWidth      = 2.5;% main plots
        defSet.LineWidth_in   = 1.5;% insets
        defSet.MarkerSize       = 8;

end

if ~exist('set','var') 
    figSet = defSet;
end
if ~exist('width','var')
    width = 10;
end
if ~exist('height','var')
    height = 10;
end

figHandle= [];
switch type
    case 'fig'
        figHandle = figure(figNum);clf
        set(figHandle, 'Units', figSet.units);
        
        % Select the preferred unit like inches, centimeters,
        % or pixels
        pos = get(figHandle, 'Position');
        pos(1) = 5;
        pos(2) = 5;
        pos(3) = width; % Select the width of the figure in [cm]
        pos(4) = height; % Select the height of the figure in [cm]
        set(figHandle, 'Position', pos);
        set(figHandle, 'PaperType', figSet.paperType);
        
        set(figHandle,'PaperPositionMode','auto')

        % Select the default font and font size
        % Note: Matlab does internally round the font size
        % to decimal pt values
        set(figHandle, 'DefaultTextFontSize', figSet.Fontsize_text); % [pt]
        set(figHandle, 'DefaultAxesFontSize', figSet.Fontsize_ax); % [pt]
        set(figHandle, 'DefaultAxesFontName', figSet.FontName);
        set(figHandle, 'DefaultTextFontName', figSet.FontName);
%         set(figHandle, 'DefaultAxesFontName', 'Arial');
%         set(figHandle, 'DefaultTextFontName', 'Arial');
        
        figSet.colors = get(gca,'ColorOrder');%

    case 'ax'
        % set(gcf,'color','none')
        set(gca,'linewidth',figSet.LineWidth_ax)
%         set(gca,'fontsize',fontsize)
        set(gca,'color','none')
        set(figHandle, 'DefaultTextFontSize', figSet.Fontsize_text); % [pt]
        set(figHandle, 'DefaultAxesFontSize', figSet.Fontsize_ax); % [pt]
        set(figHandle, 'DefaultAxesFontName', figSet.FontName);
        set(figHandle, 'DefaultTextFontName', figSet.FontName);
        
end
