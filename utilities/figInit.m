function [figHandle, figSet] = figInit(type,figNum,varargin)
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
defSet.paperType    = 'a4';
defSet.paperSize    = [21.0 29.7];
        
switch figureType
    case 'Manuscript'
        
        
        %%% axes
        defSet.axLineWidth      = 2;
        defSet.axFontsize       = 8; % ticks
        defSet.axFontsize_in    = 6; % ticks inset
        defSet.axLabFontsize    = 10; % labels/legends
        defSet.axTitleFontsize  = 14; % labels/legends
        
        %%% plots
        defSet.plLineWidth      = 2.5;% main plots
        defSet.plLineWidth_in   = 1.5;% insets
        defSet.plFontsize       = 8; % inside plot text
        defSet.MarkerSize       = 8;
        
    case 'Poster'
%         defSet.paperType    = 'a4';
%         defSet.paperSize    = [118.9 84.1];
        
        %%% axes
        defSet.axLineWidth      = 3;
        defSet.axFontsize       = 12; % ticks
        defSet.axFontsize_in    = 10; % ticks inset
        defSet.axLabFontsize    = 14; % labels/legends
        
        %%% plots
        defSet.plLineWidth      = 3.5;% main plots
        defSet.plLineWidth_in   = 2.5;% insets
        defSet.plFontsize       = 12; % inside plot text
        defSet.MarkerSize       = 12;
        defSet.axTitleFontsize  = 18; % labels/legends
end

if ~exist('set','var') 
    figSet = defSet;
end
if ~exist('width','var')
    width = 1;
end
if ~exist('height','var')
    height = 1;
end
if ~exist('margin','var')
    margin = 2.5;
end
if ~exist('fontsize','var')
    fontsize = figSet.axFontsize;
end

defsize = figSet.paperSize;
defsize = defsize - (margin * 2);

width = width * defsize(1);
height = height * defsize(2);

left = (defsize(1)- width)/2 + 3;
bottom = (defsize(2)- height)/2 + 2;
defsize = [left, bottom, width, height];

figHandle= [];
switch type
    case 'fig'
        figHandle = figure(figNum);clf
        set(figHandle,'PaperType',figSet.paperType, 'PaperUnits',figSet.units, 'Units', figSet.units);
        set(figHandle,'Position', defsize);
        set(figHandle,'PaperPositionMode','auto')
        
        figSet.colors = get(gca,'ColorOrder');%

    case 'ax'
        % set(gcf,'color','none')
        set(gca,'linewidth',figSet.axLineWidth)
        set(gca,'fontsize',fontsize)
        set(gca,'color','none')
end
