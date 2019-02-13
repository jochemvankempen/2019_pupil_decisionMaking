function figSave(figName, figPath, figFiletype, figHandle)
% Saves figure
% figSave(figName, paths.fig, figFiletype)
% Example input = figSave(figName, paths.fig, {'png', 'svg'})
%
% Jochem van Kempen, 2017

if ~exist(figPath,'dir')
    mkdir(figPath)
end

if ~iscell(figFiletype)
    figFiletype = {figFiletype};
end

if nargin<4
    figHandle = gcf;
end



for ifiletype = 1:length(figFiletype)
    
    savefilename = [figPath figName '.' figFiletype{ifiletype}];

    switch figFiletype{ifiletype}
        case 'svg'
            set(figHandle,'renderer','painters') % make sure all elements are saved as independent paths
            
            saveas(figHandle, [figPath figName '.svg'], 'svg');
            
%             plot2svg(savefilename, figHandle)
            
        case {'epsc', 'eps'}
            set(figHandle,'renderer','painters') % make sure all elements are saved as independent paths
            
            saveas(figHandle, [figPath figName '.eps'], 'epsc');
            
%             print(figHandle,['-d' 'epsc' ],[figPath figName '.eps'])

%             export_fig(savefilename, '-depsc')
%             export_fig([figPath figName '.eps'], '-depsc')

        case 'png'
            print(figHandle,['-d' figFiletype{ifiletype} ], savefilename)
            
        case 'pdf'
            
            export_fig(savefilename)
        otherwise
            saveas(figHandle, savefilename);
    end
end
