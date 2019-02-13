function plot_subplot_label(axHandle, cfg, label, displace)
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
% Plot label for figure panel, e.g. Figure 1A
%
% Jochem van Kempen, 30-01-2019
if isfield(axHandle,'OuterPosition')
    % needs work
%     pos = (axHandle.Position - axHandle.OuterPosition);

%     pos = axHandle.OuterPosition;
    pos = axHandle.Position;
else
    pos = axHandle.Position;
end

if ~exist('displace','var') || isempty(displace)
    displace = [0.07 0.05];
end

if isnumeric(label)
    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    
    label = Alphabet(label);
end


dim = [pos(1)-displace(1), pos(2)+pos(4)-displace(2),0.1,0.1];
h = annotation('textbox',dim, 'String', label);
h.FontSize = cfg.Fontsize_panel;
h.FontName = cfg.FontName;
h.EdgeColor = 'None';
% 
