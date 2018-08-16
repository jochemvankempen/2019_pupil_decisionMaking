%%% unpack varargin into separate variables
%
% Jochem van Kempen. 19-06-2018

if ~isempty(varargin)
    if length(varargin{1}) > 2
        warning('separate arguments by ;')
    end
    for ifield = 1:length(varargin{1}(:,1))
        eval([varargin{1}{ifield,1} ' =  varargin{1}{ifield,2};'])
    end
end