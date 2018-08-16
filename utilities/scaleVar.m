function outvar = scaleVar(invar, method, dim)
% scales matrix of values according to method
%
% INPUT
% invar: variable to be scaled
% method: 
%     - 'max', scales by maximum value 
%     - 'minmax', scales to the minimum and maximum value, so that it
%         will be scaled between 0 and 1 according to the formula
%         yi=(xi-min(xi)) /(max(xi)-min(xi)) 
%     - 'maxdim', scales by maximum value within a given dimension   
%     - 'minmaxdim', scales according to minmax method within given
%         dimension.
% dim: optional. Dimension of matrix according to which to scale variable
%
% jochem van kempen, 20-03-2018

if nargin < 3
    dim = [];
end

switch method
    case 'max'
        outvar = invar / max(invar(:));
    case 'minmax'
        outvar = (invar - min(invar(:))) / (max(invar(:)) - min(invar(:)));
    case {'maxdim','minmaxdim'}
        formula = [];
        for idim = 1:ndims(invar)
            if idim == dim
                formula = [formula 'i' ];
            else
                formula = [formula ':' ];
            end
            if idim ~= ndims(invar)
                formula = [formula ',' ];
            end
        end
        
        outvar = zeros(size(invar));
        for i = 1:size(invar,dim)
            
            eval(['tmp = squeeze(invar(' formula '));'])
            switch method
                case 'maxdim'
                    eval(['outvar(' formula ') = tmp / max(tmp(:));']);
                case 'minmaxdim'
                    eval(['outvar(' formula ') = (tmp / max(tmp(:))) / (max(tmp(:)) - min(tmp(:)));']);
            end
        end
            
end