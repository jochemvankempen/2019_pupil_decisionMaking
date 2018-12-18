function lp = linearProjection(erp, tt)
% lp = linearProjection(erp, tt)
%
% computes the linear projection of single trials on the average vector.
% 
% INPUT
% erp:              time x trial matrix
% tt:               time values (optional). If given, performs check of erp matrix
%                   dimensions and is used for plotting
%
% OUTPUT
% lp: linear projection for each individual trial
%
% based on papers 10.1038/78856, 10.1073/pnas.1317557111 and 10.1111/ejn.12859
% jochem van kempen 22/02/2017


% check if erp matrix has correct number of dimension
if ndims(erp)>2
    error('linear projection only works on 1 channel')
end

% check if erp matrix has size (time * trial), if not transpose
if size(erp,1) ~= length(tt)
    disp('transposing erp matrix to fit time x trial')
    erp = erp';
end
[nTime, nTrial] = size(erp);

%%% compute average vector and the norm to base the projection on
averageVector   = mean(erp,2); %column vector, average across trials
normVector      = averageVector/norm(averageVector)^2;% norm of the average vector
    
lp = zeros(nTrial,1);
%%% Calculate linear projection
for itrial = 1:nTrial
    
    % for each trial, compute linear projection by multiplying single trial
    % row vector by the norm and deviding by the length (to compensate for
    % RT dependence)
    lp(itrial,1) = (erp(:,itrial)'*normVector)/length(normVector);

end
