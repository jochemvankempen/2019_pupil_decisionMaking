function [pstring,starstring] = getSignificanceStrings(p, rounded, addboldface, prestring)
% convert p value into string for plotting
% 
% INPUT
% p: p-value
% rounded: round p-value
% addboldface: make output string bold face if p-value is < 0.05
% prestring: string that will be placed before p-value
%
% OUTPUT
% pstring: string with p-value
% starstring: string with stars indicating level of significance
%
%
% Jochem van Kempen, 2017

pstring = cell(length(p),1);
starstring = cell(length(p),1);
for ip = 1:length(p)
    if p(ip) >= 0.05
        pstring{ip} = num2str(p(ip));
        starstring{ip} = '';
    elseif (p(ip) < 0.05) && (p(ip) >=0.01)
        pstring{ip} = '< 0.05';
        starstring{ip} = '*';
    elseif (p(ip) < 0.01) && (p(ip) >=0.001)
        pstring{ip} = '< 0.01';
        starstring{ip} = '**';
    elseif (p(ip) < 0.001)
        pstring{ip} = '< 0.001';
        starstring{ip} = '***';
    end
    
    if ~rounded
        pstring{ip} = [num2str(p(ip), '= %.3f')];
        
        if (p(ip) < 0.001)
            pstring{ip} = '< 0.001';
        end
    end
    
    if nargin == 4
        pstring{ip} = [prestring pstring{ip}];
    end
    
    if addboldface
        if p(ip) < 0.05
            pstring{ip} = ['\bf' pstring{ip}];
        end
    end
end

if length(pstring) == 1
    pstring = pstring{1};
    starstring = starstring{1};
end