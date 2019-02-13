function [V]=varianceInflationCoefficients(X)
%vif() computes variance inflation coefficients  
%VIFs are also the diagonal elements of the inverse of the correlation matrix [1], a convenient result that eliminates the need to set up the various regressions
%[1] Belsley, D. A., E. Kuh, and R. E. Welsch. Regression Diagnostics. Hoboken, NJ: John Wiley & Sons, 1980.

% Variance inflation factor (VIF) quantifies how much the variance is
% inflated due to collinearity of regressor matrix columns. i_th entry in
% the output vector is the variance inflation factor for the i_th
% predictor, which indicates how much the variance of the i_th predictor is
% inflated due to collinearity.    
% https://uk.mathworks.com/matlabcentral/fileexchange/60551-vif-x

R0 = corrcoef(X); % correlation matrix
V=diag(inv(R0))';