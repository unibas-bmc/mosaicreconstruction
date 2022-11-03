function [ lambda ] = lambda_from_E( E )
%LAMBDA_FROM_E Gives wavelength [m] from energy [eV]
%   lambda = hc/E
lambda = 1.239842*10^(-6)./E;

end

