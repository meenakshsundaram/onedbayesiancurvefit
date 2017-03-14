%
%    This file is a part of the onedbayesiancurvefit a matlab library.This 
%    is free software: you can redistribute it and/or modify it under 
%    the terms of the GNU Lesser General Public License as published by the 
%    Free Software Foundation, either version 3 of the License, or 
%    (at your option) any later version.
%
%    onedbayesiancurvefit is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License
%    along with rootfinding.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2012 Meenakshi Sundaram   
%
function [mn,Sn]=posteriormodel(Phi,y,beta,alpha)
% Estimates the samples for w the coefficients of the one dimensional 
% linear regression from
%
%  Input:
%       Phi     -   basis functions of the independent variable (nsamplesxnterms)
%       y       -   dependent variable   (nsamplesx1)
%       We are trying to get w such that 
%               y=w(1)+x*w(2)+x.*x*w(3)+...
%       beta    -   precision in evaluation of y
%       alpha   -   precision in the prior of w
%
%  Output:
%       mn      -   mean of the distribution of w (ntermsx1)
%       Sn      -   covariance of the distribution w (ntermsxnterms)
%
% This approach is taken from the book:
% "Pattern Recognition and Machine Learning - By Christopher Bishop"
% This computes the posterior for estimating the fitting parameters in 
% linear regression. It assumes a prior of a normal distribution centered 
% about 0 and carrying a precision of alpha. Also it assumes that the 
% precision of the dependent variable is known to be beta.

% nsamples
[~,nterms]=size(Phi);

% This formula is given in Page. 153, Eq. 3.50, Eq. 3.51 of the book
Sninv=alpha*eye(nterms)+beta*(Phi'*Phi);
Sn=inv(Sninv);
mn=beta*Sn*(Phi'*y);


