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
function [mtest,Stest]=predictivemodel(Phi,y,beta,alpha,Phitest)
% Estimates the prediction on the test data for one dimensional linear
% regression from 
%
%  Input:
%       Phi     -   basis functions of the independent variable for training data (nsamplesxnterms)
%       y       -   dependent variable   (nsamplesx1)
%       We are trying to get w such that 
%               y=w(1)+x*w(2)+x.*x*w(3)+...
%       beta    -   precision in evaluation of y
%       alpha   -   precision in the prior of w
%       Phitest -   basis functions of the independent variable for test data (nsamplesxnterms)
%
%  Output:
%       mtest   -   mean of the distribution of dependent variable for test data (ntermsx1)
%       Stest   -   covariance of the distribution of dependent variable for test data (ntermsxnterms)
%
% This approach is taken from the book:
% "Pattern Recognition and Machine Learning - By Christopher Bishop"
% This computes the predictive distribution of the linear regression for the 
% test data. It assumes a prior of a normal distribution centered about 0 
% and carrying a precision of alpha. Also it assumes that the precision of 
% the dependent variable is known to be beta.


% nsamples
[~,nterms]=size(Phi);

% This formula is given in Page. 153, Eq. 3.50, Eq. 3.51 of the book
Sninv=alpha*eye(nterms)+beta*(Phi'*Phi);
mn=beta*(Sninv\(Phi'*y));

% This formula is given in Page. 156, Eq. 3.58, Eq. 3.59 of the book
mtest=Phitest*mn;
Stest=diag(diag(Phitest*(Sninv\Phitest'))+(1./beta));
