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
 
This library provides the following algorithms for one dimensional curve fitting using a bayesian model:
                        y= phi(x)' w + Normal(0,beta^-1)
where, 
x      represents the independent variable
phi(x) represents the vector of basis functions
w      represents the vector of coefficents for the basis functions
       with a prior given by Normal(0,alpha^-1)
       alpha    represents the precision of the prior
y      represents the dependent variable
beta   represents the precision of the dependent variable


posteriormodel.m        -   multivariate Normal distribution for the vector w
predictivemodel.m       -   multivariate Normal distribution for the input xtest
testbayesiancurvefit    -   demo program
