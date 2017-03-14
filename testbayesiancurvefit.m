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
%% clear
clear;clc;

%% Generate data
nsamples=1028;

% precision parameters 
beta=11.1; % y
alpha=1e-3; % prior of x

% Setting the seed
rng(14329874)

% Training data
x=sort(rand(nsamples,1)*2.0*pi);
y=randn(nsamples,1)*2.0/(sqrt(beta))+cos(x);
% Test data
xtest=sort(rand(nsamples,1)*2.0*pi);

% Nterms 
nterms=4;

% Centers of the Gaussian basis functions
muj=linspace(min(x),max(x),nterms);

% Gaussian basis functions standard deviation - length scale parameter
s=2;

% Training/Test Basis functions
Phitest=ones(nsamples,nterms);
Phi=ones(nsamples,nterms);
for iind=1:nterms
    % Gaussian basis function
    Phi(:,iind)=exp(-(x-muj(iind)).^2/(2.0*s*s));
    Phitest(:,iind)=exp(-(xtest-muj(iind)).^2/(2.0*s*s));
end

% Plot training data
figure(1);clf;set(gca,'fontsize',32,'linewidth',2);
scatter(x,y,'o','linewidth',1);hold on;box on;
figure(2);clf;set(gca,'fontsize',32,'linewidth',2);
scatter(x,y,'o','linewidth',1);hold on;box on;

%% Bayesian estimation of the fitting parameters
[mn,Sn]=posteriormodel(Phi,y,beta,alpha);

% 10 samples of ws from the multivariate normal distribution
wsamples=mvnrnd(mn',Sn,10)';

% Plot fit on test data
figure(1);
plot(xtest,Phitest*wsamples,'linewidth',4);
axis tight;axis square;

%% Bayesian estimation of test data
[mtest,Stest]=predictivemodel(Phi,y,beta,alpha,Phitest);

% Plot fit on test data
figure(2);
ytest=mvnrnd(mtest',Stest,nsamples)';
plot(xtest,mtest,'linewidth',4);
plot(xtest,mtest+diag(Stest),'linewidth',4);
plot(xtest,mtest-diag(Stest),'linewidth',4);
axis tight;axis square;
