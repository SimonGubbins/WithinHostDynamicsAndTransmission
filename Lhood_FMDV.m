function [logL, prior]=Lhood_FMDV(par,tClin,tStart,tStop,D,tObs,VI,mFlag)
%
% [logL, prior]=Lhood_FMDV(par,tClin,tStart,tStop,D,tObs,VI,mFlag)
%
% Matlab function for computing the log likelihood for a model linking
% within-host dynamics and transmission of FMDV in cattle based on data
% from a series of one-to-one challenge experiments
%
% Inputs:
% par - array containing the (transformed) model parameter
% tClin - time post infection (days) at which animals showed clinical signs
% tStart - time post infection (days) at which challenge began
% tStop - time post infection (days) at which challenge ended
% D - outcome of each challenge (0/1)
% tObs - observation times for virus isolation data
% VI - viral titres for each animal
% mFlag - flag indicating model to use for the probability of transmission
%         1 - proportional to titre
%         2 - proportional to log titre
%         3 - proportional to titre and transmission rate varies with
%             animal
%         4 - proportional to log titre and transmission rate varies with
%             animal
%
% Outputs:
% logL - log likelihood for parameters
% prior - log prior probability for parameters
 
%==========================================================================
% PREPARE THE INPUTS
% Determine the number of animals
nAnim=size(D,1);
 
% Extract the individual animal parameters (peak titre (Vp), time of peak
% titre (tp), growth rate (lg) and decay rate (ld))
Vp=exp(par(1:nAnim));
tp=par(nAnim+1:2*nAnim);
lg=par(2*nAnim+1:3*nAnim);
ld=par(3*nAnim+1:4*nAnim);
 
% Extract the hierarchical parameters for the within-host model
p_Vp=par(4*nAnim+1:4*nAnim+2);
p_tp=par(4*nAnim+3:4*nAnim+4);
p_lg=par(4*nAnim+5:4*nAnim+6);
p_ld=par(4*nAnim+7:4*nAnim+8);
 
% Extract the incubation period parameters
p_inc=par(4*nAnim+9:4*nAnim+10);
rho=par(4*nAnim+11);
 
% Extract the transmisson parameter(s)
if mFlag==1 || mFlag==2
    g=exp(par(end-1));
elseif mFlag==3
    g=exp(par(4*nAnim+12:5*nAnim+11));
    p_g=[par(end-2); par(end-1)];
elseif mFlag==4
    g=par(4*nAnim+12:5*nAnim+11);
    p_g=[par(end-2); par(end-1)];
end

% Extract the standard deviation for the error
sigE=par(end);
%==========================================================================
 
%==========================================================================
% COMPUTE THE PRIOR
% Individual animal parameters
prior=sum(log(gampdf(log(Vp),p_Vp(1),p_Vp(2)./p_Vp(1))))+...
      sum(log(lognpdf(tp,p_tp(1),p_tp(2))))+...
      sum(log(gampdf(lg,p_lg(1),p_lg(2)./p_lg(1))))+...
      sum(log(gampdf(ld,p_ld(1),p_ld(2)./p_ld(1))));
 
% Hierarchical parameters
prior=prior+log(exppdf(p_Vp(1),100))+...
            log(exppdf(p_Vp(2),100))+...
            log(exppdf(p_tp(1),100))+...
            log(exppdf(p_tp(2),100))+...
            log(exppdf(p_lg(1),100))+...
            log(exppdf(p_lg(2),100))+...
            log(exppdf(p_ld(1),100))+...
            log(exppdf(p_ld(2),100));

% Informative priors for the incubation  period parameters, except the
% correlation
prior=prior+log(gampdf(p_inc(1),5,1.6531/5))+...
            log(gampdf(p_inc(2),3,0.4724/3))+...
            log(unifpdf(rho,-1,1));
 
% Non-informative for the transmission parameter(s)
if mFlag==1
    prior=prior+log(exppdf(g,1e-4));
elseif mFlag==2
    prior=prior+log(exppdf(g,1));
elseif mFlag==3
    prior=prior+sum(log(normpdf(log(g),p_g(1),p_g(2))))+...
                log(normpdf(p_g(1),-9,1))+...
                log(exppdf(p_g(2),1));
elseif mFlag==4
    prior=prior+sum(log(gampdf(g,p_g(1),p_g(2)./p_g(1))))+...
                log(exppdf(p_g(1),1))+...
                log(exppdf(p_g(2),1));
end
 
% Informative for standard deviation of the error
prior=prior+log(gampdf(sigE,50,0.3./50));
%==========================================================================
 
%==========================================================================
% COMPUTE THE LOG LIKELIHOOD
% Initialise the log likelihood
logL=0;
 
% For each animal ...
for a=1:nAnim
    
% Set the observations to use
    x=(~isnan(VI(:,a)));
 
% Compute the expected viral titre at each observation
    muVI=2.*Vp(a)./...
         (exp(-lg(a).*(tObs(x)-tp(a)))+exp(ld(a).*(tObs(x)-tp(a))));
    muVI=log10(muVI);
 
% Compute the log likelihood for the shedding curve fitted to the virus
% isolation data
    cens=(VI(x,a)<=0);
    logL=logL+sum((1-cens).*log(normpdf(VI(x,a),muVI,sigE))+...
                  cens.*log(normcdf(0,muVI,sigE)));
 
% Calculate the probability of transmission at each challenge
    dt=0.01;
    intV=zeros(size(find(~isnan(D(a,:)))));
    for j=1:length(intV)
        if ~isnan(tStart(a,j))
            t=tStart(a,j):dt:tStop(a,j);
            y=2.*Vp(a)./(exp(-lg(a).*(t-tp(a)))+exp(ld(a).*(t-tp(a))));
            if mFlag==1
                intV(j)=g.*trapz(t,y);
            elseif mFlag==2
                y=log(y);
                y(y<=0)=0;
                intV(j)=g.*trapz(t,y);
            elseif mFlag==3
                intV(j)=g(a).*trapz(t,y);
            elseif mFlag==4
                y=log(y);
                y(y<=0)=0;
                intV(j)=g(a).*trapz(t,y);
            end
        end
    end
    PrT=1-exp(-intV);
 
% Compute the contribution to the log-likelihood for the challenge outcomes
    logL=logL+sum(log(binopdf(D(a,~isnan(D(a,:))),1,PrT)));
 
end
 
% Compute the contribution for the incubation period (allowing for joint
% distribution with time of peak viraemia)
mu=p_inc(1)+rho.*p_inc(2)./p_tp(2).*(log(max(0,tp))-p_tp(1));
sig=p_inc(2)*max(0,sqrt(1-rho.^2));
fC=(logncdf(tClin',mu,sig)-logncdf(tClin'-1,mu,sig));
 
% Add it to the log likelihood
logL=logL+sum(log(fC));
%==========================================================================
