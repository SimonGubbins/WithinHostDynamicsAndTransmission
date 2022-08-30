function ParEst_FMDV(dtype,mFlag,seeds,nsamp,nburnin,nthin)
%
% ParEst_FMDV(dtype,mFlag,seeds,nsamp,nburnin,nthin)
%
% Matlab function to find implement a Bayesian MCMC scheme to estimate
% parameters linking within-host dynamics and transmission of FMDV in
% cattle using a simple phenomenological model.
%
% Inputs:
% dtype - string indicating proxy used for determining infectiousness:
%         'Blood', 'OPF' (oropharyngeal fluid); or 'NF' (nasal fluid)
% mFlag - flag indicating model to use for the probability of transmission
%         1 - proportional to titre
%         2 - proportional to log titre
%         3 - proportional to titre and transmission rate varies with
%             animal
%         4 - proportional to log titre and transmission rate varies with
%             animal
% seeds - vector of seeds (must be of length nchains) to use for the random
%         number generator for each chain
% nsamp - number of samples to use when estimating parameters
% nburnin - number of samples to discard before estimating parameters
% nthin - number of samples by which to thin each chain
%
% Outputs: none (N.B. All chains are saved rather than provided as output
% arguments)

% Set the number of chains
nchains=length(seeds);

%==========================================================================
% PREPARE THE DATA
% Load the challenge data
varload=load('FMDV_TransmissionExperimentData');
 
% Extract the data needed for the analyses
tStart=varload.tStart;
tStop=varload.tStop;
tClin=varload.tClin;
D=varload.D;
 
% Load the virus isolation data
varload=load('FMDV_VirusIsolationData');
 
% Extract the data needed for the analyses
tObs=varload.tObs;
if strcmp(dtype,'Blood')==1
    VI=varload.VI_B;
elseif strcmp(dtype,'NF')==1
    VI=varload.VI_NF;
elseif strcmp(dtype,'OPF')==1
    VI=varload.VI_OPF;
end
 
% Set the number of animals in the study
nAnim=size(D,1);
%==========================================================================
 
% Set the number of parameters
if mFlag==1 || mFlag==2
    npar=4*(nAnim+2)+3+1+1;
elseif mFlag==3 || mFlag==4
    npar=4*(nAnim+2)+3+(nAnim+2)+1;
end

% Create the arrays storing the output for the chain
ParSamp=cell(1,nchains);
 
% For each chain ...
parfor chain=1:nchains
 
%==========================================================================
% INITIALISE THE CHAIN
% Initialise the random number generator
    rng(seeds(chain),'twister');
 
% Set the initial scaling factor for the proposal distribution
    sf=(2.38.^2)/npar;
    SIG=eye(npar);
 
% Set the counter for the number of accepted samples
    n_accept=0;
 
% Create the arrays storing the output for the chain
    ParSampC=zeros(nsamp/nthin,npar+2);
    iter=1;
 
% Generate the initial parameters for the chain, ensuring they generate a
% finite log likelihood and prior
    disp('Initialising chain')
    CurrL=NaN;
    prior=NaN;
    while ~isfinite(CurrL+prior)
 
% Generate an initial set of infectiousness parameters
        Vp0=10.^(unifrnd(3,5,nAnim,1));
        tp0=unifrnd(1,3,nAnim,1);
        lg0=unifrnd(5,15,nAnim,1);
        ld0=unifrnd(1,5,nAnim,1);
 
% Compute the hierarchical parameters
        p=gamfit(log(Vp0));
        p_Vp0=[p(1); p(1).*p(2)];
        p_tp0=lognfit(tp0)';
        p=gamfit(lg0)';
        p_lg0=[p(1); p(1).*p(2)];
        p=gamfit(ld0)';
        p_ld0=[p(1); p(1).*p(2)];
 
% Set the incubation period parameters
        pinc0=lognfit(tClin)';
 
% Set the correlation between generation time and incubation period
        rho0=unifrnd(-1,1);

% Set the transmission parameter(s)
        if mFlag==1
            g0=unifrnd(1e-4,1e-3);
        elseif mFlag==2
            g0=unifrnd(1e-1,1);
        elseif mFlag==3
            g0=unifrnd(1e-4,1e-3,nAnim,1);
            [m,s]=normfit(log(g0));
            p_g0=[m; s];
        elseif mFlag==4
            g0=unifrnd(1e-1,1,nAnim,1);
            p=gamfit(g0)';
            p_g0=[p(1); p(1)*p(2)];
        end
 
% Set the error standard deviation for the viral titre curve
        sig0=unifrnd(0.1,0.5);
 
% Merge the initial parameter set
        if mFlag==1 || mFlag==2
            par=[log(Vp0); tp0; lg0; ld0; ...
                 p_Vp0; p_tp0; p_lg0; p_ld0; ...
                 pinc0; rho0; log(g0); sig0];
        elseif mFlag==3
            par=[log(Vp0); tp0; lg0; ld0; ...
                 p_Vp0; p_tp0; p_lg0; p_ld0; ...
                 pinc0; rho0; log(g0); p_g0; sig0];
        elseif mFlag==4
            par=[log(Vp0); tp0; lg0; ld0; ...
                 p_Vp0; p_tp0; p_lg0; p_ld0; ...
                 pinc0; rho0; g0; p_g0; sig0];
        end

% Compute the log-likelihood
        [CurrL, prior]=Lhood_FMDV(par,tClin,tStart,tStop,D,tObs,VI,mFlag);
 
    end
%==========================================================================
 
%==========================================================================
% UPDATE THE PARAMETERS
% Sample parameter space
    disp('Sampling parameter space')
    for samp=1:nsamp+nburnin
 
% Indicate what's going on
        disp(['Chain: ' num2str(chain) ', Sample: ' num2str(samp) ';'...
              ' Accept: ' num2str(100*n_accept/samp,3) '%'])
 
% Update the variance-covariance matrix for the proposal distribution
        if samp<=nburnin && (samp<=2*npar || n_accept==0)
            SIGp=0.01*eye(npar);
        else
            SIGp=sf.*(SIG+0.01*eye(npar));
        end
 
% Generate the new set of probabilities
        par_new=par+mvnrnd(zeros(1,length(par)),SIGp)';
 
% Compute the log likelihood and prior for the new parameter set
        [NewL,prior_new]=Lhood_FMDV(par_new,tClin,tStart,tStop,D,...
                                    tObs,VI,mFlag);
 
% Test whether to accept the new parameter set
        u=unifrnd(0,1);
        if isfinite(NewL+prior_new) && ...
           u<min(1,exp((NewL+prior_new)-(CurrL+prior)))
 
% Update the counter
            n_accept=n_accept+1;
 
% Update the covariance matrix for the proposal distribution
            if n_accept==1
                pbar=mean([par par_new],2);
                SIG=cov([par'; par_new']);
            elseif samp<=nburnin && n_accept>1
                pbar_prev=pbar;
                pbar=(n_accept./(n_accept+1)).*pbar_prev+...
                     (1./(n_accept+1)).*par_new;
                SIG=((n_accept-1)./n_accept).*SIG+...
                    (1./n_accept).*(n_accept.*(pbar_prev*pbar_prev')-...
                                    (n_accept+1).*(pbar*pbar')+...
                                    (par_new*par_new'));
            end
 
% Update the chain
            CurrL=NewL;
            prior=prior_new;
            par=par_new;
 
        end
 
% Every one hundred samples during burn-in, tune the scaling factor
% for the proposal distribution to ensure an acceptance rate of 20-40%
        if samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp<0.2
            sf=sf/2;
        elseif samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp>0.4
            sf=2*sf;
        end
%==========================================================================
 
%==========================================================================
% STORE THE OUTPUT
% After burn in, save iterations of the chain, thinning as specified
        if nthin==1
            ParSampC(samp,:)=[par' prior CurrL];
        elseif samp>nburnin && mod(samp,nthin)==1
            ParSampC(iter,:)=[par' prior CurrL];
            iter=iter+1;
        end
%==========================================================================
 
    end
 
% Store the chain
    ParSamp{chain}=ParSampC;
 
end
 
%==========================================================================
% COMPUTE DIC AND pD
% Compute the deviance for each sample
Dev=[];
PS=[];
for chain=1:nchains
    Dev=[Dev; -2*ParSamp{chain}(:,end)];
    PS=[PS; ParSamp{chain}(:,1:end-2)];
end
 
% Compute the mean deviance
Dbar=mean(Dev);
 
% Compute the deviance at the posterior mean for the parameters
Dhat=-2*Lhood_FMDV(mean(PS,1)',tClin,tStart,tStop,D,tObs,VI,mFlag);
 
% Compute the DIC
DIC=2*Dbar-Dhat;
 
% Compute the effective number of parameters
pD=Dbar-Dhat;
%==========================================================================
 
% Save the outputs
save(['FMDVModel' num2str(mFlag) '_' dtype '_MCMCSamples'],...
     'ParSamp','nburnin','nsamp','DIC','pD','seeds')
