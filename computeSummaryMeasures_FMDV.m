% Matlab script to compute summary tranmsission measures (R0 and generation
% time) and theta for FMDV using each proxy measure for infectiousness

% Specify the model used: transmission proportional to titre (1,3) or log
% titre (2,4); and transmission rate common to (1,2) or varies amongst (3,4)
% animals
mFlag=4;

% Specify the measures used to infer infectiousness (Blood, NF or
% OPF)
meas={'Blood', 'NF', 'OPF'};

% Create cell arrays to store the measures (explicit and approximate)
R0=cell(1,length(meas));
Tg=cell(1,length(meas));
R0A=cell(1,length(meas));
TgA=cell(1,length(meas));
theta=cell(1,length(meas));

% Load the times of clinical onset for each animal
varload=load('FMDV_TransmissionExperimentData');
tClin=varload.tClin(1:7);

% Set the number of animals
nanim=7;

% For each measure ...
for m=1:length(meas)

% Load the MCMC output
    varload=load(['FMDVModel' num2str(mFlag) '_' meas{m} '_MCMCSamples']);
    PS=[varload.ParSamp{1}(:,1:end-2); varload.ParSamp{2}(:,1:end-2)];

% Extract the hierarchical shedding parameters
    mu_tp=PS(:,4*nanim+3);
    sig_tp=PS(:,4*nanim+4);

% Extract the incubation period parameters
    muC=PS(:,4*nanim+9);
    sigC=PS(:,4*nanim+10);
    rhoTC=PS(:,4*nanim+11);

% Create the array to store the summary transmission measures
    R0{m}=zeros(size(PS,1),nanim);
    Tg{m}=zeros(size(PS,1),nanim);
    R0A{m}=zeros(size(PS,1),nanim);
    TgA{m}=zeros(size(PS,1),nanim);
    theta{m}=zeros(size(PS,1),nanim);
       
% For each animal ...
    for j=1:nanim
        disp([meas{m} '; animal: ' num2str(j)])
        disp('  calculating R0 and Tg')

%==========================================================================
% PREPARE THE INDIVIDUAL-LEVEL PARAMETERS
% Extract the individual shedding parameters
        Vp=exp(PS(:,j));
        tp=PS(:,nanim+j);
        lg=PS(:,2*nanim+j);
        ld=PS(:,3*nanim+j);

% Extract the transmission parameter
        if mFlag==1 || mFlag==2
            g=exp(PS(:,end-1));
        elseif mFlag==3
            g=exp(PS(:,4*nanim+11+j));
        elseif mFlag==4
            g=PS(:,4*nanim+11+j);
        end
%==========================================================================
        
%==========================================================================
% COMPUTE VIRAL AND SHEDDING CURVES
% Set the times post infection to compute infectiousness (ensuring all the
% complete curves are included)
        t02=tp+(1./ld).*log(2*Vp);
        tau=0:0.01:min(100,max(ceil(t02)));

% Compute the viral curves ...
        VpR=repmat(Vp,1,length(tau));
        tpR=repmat(tp,1,length(tau));
        lgR=repmat(lg,1,length(tau));
        ldR=repmat(ld,1,length(tau));
        tauR=repmat(tau,length(Vp),1);
        V=2.*VpR./(exp(-lgR.*(tauR-tpR))+exp(ldR.*(tauR-tpR)));

% ... and the corresponding shedding
        if mFlag==1 || mFlag==3
            S=V;
        elseif mFlag==2 || mFlag==4
            S=log(V);
            S(S<0)=0;
        end

% Compute the area under the shedding curve
        intS=trapz(tau,S,2);
%==========================================================================
        
%==========================================================================
% COMPUTE R0 AND GENERATION TIME
% Compute R0 and Tg
        R0{m}(:,j)=g.*intS;
        Tg{m}(:,j)=trapz(tau,repmat(tau,length(g),1).*S,2)./intS;

% Compute the approximations for R0 and Tg
    	R0A{m}(:,j)=0.5.*g.*(1./lg+1./ld).*log(Vp).*(log(Vp)+log(4));
        TgA{m}(:,j)=tp+(1./3)*(1./ld-1./lg).*(log(Vp)+log(2));
%==========================================================================

%==========================================================================
% COMPUTE THETA
% Indicate what it's doing
        disp('  calculating theta')
        
% Compute the CDF for the time of onset of clinical signs, conditional on
% (i) the time of peak titre and (ii) the interval when onset occurred
        mu=muC+rhoTC.*(sigC./sig_tp).*(log(tp)-mu_tp);
        sig=sigC.*sqrt(1-rhoTC.^2);
        F=max(0,min(1,(logncdf(repmat(tau,length(mu),1),...
                               repmat(mu,1,length(tau)),...
                               repmat(sig,1,length(tau)))-...
                       repmat(logncdf(tClin(j)-1,mu,sig),1,length(tau)))./...
                      repmat(logncdf(tClin(j),mu,sig)-logncdf(tClin(j)-1,mu,sig),...
                             1,length(tau))));

% Compute the proportion of infectiousness prior to onset
        theta{m}(:,j)=trapz(tau,S.*(1-F),2)./intS;
%==========================================================================

    end
end

% Save the results
save(['FMDVModel' num2str(mFlag) '_SummaryTransmissionMeasures'],...
     'R0','Tg','R0A','TgA','theta','meas')

% Tidy up
close('all')
clear
