% Matlab script to compute summary transmission measures (R0 and
% generation time) for swine influenza virus

% Specify the model used: transmission proportional to titre (1,3) or log
% titre (2,4); and transmission rate common to (1,2) or varies amongst (3,4)
% animals
mFlag=4;

% Set the number of animals
nanim=11;

% Load the MCMC output
varload=load(['SwIVModel' num2str(mFlag) '_MCMCSamples']);
PS=[varload.ParSamp{1}(:,1:end-2); varload.ParSamp{2}(:,1:end-2)];

% Create the arrays to store the measures (explicit and approximate)
R0=zeros(size(PS,1),nanim);
Tg=zeros(size(PS,1),nanim);
R0A=zeros(size(PS,1),nanim);
TgA=zeros(size(PS,1),nanim);
       
% For each animal ...
for j=1:nanim
    disp(['animal: ' num2str(j)])

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
        g=exp(PS(:,4*nanim+8+j));
    elseif mFlag==4
        g=PS(:,4*nanim+8+j);
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
% Compute R0 and Tg based on the explicit formulae
    R0(:,j)=g.*intS;
    Tg(:,j)=trapz(tau,repmat(tau,length(g),1).*S,2)./intS;

% Compute the approximations for R0 and Tg
    R0A(:,j)=0.5.*g.*(1./lg+1./ld).*log(Vp).*(log(Vp)+log(4));
    TgA(:,j)=tp+(1./3)*(1./ld-1./lg).*(log(Vp)+log(2));
%==========================================================================

end

% Save the results
save(['..\SwIVModel' num2str(mFlag) '_SummaryTransmissionMeasures'],...
     'R0','Tg','R0A','TgA')

% Tidy up
close('all')
clear
