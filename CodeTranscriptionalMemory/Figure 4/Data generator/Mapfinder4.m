ClonekeyN = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Jam2_S1_STRSa1','Jarid2_S1_STRSa1','Rbpj_S1_STRSa1','Nono_S1_STRSa1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1','Tfpd2_S1_STRSa1'};

zN = 4%:length(Clonekey)
load(char(ClonekeyN(zN)))
%%
%%
Thresh = 10000;
X = mean(X(Thresh:5:end,1:7));

X([2,4]) = [];

X(3) = X(3)^2;
X(4) = X(4)^2;
X(1:4) = log(X(1:4));

Xtemp = X;
X(1) = Xtemp(2);
X(2) = Xtemp(5);
Xrep = [Xtemp(1);Xtemp(3);Xtemp(4)];
Xallvars = repmat(Xrep,round(size(yTOT,2)),1);
X = [Xtemp(2);0.7;0;0;0;0;0;0;0;Xallvars];
R0R0coeff = 0;

% [yTOT] = Randomise_yTOT(yTOT);


itlimit1 = 120;

logpdf = @(X) logPosteriorReworked_noprior(X,x,yTOT,gamma,R0R0coeff);
startpoint = X;
smp = hmcSampler(logpdf,startpoint,'NumSteps',100,'CheckGradient',false,'StepSize',0.001,'JitterMethod','jitter-both');
[MAPpars,fitInfo] = estimateMAP(smp,'VerbosityLevel',1,'IterationLimit',itlimit1);

%%

meanextrct = MAPpars(10:3:end);
sigmaD1extrct = MAPpars(11:3:end);
sigmaD2extrct = MAPpars(12:3:end);
normedsigmaD1 = sqrt(exp(sigmaD1extrct));
normedsigmaD2 = sqrt(exp(sigmaD2extrct));
mu_mu_guess = mean(meanextrct);
mu_D1_guess = mean(log(normedsigmaD1));
mu_D2_guess = mean(log(normedsigmaD2));
sigma_mu_guess = std(meanextrct,1); % use biased estimator - dividing by N not N-1
sigma_D1_guess = std(log(normedsigmaD1),1);
sigma_D2_guess = std(log(normedsigmaD2),1);
corrguess = corr(meanextrct(1:2:end),meanextrct(2:2:end));

X = MAPpars;
X(3) = mu_mu_guess;
X(4) = mu_D1_guess;
X(5) = mu_D2_guess;
X(6) = log(sigma_mu_guess);
X(7) = log(sigma_D1_guess);
X(8) = log(sigma_D2_guess);
X(9) = atanh(corrguess);

%%


logpdf = @(X) logPosteriorReworked(X,x,yTOT,gamma,R0R0coeff);
startpoint = X;
smp = hmcSampler(logpdf,startpoint,'NumSteps',100,'CheckGradient',false,'StepSize',0.01,'JitterMethod','jitter-both');

NumChains = 4;
chains = cell(NumChains,1);
Burnin = 200;
NumSamples = 2500;%150;
for c = 1:NumChains
    if (c == 1)
        level = 1;
    else
        level = 0;
    end
    chains{c} = drawSamples(smp,'Start',X, ...
        'Burnin',Burnin,'NumSamples',NumSamples, ...
        'VerbosityLevel',level,'NumPrint',1);
end


save([char(ClonekeyN(zN)),'_HMCprelim'])