#########################################################
# Bayesian Forward Continuation Ratio Model Example Code 
# Anna Eames Seffernick
# August 13, 2023
########################################################

# Note you must have JAGS installed. See https://mcmc-jags.sourceforge.io

# Load required packages
library(runjags)
library(VGAM)
library(parallel)
library(tidyr)
library(GEOquery)

seed1.tmp <- 22204
seed2.tmp <- 21229

# Load the data, called sm.sim
load("BLASSO-FCR_SmallData.RData")
# This data has 20 features, 
# the first 10 of which are truly associated with the outcome

################
# Fit the Model
################

## Restructure the data

temp.data.m <- sm.sim
Xmat <- as.matrix(temp.data.m[,-c(1,2,3)])
temp.data <- as.data.frame(temp.data.m)
y <- temp.data$OS5
relapsed <- temp.data$censor

levels <- sort(unique(y))
K <- length(unique(y))
Ymat <- matrix(0, nrow=length(y), ncol=K)

c.y <- ifelse(relapsed == 1, 0, 1) * y
# y.event = time point subj experienced the event; =0 if subj was censored
y.event <- y * relapsed
Tmat <- matrix(0, nrow=length(y), ncol=K)
for(j in levels){
  Ymat[which(y.event == j), which(levels == j)] <- 1
  Tmat[which(y.event == j), which(levels == j)] <- 1
  if(j!=K) {
    Tmat[which(c.y == (j + 1)), which(levels == j)] <- 1
  }
}
YColsum <- colSums(Ymat)

alpha.vec <- numeric()
pi.0 <- table(y)/length(y)
tab <- table(y)

Cum.Tmat <- matrix(0, nrow=nrow(Tmat), ncol=K)
for(i in 1:(K-1)) {
  alpha.vec[i] <- log(-log(1 - (tab[i] / sum(tab[i:K]))))
  Cum.Tmat[, i] <- rowSums(Tmat[, i:K])
}
Cum.Tmat[, K] <- Tmat[, K]
alpha.vec[K] <- log(-log(1 - .99))

data1_wide <- cbind.data.frame(Ymat, Xmat)
colnames(data1_wide)[1:5] <- c("Y1", "Y2", "Y3", "Y4", "Y5")
dontchange <- colnames(data1_wide)[6:length(colnames(data1_wide))]
data1_long <- gather(data1_wide, key="Y", value="Succ", -all_of(dontchange))
data1_wide_T <- cbind.data.frame(Cum.Tmat, Xmat)
colnames(data1_wide_T)[1:5] <- c("CumT1", "CumT2", "CumT3", "CumT4", "CumT5")
dontchange2 <- colnames(data1_wide_T)[6:length(colnames(data1_wide_T))]
data1_long_T <- gather(data1_wide_T, key="CumT", value="Trials", -all_of(dontchange2))
data1_long_comb <- cbind.data.frame(data1_long, data1_long_T[,grep("Trials", colnames(data1_long_T))])
data_long_cens1 <- data1_long_comb
colnames(data_long_cens1)[length(colnames(data_long_cens1))] <- "Trials"
data_long_cens1_sm <- data_long_cens1[-which(data_long_cens1$Trials==0),]
data_long_cens1_sm$cp1 <- ifelse(data_long_cens1_sm$Y=="Y1", 1, 0)
data_long_cens1_sm$cp2 <- ifelse(data_long_cens1_sm$Y=="Y2", 1, 0)
data_long_cens1_sm$cp3 <- ifelse(data_long_cens1_sm$Y=="Y3", 1, 0)
data_long_cens1_sm$cp4 <- ifelse(data_long_cens1_sm$Y=="Y4", 1, 0)
data_long_cens1_sm$cp5 <- ifelse(data_long_cens1_sm$Y=="Y5", 1, 0)


# Write the model file
Model = "model{
  for(i in 1:N){
    mu[i] <- inprod(betgma[], X[i,])
    pi[i] <- icloglog(mu[i])
    Y[i] ~ dbern(pi[i])
  }
  for(l in 1:K){
    beta0[l] ~ dnorm(0, 0.1)
  }
  beta[1:K] <- sort(beta0[1:K])
  for(k in (K+1):P){
     beta[k] ~ ddexp(0, lambda)
  }
  lambda ~ dgamma(0.1, 0.1)
  for(s in 1:K){
    gamma[s] <- 1
    betgma[s] <- beta[s]*gamma[s]
  }
  for(m in (K+1):P){
    betgma[m] <- beta[m]*gamma[m]
    gamma[m] ~ dbern(0.1)
  }
  
  #inits# beta, lambda, .RNG.seed, .RNG.name
  #monitor# beta, lambda, gamma, betgma
  #modules# glm on
}"


JAGSFILE="BLassoDisSurv.bug"
cat(Model, file=JAGSFILE)

##Set Model Parameters
Xmat <- data_long_cens1_sm[,c(grep("cp", colnames(data_long_cens1_sm)), grep("at", colnames(data_long_cens1_sm)))]
N <- dim(data_long_cens1_sm)[1]
P <- dim(Xmat)[2]

#Data List
dataList <- list("Y" = data_long_cens1_sm$Succ, "X" = as.matrix(Xmat),
                 "N" = N, "P" = P, "K"=K)
# Parameters to be monitored
parameters <- c("beta", "gamma", "lambda", "betgma")

# what are we supposed to guess for alpha3?
seed1 <- seed1.tmp + 1
set.seed(seed1)
inits1 <- list("beta0" = alpha.vec,"beta" = c(NA, NA, NA, NA, NA, rep(0.0, P-5)), "lambda" = rgamma(1, shape=0.1, rate=0.1))
inits2 <- list("beta0" = alpha.vec+0.2, "beta" = c(NA, NA,NA, NA, NA, rep(0.0, P-5)), "lambda" = rgamma(1, shape=0.1, rate=0.1))
inits3 <- list("beta0" = alpha.vec+0.05, "beta" = c(NA, NA,NA, NA, NA, rep(0.0, P-5)), "lambda" = rgamma(1, shape=0.1, rate=0.1))
inits.list <- list(inits1, inits2, inits3)
names(inits.list) <- c("chain1", "chain2", "chain3")

.RNG.seed <- function(chain){
  return(switch(chain, "1" = seed1+1, "2" = seed1+2, "3"=seed1+3))
}

.RNG.name <- function(chain){
  return(switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill", "3"="base::Super-Duper"))
}

library(parallel)
cl <- makeCluster(10)
seed2 <- seed2.tmp + 1
set.seed(seed2)
model.fit.parallel.post.BLI <- run.jags(model=Model, data=dataList, n.chains=3, inits=inits.list, burnin=500, adapt=500, sample=3333, thin=3, method="parallel", cl=cl)
stopCluster(cl)


############################
# Explore Posterior Samples
############################

# load libraries
library(rjags)
library(dclone)
library(ordinalgmifs)
library(glmnetcr)
library(Biobase)
library(BiocGenerics)
library(affy)


# Define Functions to summarize posterior distribution
# Bayes Factor functions
prior.odds.beta.gamma.func <- function(a, b, epsilon, pgamma1){
  pgamma0 <- 1 - pgamma1
  prior.odds <- (pgamma1*(b^a)*gamma(a+1))/(pgamma1*((b+epsilon)^a - b^a)*gamma(a+1) + pgamma0*a*((b+epsilon)^a)*gamma(a))
  prior.odds
}
prior.odds.beta.func <- function(a, b, epsilon){
  prior.odds <- (b^a)/((b+epsilon)^(a) - b^a)
  prior.odds
}
prior.odds.gamma.func <- function(pgamma1){
  pgamma1/(1-pgamma1)
}
post.odds.beta.gamma.func <- function(mcmcChainName, epsilon){
  post.odds <- c()
  for(l in 1:dim(mcmcChainName)[2]){
    param.test <- mcmcChainName[,l]
    P.Gr.Eps <- sum(ifelse(abs(param.test) > epsilon, 1, 0))/dim(mcmcChainName)[1]
    P.LE.Eps <- sum(ifelse(abs(param.test) <= epsilon, 1, 0))/dim(mcmcChainName)[1]
    post.odds[l] <- P.Gr.Eps/P.LE.Eps
  }
  post.odds
}
post.odds.gamma.func <- function(mcmcChainName){
  post.odds <- c()
  for(j in 1:dim(mcmcChainName)[2]){
    param.test <- mcmcChainName[,j]
    P.Gamma.1 <- mean(param.test)
    P.Gamma.0 <- 1 - mean(param.test)
    post.odds[j] <- P.Gamma.1/P.Gamma.0
  }
  post.odds
}
bayes.factor.func <- function(post.odds, prior.odds){
  BF <- c()
  for(k in 1:length(post.odds)){
    BF[k] <- post.odds[k]/prior.odds
  }
  BF
}

# FDR Function
fdr.func <- function(num.dis, num.false.dis){
  ifelse(num.dis==0, NA, num.false.dis/num.dis)
}
# TPR Function
tpr.func <- function(num.pos, num.true.dis){
  num.true.dis/num.pos
}
# TNR Function
tnr.func <- function(num.neg.dis.cor, num.neg){
  num.neg.dis.cor/num.neg
}
# PPV Function
ppv.func <- function(num.true.dis, num.dis){
  ifelse(num.dis==0, NA, num.true.dis/num.dis)
}
# NPV Function
npv.func <- function(num.neg.dis.cor, num.neg.dis){
  ifelse(num.neg.dis==0, NA, num.neg.dis.cor/num.neg.dis)
}

# Look at results
sig.feats <- 1:10
not.sig.feats <- 11:20
mcmcList <- as.mcmc.list(model.fit.parallel.post.BLI)
results.m <- as.matrix(mcmcList)
mcmcAlpha <- results.m[,1:5]
mcmcBeta <- results.m[,6:25]
mcmcLambda <- results.m[,26]
mcmcGamma <- results.m[,32:51]
mcmcBetaGamma <- results.m[,57:76]
sum.fit <- summary(model.fit.parallel.post.BLI)
sum.fit
# Extract needed components - HPDI
HPDI.dis <- which(sign(sum.fit[6:25,grep("Lower95", colnames(sum.fit))]) == sign(sum.fit[6:25,grep("Upper95", colnames(sum.fit))]))
HPDI.neg <- which(sign(sum.fit[6:25,grep("Lower95", colnames(sum.fit))]) != sign(sum.fit[6:25,grep("Upper95", colnames(sum.fit))]))
# Extract needed components - Bayes Factor
# Bayes Factor Method
a <- 0.1
b <- 0.1
epsilon <- 0.1
pgamma1 <- 0.1
bg.prior <- prior.odds.beta.gamma.func(a, b, epsilon, pgamma1)
beta.prior <- prior.odds.beta.func(a, b, epsilon)
gamma.prior <- prior.odds.gamma.func(pgamma1)
bg.post <- post.odds.beta.gamma.func(mcmcBetaGamma, epsilon)
beta.post <- post.odds.beta.gamma.func(mcmcBeta, epsilon)
gamma.post <- post.odds.gamma.func(mcmcGamma)
BF.beta <- bayes.factor.func(beta.post, beta.prior)
BF.betagma <- bayes.factor.func(bg.post, bg.prior)
BF.gamma <- bayes.factor.func(gamma.post, gamma.prior)
BF.bg.dis <- which(BF.betagma > 5)
BF.bg.neg <- which(BF.betagma < 5)
BF.beta.dis <- which(BF.beta > 5)
BF.beta.neg <- which(BF.beta < 5)
BF.gamma.dis <- which(BF.gamma > 5)
BF.gamma.neg <- which(BF.gamma < 5)
post.dis <- which(colMeans(mcmcGamma) >= 0.5)
post.neg <- which(colMeans(mcmcGamma) < 0.5)
# Save Variable Selection Results
v1 <- c("Method", "Discoveries", "FDR", "TPR", "TNR", "PPV", "NPV")
# Save HPDI variable selection results
v2 <- c("HDPI", length(HPDI.dis),  fdr.func(length(HPDI.dis), length(which(HPDI.dis %in% not.sig.feats))), tpr.func(length(sig.feats), length(which(HPDI.dis %in% sig.feats))), tnr.func(length(HPDI.neg) - length(which(HPDI.neg %in% sig.feats)), length(HPDI.dis) + length(HPDI.neg)-length(sig.feats)), ppv.func(length(which(HPDI.dis %in% sig.feats)), length(HPDI.dis)), npv.func(length(HPDI.neg) - length(which(HPDI.neg %in% sig.feats)), length(HPDI.neg)))
# Save BF for betavariable selection results
v3 <- c("beta*gamma BF > 5", length(BF.bg.dis), fdr.func(length(BF.bg.dis), length(which(BF.bg.dis %in% not.sig.feats))), tpr.func(length(sig.feats), length(which(BF.bg.dis %in% sig.feats))), tnr.func(length(BF.bg.neg) - length(which(BF.bg.neg %in% sig.feats)), length(BF.bg.dis) + length(BF.bg.neg)-length(sig.feats)), ppv.func(length(which(BF.bg.dis %in% sig.feats)), length(BF.bg.dis)), npv.func(length(BF.bg.neg) - length(which(BF.bg.neg %in% sig.feats)), length(BF.bg.neg)))
v4 <- c("beta BF > 5", length(BF.beta.dis), fdr.func(length(BF.beta.dis), length(which(BF.beta.dis %in% not.sig.feats))), tpr.func(length(sig.feats), length(which(BF.beta.dis %in% sig.feats))), tnr.func(length(BF.beta.neg) - length(which(BF.beta.neg %in% sig.feats)), length(BF.beta.dis) + length(BF.beta.neg)-length(sig.feats)), ppv.func(length(which(BF.beta.dis %in% sig.feats)), length(BF.beta.dis)), npv.func(length(BF.beta.neg) - length(which(BF.beta.neg %in% sig.feats)), length(BF.beta.neg)))
v5<- c("gamma BF > 5", length(BF.gamma.dis), fdr.func(length(BF.gamma.dis), length(which(BF.gamma.dis %in% not.sig.feats))), tpr.func(length(sig.feats), length(which(BF.gamma.dis %in% sig.feats))), tnr.func(length(BF.gamma.neg) - length(which(BF.gamma.neg %in% sig.feats)), length(BF.gamma.dis) + length(BF.gamma.neg)-length(sig.feats)), ppv.func(length(which(BF.gamma.dis %in% sig.feats)), length(BF.gamma.dis)), npv.func(length(BF.gamma.neg) - length(which(BF.gamma.neg %in% sig.feats)), length(BF.gamma.neg)))
v6 <- c("Post Mean", length(post.dis),  fdr.func(length(post.dis), length(which(post.dis %in% not.sig.feats))), tpr.func(length(sig.feats), length(which(post.dis %in% sig.feats))), tnr.func(length(post.neg) - length(which(post.neg %in% sig.feats)), length(post.dis) + length(post.neg)-length(sig.feats)), ppv.func(length(which(post.dis %in% sig.feats)), length(post.dis)), npv.func(length(post.neg) - length(which(post.neg %in% sig.feats)), length(post.neg)))

# Create and view data frame of results
VI.temp.df <- rbind.data.frame(v2, v3, v4, v5, v6)
colnames(VI.temp.df) <- v1
VI.temp.df
#             Method Discoveries FDR TPR TNR PPV NPV
#               HDPI          10   0   1   1   1   1
#  beta*gamma BF > 5          10   0   1   1   1   1
#        beta BF > 5          10   0   1   1   1   1
#       gamma BF > 5          10   0   1   1   1   1
#          Post Mean          10   0   1   1   1   1
# All 10 significant features were discovered 


#######################
# Session Info
#######################

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6.3
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  splines   stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] affy_1.72.0               glmnetcr_1.0.6            glmnet_4.1-4              ordinalgmifs_1.0.7       
# [5] dclone_2.3-0              Matrix_1.4-1              rjags_4-13                coda_0.19-4              
# [9] tidyr_1.2.0               VGAM_1.1-7                runjags_2.2.1-7           flexsurv_2.2             
# [13] GEOquery_2.62.2           riskRegression_2023.03.22 MatrixModels_0.5-0        coin_1.4-2               
# [17] readxl_1.4.0              data.table_1.14.2         knitr_1.39                penalized_0.9-52         
# [21] DescTools_0.99.49         robustbase_0.99-0         Rfit_0.24.2               cmprsk_2.2-11            
# [25] survival_3.3-1            lawstat_3.6               Biobase_2.54.0            BiocGenerics_0.40.0      
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.4.1       Hmisc_5.1-0           mstate_0.3.2          listenv_0.8.0         ggplot2_3.3.6        
# [6] TH.data_1.1-1         inline_0.3.19         digest_0.6.29         foreach_1.5.2         htmltools_0.5.2      
# [11] fansi_1.0.3           magrittr_2.0.3        checkmate_2.1.0       cluster_2.1.3         tzdb_0.3.0           
# [16] limma_3.50.3          globals_0.15.1        readr_2.1.2           mets_1.3.2            RcppParallel_5.1.5   
# [21] matrixStats_0.62.0    Kendall_2.2.1         sandwich_3.0-2        prettyunits_1.1.1     colorspace_2.0-3     
# [26] rbibutils_2.2.8       xfun_0.31             dplyr_1.0.9           callr_3.7.0           crayon_1.5.1         
# [31] libcoin_1.0-9         Exact_3.2             zoo_1.8-10            iterators_1.0.14      glue_1.6.2           
# [36] gtable_0.3.0          zlibbioc_1.40.0       pkgbuild_1.3.1        rstan_2.21.5          shape_1.4.6          
# [41] future.apply_1.9.0    rms_6.7-0             DEoptimR_1.1-0        SparseM_1.81          scales_1.2.0         
# [46] mvtnorm_1.1-3         DBI_1.1.3             Rcpp_1.0.9            htmlTable_2.4.1       foreign_0.8-82       
# [51] proxy_0.4-27          preprocessCore_1.56.0 deSolve_1.32          Formula_1.2-4         StanHeaders_2.21.0-7 
# [56] lava_1.7.2.1          prodlim_2019.11.13    htmlwidgets_1.5.4     httr_1.4.3            modeltools_0.2-23    
# [61] ellipsis_0.3.2        loo_2.5.1             pkgconfig_2.0.3       nnet_7.3-17           utf8_1.2.2           
# [66] tidyselect_1.1.2      rlang_1.0.3           munsell_0.5.0         cellranger_1.1.0      tools_4.1.2          
# [71] cli_3.3.0             generics_0.1.3        evaluate_0.15         stringr_1.4.0         fastmap_1.1.0        
# [76] yaml_2.3.5            processx_3.7.0        timereg_2.0.5         purrr_0.3.4           rootSolve_1.8.2.3    
# [81] future_1.26.1         nlme_3.1-158          quantreg_5.93         xml2_1.3.3            compiler_4.1.2       
# [86] rstudioapi_0.13       affyio_1.64.0         e1071_1.7-11          tibble_3.1.7          statmod_1.4.36       
# [91] stringi_1.7.6         ps_1.7.1              lattice_0.20-45       vctrs_0.4.1           pillar_1.7.0         
# [96] lifecycle_1.0.1       BiocManager_1.30.18   Rdpack_2.3.1          lmom_2.9              R6_2.5.1             
# [101] muhaz_1.2.6.4         gridExtra_2.3         parallelly_1.32.0     gld_2.6.6             codetools_0.2-18     
# [106] polspline_1.1.22      boot_1.3-28           MASS_7.3-57           assertthat_0.2.1      withr_2.5.0          
# [111] multcomp_1.4-19       expm_0.999-7          hms_1.1.1             quadprog_1.5-8        grid_4.1.2           
# [116] rpart_4.1.16          class_7.3-20          rmarkdown_2.14        numDeriv_2016.8-1.1   base64enc_0.1-3      
