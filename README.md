# TRIO_RVEMVS
The codes in this folder is used for implementing method TRIO_RVEMVS which is proposed from
paper "TRIO RVEMVS: A Bayesian framework for rare variant association analysis with expectation-maximization variable selection using family trio data".
TRIO_RVEMVS is an expectation maximization variable selection (EMVS) method to simultaneously detect common and rare variants at the individual variant level using family trio data.

# code description
rv_emvs.prepare.R: this file is used to pre-process the trio data set for rv_emvs; 
                   
                   rv_emvs.prepare('example.trios.csv')
                
rv_emvs.R: the main function of TRIO_RVEM
           
           rv_emvs(prep.list[[1]], maf = prep.list[[2]] , map = map.example , common.inclusion = 0.1, rare.inclusion = 0.5,
                                 common.exclusion = (log(sqrt(1.015/0.985))/qnorm(.975))^2,
                                 rare.exclusion = (log(sqrt(1.01/0.99))/qnorm(.975))^2,
                                 threshold = 0.05, theta0 = 0.5,        # the common  inclusion probability
                                 gth0 = 0.5,
                                 pi0 = 0.5                           # individual rare variants inclusion probability
                                 ) 

rv_emvs.regularization.R: code for tunning exclusion parameter tuning using regularization plot, which generates "example regularization plot.pdf".

# reference
Duo Yu, et al. "TRIO RVEMVS: A Bayesian framework for rare variant association analysis with expectation-maximization variable selection using family trio data".

# contact information
duoyu@mcw.edu
