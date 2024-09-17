# main function of rv_emvs
# Duo.Yu@uth.tmc.edu

require(foreach)

# seperate common and rare variants given trio data set
common_seperate = function(prep.trio, maf, threshold){
  # prep_trio_file: prepared trio file
  # maf_file: the maf information for each SNP accross 500 data set
  # return: the common varints names
  trio_maf = maf[which(names(maf) %in% colnames(prep.trio))]
  names(trio_maf[which(trio_maf >= threshold)])
}

# sum by strata, repeated each the number of controls if necessary
sumbys = function(x,strata,r=TRUE){
  ust = unique(strata)
  z = matrix(NA,nrow = length(ust),ncol = 1)
  cont = 0
  for (i in ust )
  {
    cont = cont + 1
    z[cont,1] = sum(x[which(strata == i)])
  }
  if (r)
  {z = rep(z,each = (length(x)/length(ust)))
  }
  return(z)
}

# inclusion paramter p_star for common variants
inclufun_common = function(beta0,theta0,t, common.inclusion, common.exclusion)
{
  p_star = (dnorm(as.numeric(beta0), mean = 0,sd = sqrt(common.inclusion))*theta0)^t /
    ((dnorm(as.numeric(beta0), mean = 0, sd = sqrt(common.inclusion))*theta0)^t +
        (dnorm(as.numeric(beta0), mean = 0, sd = sqrt(common.exclusion))*(1 - theta0))^t)
  P = (1 - p_star)/common.exclusion + p_star/common.inclusion
  return(list(p_star = p_star,P = P))
}

# inclusion paramter g_star for group
inclufun_group = function(alpha0, gth0, pi0, region,t, rare.map, rare.inclusion, rare.exclusion)
{
  almark =  alpha0
  g_star =  c()
  for (i in 1:region)
  {
    alpha0 = almark[which(rare.map$gene == i)]
    cup = exp(log((1 - gth0)/gth0) - sum(log(1 - pi0*(1 - dnorm(as.numeric(alpha0), mean = 0,
      sd = sqrt(rare.inclusion))/dnorm(as.numeric(alpha0), mean = 0,sd = sqrt(rare.exclusion))))))
    g_star = rbind(g_star,1/(1 + cup^t))
  }
  return(g_star)
}

# inclusion parameter q_star for individual variants
inclufun_rare = function(alpha0, pi0,t, rare.inclusion, rare.exclusion)
{
  q_star = (dnorm(as.numeric(alpha0), mean = 0,sd = sqrt(rare.inclusion))*pi0)^t /
    ((dnorm(as.numeric(alpha0), mean = 0,sd = sqrt(rare.inclusion))*pi0)^t + (dnorm(as.numeric(alpha0),
      mean = 0,sd = sqrt(rare.exclusion))*(1 - pi0))^t)
  Q_rare = (1 - q_star)/rare.exclusion + q_star/rare.inclusion

  return(list(q_star = q_star, Qr = Q_rare))
}

Qlogistic = function(data,beta_test,alpha_test,theta,gth,pii,p_star,g_star,q_star,
  strata, common.exclusion, common.inclusion, rare.exclusion, rare.inclusion,rare.map, b1, b2, b3, Z)
{
  beta = c(beta_test,alpha_test)

  Q = -sum(log(1 + sumbys(exp(-data %*% beta),strata, r = 0))) -
    1/2*t(beta_test) %*% (beta_test*((1 - p_star)/common.exclusion + p_star/common.inclusion)) -
    # update error common.exclusion and common.inclusion to rare.exclusion and rare.inclusion. date: 2.9. 2019
    1/2*t(alpha_test) %*% (alpha_test*((1 - q_star)/rare.exclusion + q_star/rare.inclusion)) +
    (sum(p_star))*log(theta) + (b1 - 1 + length(p_star) - sum(p_star))*log(1 - theta) +
    (b2 + length(g_star) - sum(g_star) - 1)*log(1 - gth) + sum(g_star)*log(gth) +
    # t(Z) %*% q_star only sum the inclusion probability inside each gene, return vector of 12 elem
    sum(g_star*(t(Z) %*% q_star))*log(pii) + (sum(as.vector(g_star)*table(rare.map$gene)) +
        b3 - 1 - sum(g_star*(t(Z) %*% q_star)))*log(1 - pii)
  return(Q)
}

Rlogistic = function(p_star,g_star, q_star,theta, gth, pii, rare.map, b1, b2, b3, Z)
{
  R = (sum(p_star))*log(theta) + (b1 - 1 + length(p_star) -
      sum(p_star))*log(1 - theta) + (b2 + length(g_star) - sum(g_star) - 1)*log(1 - gth) +
    sum(g_star)*log(gth) + sum(g_star*(t(Z) %*% q_star))*log(pii) +
    (sum(as.vector(g_star)*table(rare.map$gene)) + b3 - 1 - sum(g_star*(t(Z) %*% q_star)))*log(1-pii)
  return(R)
}

# gradient by family
grad.diff3M = function(x,beta,strata,sub_fam)
{
  U = exp(-x[strata == sub_fam,] %*% beta)/(1 + sumbys(exp(-x[strata == sub_fam,] %*% beta),rep(1,3),
    r = 3))
  res = as.vector(-t(x[strata == sub_fam,]) %*% U)
  return(res)}


#### SDCA in mini-batch setting

SDCA1 = function(num.var, num.sample, x.mat, P,seeds, strata) {


  m = 2 ## mini batch size
  Nbatch = num.sample/m

  alpha1 = matrix(0,ncol(x.mat),num.sample)
  alpha0 = matrix(0,ncol(x.mat),num.sample)
  bar_alpha = matrix(0,ncol(x.mat),1)
  # v.alpha = c()
  # cost = function(beta,P,x.mat, strata){
  #  1/num.sample*sum(1 + sumbys(exp(-x.mat %*% beta),strata,r = 0)) + 1/2/num.sample*sum((beta)^2*P)
  #}
  ################ INITIAL VALUES ############################

  beta = matrix(0,ncol(x.mat),1)
  #vcost = cost(beta,P,x.mat)
  s = 0.01  ## the parameter  used in  original paper
  u = c()
  epoch = 1
  while (epoch < 4)
  {
    tot.it = 1
    whole.col = c(1:num.sample)


    while (tot.it < (Nbatch + 1)) {

      u = (1 - s) * beta + s * num.sample*bar_alpha/P
      set.seed(seeds[epoch])
      random.seq = sample(whole.col, m,replace = FALSE)   # ramdom select a mini batch
      whole.col = whole.col[-which(whole.col %in% random.seq)]
      #used.col = cbind(used.col,random.seq)
      alpha1[,random.seq] = (1 - s) * alpha0[,random.seq] -
        s*(foreach(i = random.seq,.combine = cbind) %do% grad.diff3M(x.mat,u,strata = strata,i))
      #print( alpha1[,random.seq])
      bar_alpha = bar_alpha + 1/num.sample*apply(alpha1[,random.seq] - alpha0[,random.seq],1,sum)

      beta = (1 - s)*beta + s*num.sample*bar_alpha/P

      tot.it = tot.it + 1

    }

    epoch = epoch + 1
  }
  #print(cost(beta,P,x.mat))
  beta
}


# EM algorithm
myEM = function(Sub_Gen, Nfam, theta0, gth0, pi0, tol, maxit, S, rare.map, strata, seeds, common.inclusion, rare.inclusion, common.exclusion, rare.exclusion){

  L = ncol(Sub_Gen) - S                    #  number of rare variants
  region = length(unique(rare.map$gene))   #  number of region
  beta_old  = matrix(rep(0.5,S), S, 1)     # coefficients of common variants
  alpha_old = matrix(rep(0.5,L), L, 1)     # coefficients of rare variants

  b1 = S
  b2 = region
  b3 = L
  p_star0 = c(rep(0,S))                    # Current probability of inclusion for common variants
  g_star0 = c(rep(0,region))               # current probability of inclusion for group
  q_star0 = c(rep(0,L))                    # current probability of inclusion for rare variants

  # begin loop to run EM algorthm
  t = 0.1                      # deterministic anealing
  counter = 0
  while (counter <= maxit & t <= 1) {

    counter <- counter + 1
    P = c(inclufun_common(beta0 = beta_old,theta0 = theta0,t = t, common.inclusion = common.inclusion, common.exclusion = common.exclusion)$P,
      inclufun_rare(alpha0 = alpha_old, pi0 = pi0, t = t, rare.inclusion = rare.inclusion, rare.exclusion = rare.exclusion)$Qr)
    p_star = inclufun_common(beta0 = beta_old,theta0 = theta0,t = t, common.inclusion = common.inclusion, common.exclusion = common.exclusion)$p_star
    g_star = inclufun_group(alpha0 = alpha_old,gth0 = gth0,pi0 = pi0,region = region,t = t,
      rare.map = rare.map, rare.inclusion = rare.inclusion, rare.exclusion = rare.exclusion)
    q_star = inclufun_rare(alpha0 = alpha_old, pi0 = pi0, t = t,rare.inclusion = rare.inclusion, rare.exclusion = rare.exclusion)$q_star
    out_p =  as.matrix(c(p_star,q_star))
    rownames(out_p) = colnames(Sub_Gen)
    omega_new = SDCA1(num.var = ncol(Sub_Gen), num.sample = Nfam, x.mat = Sub_Gen, P = P, seeds = seeds[1:3], strata = strata)
    beta_new =  omega_new[1:S]
    alpha_new = omega_new[-c(1:S)]
    theta1 <- (sum(p_star) )/(S + b1 - 1)
    gth1 = (sum(g_star))/(region + b2 - 1)

    Z1 = diag(1,region,region)                                # max(rare_map$gene) by max(rare_map$gene) diagonal
    Z = Z1[rep(seq_len(nrow(Z1)), table(rare.map$gene)),]     # number of rare variants by number of genes

    pi1 = (sum(g_star*(t(Z) %*% q_star)))/(sum(as.vector(g_star)*table(rare.map$gene)) + b3 - 1)
    #print(pi1)

    #  Make decision on covergence using difference in likelihoods.
    # Finds Q and R at current and previous steps

    test <- Qlogistic(Sub_Gen,beta_new,alpha_new, theta1,gth1, pi1, p_star0,g_star0, q_star0,
      strata, common.exclusion, common.inclusion, rare.exclusion, rare.inclusion, rare.map, b1, b2, b3, Z) -
      Rlogistic(p_star0,g_star0,q_star0,theta1,gth1, pi1,rare.map, b1, b2, b3, Z)
    test2 <- Qlogistic(Sub_Gen,beta_old,alpha_old, theta0,gth0,pi0, p_star0,g_star0,q_star0,
      strata, common.exclusion, common.inclusion, rare.exclusion, rare.inclusion,rare.map, b1, b2, b3, Z) -
      Rlogistic(p_star0,g_star0, q_star0,theta0,gth0, pi0, rare.map, b1, b2, b3, Z)
    #print(abs(test-test2))
    # Convergence diagnosis
    if (abs(test - test2) < tol)
    {                                   # calculate the l1 norm


      beta_old = beta_new
      alpha_old = alpha_new
      theta0 = theta1
      pi0 = pi1
      gth0 = gth1
      p_star0 = p_star
      g_star0 = g_star
      q_star0 = q_star
      t = t + 0.1

    } else {


      beta_old = beta_new
      alpha_old = alpha_new
      theta0 = theta1
      pi0 = pi1
      gth0 = gth1
      p_star0 = p_star
      g_star0 =  g_star
      q_star0 = q_star
    }

  }

  if (t > 1)
  {
    cat("\nSuccessfully Converged\n")
    return(cbind(c(beta_new,alpha_new),out_p))
  } else{
    print("Convergence Failed")
    return(cbind(c(beta_new,alpha_new),out_p))
  }

}


rv_emvs <- function(prep.trio, maf, map, common.inclusion = 0.1, rare.inclusion = 0.5,
  common.exclusion = (log(sqrt(1.015/0.985))/qnorm(.975))^2,
  rare.exclusion = (log(sqrt(1.01/0.99))/qnorm(.975))^2,
  threshold = 0.05, theta0 = 0.5,        # the common  inclusion probability
  gth0 = 0.5,
  pi0 = 0.5                           # individual rare variants inclusion probability
) {
  # function testing:
  #  prep.trio = prep.list[[1]]
  #  maf = prep.list[[2]]
  #  map = read.csv('map.example.csv')
  #  common.inclusion = 0.1
  #  rare.inclusion = 0.5
  #  common.exclusion.range = c((log(sqrt(1.005/0.995))/qnorm(.975))^2,(log(sqrt(1.015/0.985))/qnorm(.975))^2)
  #  common.regular.number = 20
  #  threshold = 0.05; theta0 = 0.5        # the common  inclusion probability
  ## gth0 = 0.5
  #  pi0 = 0.5

  # separate common and rare variants by MAF, default for distinguish common and rare is maf >= 0.05
  common_variants = common_seperate(prep.trio = prep.trio,maf = maf, threshold = threshold)
  # re-arrange the order the columns
  eff_gen = cbind(prep.trio[,common_variants],prep.trio[,colnames(prep.trio)[which(colnames(prep.trio) %in% common_variants == F)]])
  # label the gene(region) only for rare variants
  rare_map = map[which(map$SNP %in% colnames(eff_gen)[(length(common_variants) + 1):ncol(eff_gen)]),]


  ## general initial parameters
  Nfam = nrow(prep.trio)/3
  strata = rep(1:Nfam,each = 3)

  set.seed(606)
  seeds =  sample(200:3000, 75,replace = FALSE)

  result = myEM(Sub_Gen = as.matrix(eff_gen),Nfam = Nfam, theta0 = theta0, gth0 = gth0, pi0 = pi0, tol = 1,
                maxit = 100, S = length(common_variants), rare.map = rare_map, strata = strata, seeds = seeds[1:3],
                common.inclusion, rare.inclusion, common.exclusion, rare.exclusion)

  colnames(result) = c('Coefficient','Probability')
  result = as.data.frame(result)

  return(result)
}