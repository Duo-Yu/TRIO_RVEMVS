# this file used to pre-porcess the trio data set for rv_emvs
rv_emvs.prepare = function(file.name){
  # file.name = 'example.trios.csv'
  require('data.table')
  require('foreach')
  trio.data = fread(file.name)
  trio.data = as.data.frame(trio.data)
  Nfam = nrow(trio.data)/6
  Hapid = trio.data$Hap_id
  # get cases and pesudo controls: change from 6 rows for each family to 8 rows.
  expand_hapid = c()
  for (i in 1:Nfam) {
    expand_hapid[8*(i - 1) + 1] = Hapid[6*(i - 1) + 5]          # haplotype 1 for the case child
    expand_hapid[8*(i - 1) + 2] = Hapid[6*(i - 1) + 6]          # haplotype 2 for the case child

    expand_hapid[8*(i - 1) + 3] = Hapid[6*(i - 1) +
                            which(Hapid[c(6*(i - 1) + 1,6*(i - 1) + 2)] %in%
    expand_hapid[c(8*(i - 1) + 1,8*(i - 1) + 2)] == FALSE)]     # haplotype not transmitted to case

    expand_hapid[8*(i - 1) + 4] = Hapid[6*(i - 1) + 2 +
                            which(Hapid[c(6*(i - 1) + 3,6*(i - 1) + 4)] %in%
    expand_hapid[c(8*(i - 1) + 1,8*(i - 1) + 2)] == TRUE)]     # haplotype transmitted to case

    expand_hapid[8*(i - 1) + 5] = Hapid[6*(i - 1) +
                            which(Hapid[c(6*(i - 1) + 1,6*(i - 1) + 2)] %in%
    expand_hapid[c(8*(i - 1) + 1,8*(i - 1) + 2)] == FALSE)]   # haplotype not transmitted to case

    expand_hapid[8*(i - 1) + 6] = Hapid[6*(i - 1) + 2 +
                            which(Hapid[c(6*(i - 1) + 3,6*(i - 1) + 4)] %in%
    expand_hapid[c(8*(i - 1) + 1,8*(i - 1) + 2)] == FALSE)]  # haplotype not transmitted to case

    expand_hapid[8*(i - 1) + 7] = Hapid[6*(i - 1) +
                            which(Hapid[c(6*(i - 1) + 1,6*(i - 1) + 2)] %in%
    expand_hapid[c(8*(i - 1) + 1,8*(i - 1) + 2)] == TRUE)]

    expand_hapid[8*(i - 1) + 8] = Hapid[6*(i - 1) + 2 +
                            which(Hapid[c(6*(i - 1) + 3,6*(i - 1) + 4)] %in%
    expand_hapid[c(8*(i - 1) + 1,8*(i - 1) + 2)] == FALSE)]
  }

  cacon_hap = trio.data[match(expand_hapid,trio.data$Hap_id),]
  cacon_hap = as.data.frame(cacon_hap[,-c(1:5)])

  Ngen = nrow(cacon_hap)/2
  holder = c(rep(1:Ngen,each = 2))
  cacon_hap = cbind(holder,cacon_hap)
  cacon_gen = aggregate(cacon_hap, by = list(Category = cacon_hap$holder), FUN = sum)
  cacon_gen = cacon_gen[,colnames(cacon_gen) %in% c('Category','holder') == F]  # remove holder
  cacon_gen = as.matrix(cacon_gen)
  remove(cacon_hap)

  ##removing the non-polymorphic SNPs
  cov = colnames(cacon_gen[,which(apply(cacon_gen,2,sum) > 0)])
  eff_gen = cacon_gen[,cov]
  maf = apply(cacon_gen,2,sum)/4/Nfam/2
  ## constructing X matrix, each row of X is the subtraction of case and controls
  strata1 = rep(1:Nfam,each = 4)
  Sub_Gen = eff_gen[FALSE,]
  Sub_Gen = foreach(i = unique(strata1),.combine = rbind) %do% apply(eff_gen[strata1 == i,],2,
                        function(u) rbind(u[1] - u[2], u[1] - u[3], u[1] - u[4]))
  Sub_Gen = as.matrix(Sub_Gen)

  return(list(Sub_Gen,maf))
}

# test = rv_emvs.prepare('example.trios.csv')
