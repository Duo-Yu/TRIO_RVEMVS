
rv_emvs.regularization <- function(prep.trio, maf, map, common.inclusion = 0.1, rare.inclusion = 0.5,
  common.exclusion = NULL, common.exclusion.range,common.regular.number,
  rare.exclusion.range, rare.regular.number,
  threshold = 0.05, theta0 = 0.5,        # the common  inclusion probability
  gth0 = 0.5,
  pi0 = 0.5                           # individual rare variants inclusion probability
  ) {

  if (is.null(common.exclusion)) {
    exclusion = seq(common.exclusion.range[1], common.exclusion.range[2],length.out = common.regular.number)
    regular_result = colnames(prep.trio)
    for (V2 in exclusion) {
      V0 = V2
      result_temp = rv_emvs(prep.trio = prep.trio,maf = maf, map = map,
                            common.exclusion = V0,
                            rare.exclusion = V2,
                            threshold = 0.05, theta0 = 0.5,        # the common  inclusion probability
                            gth0 = 0.5,
                            pi0 = 0.5)

    #colnames(result_temp) = c('Coeffs','Probability')
    #result_temp = as.data.frame(result_temp)
    regular_result = cbind(regular_result,result_temp[,1])
  }
  colnames(regular_result) =  c("SNP",1:common.regular.number)
  plot_data = as.data.frame(cbind(exclusion,t(apply(regular_result[,-c(1)], 2, as.numeric))))
  colnames(plot_data) = c("exclusion.param",as.vector(regular_result[,1]))
  return(plot_data)
  } else{
    exclusion = seq(rare.exclusion.range[1], rare.exclusion.range[2],length.out = rare.regular.number)
    regular_result = colnames(prep.trio)
    for (V2 in exclusion) {
      result_temp = rv_emvs(prep.trio = prep.trio,maf = maf, map = map,
        common.exclusion = common.exclusion,
        rare.exclusion = V2,
        threshold = 0.05, theta0 = 0.5,        # the common  inclusion probability
        gth0 = 0.5,
        pi0 = 0.5)

      #colnames(result_temp) = c('Coeffs','Probability')
      #result_temp = as.data.frame(result_temp)
      regular_result = cbind(regular_result,result_temp[,1])
    }
    colnames(regular_result) =  c("SNP",1:common.regular.number)
    plot_data = as.data.frame(cbind(exclusion,t(apply(regular_result[,-c(1)], 2, as.numeric))))
    colnames(plot_data) = c("exclusion.param",as.vector(regular_result[,1]))
    return(plot_data)
  }
}

# library(ggplot2)
# library(latex2exp)
# melted  = melt(test_plot,id ='exclusion.param')
# ggplot(data = melted,aes(x = exclusion.param, y = value, group = variable)) + geom_point(size = 1) + geom_line() +
#   xlab(TeX("exclusion parameter")) + ylab(TeX("$\\hat{\\beta}$")) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


