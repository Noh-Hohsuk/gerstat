#' Conduct a range of germination rate comparison tests.
#'
#' The function executes a germination rate comparison test based on
#' the specified method. If the difference in germination rate comparison is
#'  statistically significant at the designated significance level, it generates a compact
#'  letter display that groups treatments with similar germination rates,
#'  utilizing the specified multiple comparison method.
#'
#' @param gdata A data-frame containing the germination data
#' ('gdata,' is expected to be a data frame containing three variables:
#' 'treatment,' 'num_seed,' and 'ger_seed'. The 'treatment' variable is a factor
#' variable representing the treatment level, 'num_seed' should indicate
#' the number of seeds used in each iteration of the treatment,
#' and 'ger_seed' should indicate the number of seeds that germinated.)
#'
#' @param method A statistical method for germination rate comparison
#' (The options for method are "LRT" (likelihood ratio test), "ANOVA"
#' (analysis of variance), "AS-ANOVA" (ANOVA with arcsine square root transformed response),
#'  "logistic" (logistic regression), or "KW" (Kruskal-Wallis test).
#'
#' @param p_adjust_mtd P-value correction method for multiple comparison
#' The options for the correction method are
#' c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none").
#' For ANOVA with Tukey's HSD, the method should be "single-step"
#'
#' @param ctr_lv significance level for overall test and control level for
#' multiple testing (we assume that both levels are the same)
#' @return p_value the resulting P-value from the test
#' @return cld compact letter display from the multiple testing
#' (When the p-value is smaller than the specified significance level,
#'  the cld becomes NULL)
#' @examples
#' # Data generation
#' data(nandina)
#' ger_test(gdata=nandina,method="LRT")
#' @export
ger_test<-function(gdata=gdata,method,p_adjust_mtd="holm",ctr_lv=0.05){
  if (!is.factor(gdata$treatment)){
    gdata$treatment<-as.factor(gdata$treatment)
  }

  if (method=="LRT"){
    return(test_pANOVA_multiple(gdata=gdata,p_adjust_mtd=p_adjust_mtd,ctr_lv=ctr_lv))
  }

  if (method=="ANOVA"){
    return(test_ANOVA(gdata=gdata,p_adjust_mtd=p_adjust_mtd,ctr_lv=ctr_lv))
  }

  if (method=="AS_ANOVA"){
    return(test_AS_ANOVA(gdata=gdata,p_adjust_mtd=p_adjust_mtd,ctr_lv=ctr_lv))
  }

  if (method=="logistic"){
    return(test_logistic(gdata=gdata,p_adjust_mtd=p_adjust_mtd,ctr_lv=ctr_lv))
  }

  if (method=="KW"){
    return(test_KW(gdata=gdata,p_adjust_mtd=p_adjust_mtd,ctr_lv=ctr_lv))
  }

}

test_pANOVA<-function(gdata){

  p_null<-sum(gdata$ger_seed)/sum(gdata$num_seed)
  p_alter_mat<- dplyr::summarise(dplyr::group_by(gdata,treatment),g_ratio = sum(ger_seed)/sum(num_seed))
  level_vec<-levels(gdata$treatment) # level vector of the treatment

  # negative log likelihood under the null hypothesis
  nloglik_null<-0
  for (i in level_vec){
    temp_mat<-dplyr::select(dplyr::filter(gdata,treatment==i),num_seed,ger_seed)
    for (j in 1:dim(temp_mat)[1]){
      nloglik_null<-nloglik_null+(-1)*log(dbinom(temp_mat[j,2], size=temp_mat[j,1], prob=p_null) )
    }
  }

  # negative log likelihood under the alternative hypothesis
  nloglik_alter<-0
  for (i in level_vec){
    temp_mat<- dplyr::select(dplyr::filter(gdata,treatment==i),num_seed,ger_seed)
    p_alter_i<- dplyr::pull(dplyr::select(dplyr::filter(p_alter_mat,treatment==i),g_ratio))
    for (j in 1:dim(temp_mat)[1]){
      nloglik_alter<-nloglik_alter+(-1)*log(dbinom(temp_mat[j,2], size=temp_mat[j,1], prob=p_alter_i))
    }
  }

  stat<-2*(nloglik_null-nloglik_alter)
  df_test <- length(dplyr::pull(dplyr::select(p_alter_mat,g_ratio)))-1
  return(1-pchisq(stat, df_test, ncp = 0, lower.tail = TRUE, log.p = FALSE))
}


test_pANOVA_multiple<-function(gdata=gdata,p_adjust_mtd,ctr_lv){

  pval<-test_pANOVA(gdata)
  if (pval>=ctr_lv) {
    return(list(p_value = pval,cld=NULL))
  }

  comp_mat<-combn(levels(gdata$treatment), 2)
  res_pval<-numeric(dim(comp_mat)[2])
  for (i in 1:dim(comp_mat)[2]){
    sub_data<-dplyr::filter(gdata,(treatment==comp_mat[1,i])|(treatment==comp_mat[2,i]))
    sub_data<-droplevels(sub_data)
    res_pval[i]<-test_pANOVA(sub_data)
  }


  p_adjust<-p.adjust(res_pval,p_adjust_mtd)

  lev_vec<-levels(gdata$treatment)
  dim_num<-length(lev_vec)-1
  pval_mat<-matrix(NA,dim_num,dim_num)

  for (j in 1:dim_num)
  {
    f_ind<-1+(dim_num+1)*(j-1)-(j-1)*j/2
    l_ind<-(dim_num+1)*j-j*(j+1)/2
    pval_mat[(j:dim_num),j]<-p_adjust[f_ind:l_ind]
  }

  dimnames(pval_mat)<-list(lev_vec[2:(dim_num+1)],lev_vec[1:(dim_num)])
  lvl_order<-levels(gdata$treatment)
  comps <-na.omit(tidyr::pivot_longer(tibble::as_tibble(pval_mat, rownames="row"),-row, names_to="col", values_to="p"))
  mycomps <- as.matrix(comps[,1:2])
  signif <- comps$p < ctr_lv
  cld_res<-multcomp:::insert_absorb(signif,decreasing=FALSE,comps=mycomps,lvl_order=lvl_order)

  return(list(p_value = pval,cld=cld_res$Letters))
}

test_ANOVA<-function(gdata=gdata,p_adjust_mtd,ctr_lv){
  prop_ger<-gdata$ger_seed/gdata$num_seed
  gdata_ano<-cbind(gdata,prop_ger)
  anova_res<-aov(prop_ger~treatment,data=gdata_ano)
  pval<-summary(anova_res)[[1]][["Pr(>F)"]][1]

  if (pval>=ctr_lv) {
    return(list(p_value = pval,cld=NULL))
  }

  glht.contr <-multcomp::glht(anova_res,linfct = multcomp::mcp(treatment="Tukey"))
  cld_res<-multcomp::cld(summary(glht.contr, test=multcomp::adjusted(p_adjust_mtd)))
  return(list(p_value = pval,cld=cld_res))
}

test_AS_ANOVA<-function(gdata=gdata,p_adjust_mtd,ctr_lv){
  as_prop_ger<-asin(sqrt((gdata$ger_seed+3/8)/(gdata$num_seed+3/4)))
  gdata_as_ano<-cbind(gdata,as_prop_ger)
  as_anova_res<-aov(as_prop_ger~treatment,data=gdata_as_ano)
  pval<-summary(as_anova_res)[[1]][["Pr(>F)"]][1]

  if (pval>=ctr_lv) {
    return(list(p_value = pval,cld=NULL))
  }

  glht.contr <-multcomp::glht(as_anova_res,linfct = multcomp::mcp(treatment="Tukey"))
  cld_res<-multcomp::cld(summary(glht.contr, test=multcomp::adjusted(p_adjust_mtd)))
  return(list(p_value = pval,cld=cld_res))
}

test_logistic<-function(gdata=gdata,p_adjust_mtd,ctr_lv){
  bmod <- glm(cbind(ger_seed,num_seed-ger_seed) ~ treatment, family=binomial,gdata)
  stat<-bmod$null.deviance-bmod$deviance
  pval<-1-pchisq(stat,length(bmod$coefficients)-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)

  if (pval>=ctr_lv) {
    return(list(p_value = pval,cld=NULL))
  }

  glht.contr <-multcomp::glht(bmod,linfct = multcomp::mcp(treatment="Tukey"))
  cld_res<-multcomp::cld(summary(glht.contr, test=multcomp::adjusted(p_adjust_mtd)))
  return(list(p_value = pval,cld=cld_res))
}

test_KW<-function(gdata=gdata,p_adjust_mtd,ctr_lv){
  prop_ger<-gdata$ger_seed/gdata$num_seed
  gdata_KW<-cbind(gdata,prop_ger)
  DT=kruskal.test(prop_ger ~ treatment, data = gdata_KW)
  pval<-DT$p.value

  if (pval>=ctr_lv) {
    return(list(p_value = pval,cld=NULL))
  }

  DT_multi = PMCMRplus::kwAllPairsDunnTest(prop_ger ~ treatment, data=gdata_KW,  p.adjust.method = p_adjust_mtd)
  DTT =rcompanion::PMCMRTable(DT_multi)

  comps <-na.omit(tidyr::pivot_longer(tibble::as_tibble(DT_multi$p.value, rownames="row"),-row, names_to="col", values_to="p"))
  mycomps <- as.matrix(comps[,1:2])
  signif <- comps$p < .05

  cld_res<-multcomp:::insert_absorb(signif,decreasing=FALSE,comps=mycomps,lvl_order=levels(gdata$treatment))
  return(list(p_value = pval,cld=cld_res$Letters))
}

