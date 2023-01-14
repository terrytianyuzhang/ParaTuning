rm(list=ls()); gc()
options(stringsAsFactors = F)



load(file="simulations-Results-6xx-7xx-8xx_01-13-2023.RData")

series=700  # specify series 600, 700, 800
TST="YRI"   # specify trest ancestry

# determine cis and trans ancestry (first one is cis)
TRN=c(TST,c("CEU","YRI")[c("CEU","YRI") != TST])

# collect AUC for all tests
AUC.df=NULL
for(sim in 0:9){
  sim=paste0("Sim-",series+sim)
  if(TST == "CEU"){
    AUC.df=rbind.data.frame(AUC.df,
                            data.frame(pnt.cis=RE[[sim]][["pnt"]][RE[[sim]][["pnt"]][,"TRN"] == TRN[1] & RE[[sim]][["pnt"]][,"TST"]==TST,"auc"],
                                       pnt.trans=RE[[sim]][["pnt"]][RE[[sim]][["pnt"]][,"TRN"] == TRN[2] & RE[[sim]][["pnt"]][,"TST"]==TST,"auc"],
                                       wpnt=RE[[sim]][["wpnt"]][RE[[sim]][["wpnt"]][,"TST"]==TST,"auc"],
                                       ls.cis=RE[[sim]][["ls"]][RE[[sim]][["ls"]][,"TRN"] == TRN[1] & RE[[sim]][["ls"]][,"TST"]==TST,"auc"],
                                       ls.trans=RE[[sim]][["ls"]][RE[[sim]][["ls"]][,"TRN"] == TRN[2] & RE[[sim]][["ls"]][,"TST"]==TST,"auc"],
                                       wls=RE[[sim]][["wls"]][RE[[sim]][["wls"]][,"TST"]==TST,"auc"],
                                       jls=RE[[sim]][["jls"]][RE[[sim]][["jls"]][,"TST"]==TST,"auc"]
                            )
    )
  }else{
    AUC.df=rbind.data.frame(AUC.df,
                            data.frame(pnt.cis=RE[[sim]][["pnt"]][RE[[sim]][["pnt"]][,"TRN"] == TRN[1] & RE[[sim]][["pnt"]][,"TST"]==TST,"auc"],
                                       pnt.trans=RE[[sim]][["pnt"]][RE[[sim]][["pnt"]][,"TRN"] == TRN[2] & RE[[sim]][["pnt"]][,"TST"]==TST,"auc"],
                                       wpnt=RE[[sim]][["wpnt"]][RE[[sim]][["wpnt"]][,"TST"]==TST,"auc"],
                                       ls.cis=RE[[sim]][["ls"]][RE[[sim]][["ls"]][,"TRN"] == TRN[1] & RE[[sim]][["ls"]][,"TST"]==TST,"auc"],
                                       ls.trans=RE[[sim]][["ls"]][RE[[sim]][["ls"]][,"TRN"] == TRN[2] & RE[[sim]][["ls"]][,"TST"]==TST,"auc"],
                                       wls=RE[[sim]][["wls"]][RE[[sim]][["wls"]][,"TST"]==TST,"auc"],
                                       jls=RE[[sim]][["jls"]][RE[[sim]][["jls"]][,"TST"]==TST,"auc"],
                                       tl=RE[[sim]][["tl"]][RE[[sim]][["tl"]][,"TST"]==TST,"auc"]
                            )
    )
    
  }
}

# mean and standard error
(mn=round(colMeans(as.matrix(AUC.df)),5))
(se=round(matrixStats::colSds(as.matrix(AUC.df))/sqrt(nrow(AUC.df)),5))

# collect information on the hyper parameters
HYPER.df=NULL
for(sim in 0:9){
  sim=paste0("Sim-",series+sim)
  HYPER.df=rbind.data.frame(HYPER.df,
                          data.frame(pnt.r2.cis=RE[[sim]][["pnt"]][RE[[sim]][["pnt"]][,"TRN"] == TRN[1] & RE[[sim]][["pnt"]][,"TST"]==TST,"r2"],
                                     pnt.r2.trans=RE[[sim]][["pnt"]][RE[[sim]][["pnt"]][,"TRN"] == TRN[2] & RE[[sim]][["pnt"]][,"TST"]==TST,"r2"],
                                     pnt.p.cis=RE[[sim]][["pnt"]][RE[[sim]][["pnt"]][,"TRN"] == TRN[1] & RE[[sim]][["pnt"]][,"TST"]==TST,"pvalue"],
                                     pnt.p.trans=RE[[sim]][["pnt"]][RE[[sim]][["pnt"]][,"TRN"] == TRN[2] & RE[[sim]][["pnt"]][,"TST"]==TST,"pvalue"],
                                     wpnt.w.trn=RE[[sim]][["wpnt"]][RE[[sim]][["wpnt"]][,"TST"]==TST,paste0(TST,".wght")],
                                     ls.shrink.cis=RE[[sim]][["ls"]][RE[[sim]][["ls"]][,"TRN"] == TRN[1] & RE[[sim]][["ls"]][,"TST"]==TST,"s"],
                                     ls.shrink.trans=RE[[sim]][["ls"]][RE[[sim]][["ls"]][,"TRN"] == TRN[2] & RE[[sim]][["ls"]][,"TST"]==TST,"s"],
                                     ls.lambda.cis=RE[[sim]][["ls"]][RE[[sim]][["ls"]][,"TRN"] == TRN[1] & RE[[sim]][["ls"]][,"TST"]==TST,"lambda"],
                                     ls.lambda.trans=RE[[sim]][["ls"]][RE[[sim]][["ls"]][,"TRN"] == TRN[2] & RE[[sim]][["ls"]][,"TST"]==TST,"lambda"],
                                     wls.w.trn=RE[[sim]][["wls"]][RE[[sim]][["wls"]][,"TST"]==TST,paste0(TST,".wght")],
                                     jls.shrink=RE[[sim]][["jls"]][RE[[sim]][["jls"]][,"TST"]==TST,"shrink"],
                                     jls.lambda=RE[[sim]][["jls"]][RE[[sim]][["jls"]][,"TST"]==TST,"lambda"],
                                     jls.gamma=RE[[sim]][["jls"]][RE[[sim]][["jls"]][,"TST"]==TST,"gamma"]
                          )
  )
}
(mn=round(colMeans(as.matrix(HYPER.df)),5))
(se=round(matrixStats::colSds(as.matrix(HYPER.df))/sqrt(nrow(HYPER.df)),5))

