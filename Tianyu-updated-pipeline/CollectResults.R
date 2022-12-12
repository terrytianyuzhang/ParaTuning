#rm(list=ls()); gc()
options(stringsAsFactors = F)
auc.dir <- paste0(main.dir,"Data/Simulation-AUC")
dir.create(auc.dir)

# i.sim <- 800
for(i.sim in 800:801){
  #### load the parameters for the simulation
  load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
  main.dir=params$run.info$main.dir
  work.dir=params$run.info$work.dir
  
  train.sets=names(params)[grep(".TRN",names(params))]
  tune.sets=names(params)[grep(".TUNE",names(params))]
  test.sets=names(params)[grep(".TST",names(params))]
  
  
  #######PT###########
  for(tune.set in tune.sets){
  # tune.set <- tune.sets[1]
  # train.dir=paste0(work.dir,train.set,"/")
  tune.dir=paste0(work.dir,tune.set,"/")
  # test.dir=paste0(work.dir,test.set,"/")
  
  cp.command <- paste0("cp ", 
                       tune.dir, "PRS/weightedPT-AUC.RData ", 
                       auc.dir, "/No", i.sim, tune.set, "weightedPT-AUC.RData")
  system(cp.command)
  }

  #######Lassosum###########
  for(test.set in test.sets){
    # tune.set <- tune.sets[1]
    # train.dir=paste0(work.dir,train.set,"/")
    # tune.dir=paste0(work.dir,tune.set,"/")
    test.dir=paste0(work.dir,test.set,"/")
    
    cp.command <- paste0("cp ", 
                         test.dir,"Lassosum/",test.set,"/WEIGHTED-",test.set,"-lassosum.AUC.RData ", 
                         auc.dir, "/No", i.sim, test.set, "weightedLassosum-AUC.RData")
    print(cp.command)
    system(cp.command)
  }

  #######JointLassosum######
  GAMMA <- c(0.2,0.5,0.8)
  for(gamma in GAMMA){
    for(test.set in test.sets){
      cp.command <- paste0("cp ", 
                           work.dir, "JointLassoSum/JointLassosum-",sprintf("-gamma-%.2f",gamma),"-", test.set,"-AUC.Rdata ", 
                           auc.dir, "/No", i.sim, sprintf("-gamma-%.2f",gamma),"-", test.set, "JointLassosum-AUC.RData")
      print(cp.command)
      system(cp.command)
    }
    
  }
  
}#i.sim
