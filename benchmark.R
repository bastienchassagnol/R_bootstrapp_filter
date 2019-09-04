# check.packages function: install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}



  compute_struct_vector<-function (golden_bn, bn_learnt, score_measured) {
    vector_score=c()
    num_nodes=length(nodes(golden_bn))
    
    compte=compare(golden_bn,bn_learnt)
    compte=c(compte,tn=num_nodes*(num_nodes-1)-compte$tp-compte$fp-compte$fn)
    compte=c(compte,recall=compte$tp/(compte$tp+compte$fn), precision=compte$tp/(compte$tp+compte$fp),compte$tn/(compte$tn+compte$fp))
    compte=c(compte,fscore=2*(compte$recall*compte$precision)/(compte$recall+compte$precision))
    compte=c(compte,dist2opt=sqrt((1-compte$precision)^2+ (1-compte$precision)^2))
    compte=c(compte,specificity=compte$tn/(compte$tn+compte$fp))
    
    for (struct_score in score_measured) {
      individual_score=switch(struct_score,
             "shd" = shd(golden_bn,bn_learnt),
             "hamming" = hamming(golden_bn,bn_learnt),
             "tp" = compte$tp,
             "fp" = compte$fp,
             "fn" = compte$fn,
             "tn" = compte$tn,
             "recall" = compte$recall,
             "precision"=compte$precision,
             "specificity"=compte$specificity,
             "fscore"=compte$fscore,
             "dist2opt"=compte$dist2opt,
      )
      vector_score=c(vector_score, individual_score) 
    }
    return (as.list(setNames(vector_score, score_measured)))
  }
  
  
  
  
  compute_proba_vector<- function (bn_learnt, data.set, score_measured,type="train") {
    
      name_score=c()
      vector_score=c()
      for (proba_score in score_measured) {
        
        individual_score=switch(strsplit(proba_score,split="_")[[1]][1],
                                "bic" = BIC(bn_learnt,data.set),
                                "aic"= AIC(bn_learnt,data.set),
                                "loglik" = score(bn_learnt,data.set,type="loglik"),
                                "bde" = score(bn_learnt,data.set,type="bde")
        )
        name_score=c(name_score,paste(proba_score,type,sep="_"))
        vector_score=c(vector_score, individual_score) 
      }
      temp=as.data.frame(as.list(setNames(vector_score, name_score)))
      
      return (as.data.frame(as.list(setNames(vector_score, name_score))))
  }
  
  
  compute_structural_scores<- function(golden_bn, size,restrict,score_measured ) {
    
    temp_database=rbn(golden_bn, n = size)
    start.time <- Sys.time()
    if (restrict=="hpc") {
      bn_learnt=h2pc(temp_database, maximize = "tabu", maximize.args = list(tabu = 100, max.tabu = 20),alpha=0.05,test="x2",pc.method="fdr.iamb")
          }
    else {
      bn_learnt=mmhc(temp_database, maximize = "tabu", maximize.args = list(tabu = 100, max.tabu = 20),alpha=0.05,test="x2" )
    }
    
    end.time <- Sys.time()
    time.taken <- as.numeric(end.time - start.time)
    structural_result=compute_struct_vector(bn.net(golden_bn),bn_learnt,score_measured)
    if ("time" %in% score_measured) {
      structural_result=c(structural_result,time=time.taken)
    }
    rm(temp_database)
    
    return (as.data.frame(structural_result))
  }
    
  
  compute_probalistic_scores<- function(golden_bn, size,restrict,score_measured,n_splits=10){
    
    #generate folds
    flds <- createFolds(seq(1,size,1), k = n_splits, list = TRUE, returnTrain = FALSE)
    #generate df
    temp_database=rbn(golden_bn, n = size)
  
    for (index_fold in seq_along(flds)) {
      train_data=temp_database[-flds[[index_fold]],]
      if (restrict=="hpc") {
        bn_learnt=h2pc(train_data, maximize = "tabu", maximize.args = list(tabu = 100, max.tabu = 20),alpha=0.05,test="x2",pc.method="fdr.iamb")
      }
      else {
        bn_learnt=mmhc(train_data, maximize = "tabu", maximize.args = list(tabu = 100, max.tabu = 20),alpha=0.05,test="x2" )
              }
      test_data=temp_database[flds[[index_fold]],]
      train_score=compute_proba_vector(bn_learnt,train_data,score_measured,type="train")
      test_score=compute_proba_vector(bn_learnt,test_data,score_measured,type="test")
      fused_score=cbind(train_score,test_score)
      if (index_fold==1) {
        #define df_proba
        df_proba=fused_score
      }
      else {
        df_proba=rbind(df_proba,fused_score) 
      }
    }
    rm(temp_database)
    return (as.data.frame(as.list(colMeans(df_proba))))
  }
    
  
  compute_average_distance<- function(df_scores,golden_bn, nsamples,size,algorithm,proba_score,struct_score,n_splits) {
    #matrix of shape (repetitions * scores) to get the mean of each type of score
    #algorithm:(type,agrs,kwargs)
    scoring_df=df_scores
    for (repetition in seq(1,nsamples,1)) {
      if (length(proba_score>0)) {
        temp_proba=compute_probalistic_scores(golden_bn, size,algorithm,proba_score, n_splits=n_splits)
        
      }
      
      if (length(struct_score>0)) {
        
        temp_struct=compute_structural_scores(golden_bn, size,algorithm,struct_score)
      }
      temp_score=cbind(temp_proba,temp_struct)
      scoring_df=rbind(scoring_df,temp_score)
    } 
      return (scoring_df)
  }
   
  
  learn_scores<- function(name_bn,sample_size,score_measured=c('dist2opt'),algorithms=c('hpc'),n_splits=10,nsamples=30) {
     
    structural_scores=c('recall','precision','fscore','dist2opt','specificity','hamming','shd','time','number_tests')
    proba_scores=c('bic','aic','loglik',"bde")
    
    bn=read.bif(file.path("true_graphes_structures",name_bn))
    struct.vector=c()
    proba.vector=c()
    for (score in score_measured) {
      if (!(score %in% c(structural_scores,proba_scores))) {
        stop(cat("distance score still not implemented, list of of possible computations is ",c(structural_scores,proba_scores)))
      }
      if (score %in% structural_scores)   struct.vector=c(struct.vector,score)
      if (score %in% proba_scores)      proba.vector=c(proba.vector,score)
    }
    struct.vector=sort(struct.vector)
    
    proba_complet=c()
    for (score in proba.vector) {
      proba_complet=c(proba_complet,paste0(score,"_train"))
      proba_complet=c(proba_complet,paste0(score,"_test"))
    }
    
    restrain_algorithms=c('si.hiton.pc','hpc', 'mmpc', 'iamb.fdr', 'gs','pc.stable')
    #assert that algorithms are computed
    for (algo in algorithms) {
      if (! (algo %in% restrain_algorithms)) {
        stop("algorithm still not implemented")
      }
    } 
     
    name_colonnes=c(order(proba_complet),struct.vector)
    #matrix scoring all scores measured for each database
    df_score_final=data.frame(matrix(nrow = 0, ncol = length(name_colonnes)))
    colnames(df_score_final) <- name_colonnes
    for (size in sample_size) {
      cat("\n We are at size ", size,"\n")
      for (algorithm in algorithms) {
        cat("nous en sommes a l'algo ", algorithm, "pour la taille suivante de database ", size,"\n")
        df_score_final=compute_average_distance(df_score_final,bn, nsamples,size,algorithm,proba_score=proba.vector,struct_score=struct.vector,n_splits=n_splits)
      }
    }
    name_file=paste('scores_total_',sub("^([^.]*).*", "\\1", name_bn),"csv",sep=".")
    write.csv(alarm,file=file.path("scores",name_file),row.names=FALSE,quote=FALSE)
    
  
    return (df_score_final)
  }    
  
  
  #install packages
  setwd("C:/Users/Bastien/Documents/stage_projet_bayesian_network-20190516T081042Z-001/ancien_programme_r")
  packages<-c("caret","parallel","snow","devtools")
  check.packages(packages)
  require(bnlearn)
  
  liste_graphes=list.files(path = "true_graphes_structures")
  score_measured=c('recall','precision','dist2opt','specificity','hamming','shd','aic','bic')
  algorithms=c('hpc','mmpc')
  sample_size=c(500,1000,2000,5000,10000,20000,50000)
  for (graph in liste_graphes) {
    cat("we learn following graph ", graph)
    learn_scores(graph,sample_size,score_measured=score_measured,algorithms=algorithms,n_splits=10,nsamples=30)
  }
  
  
  
  

  
  
  
  
  