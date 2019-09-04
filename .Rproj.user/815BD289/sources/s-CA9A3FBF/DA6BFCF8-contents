



par(mfrow = c(1, 1))
score_perso <- function(dag_reference,dag_mmhc,nb_variables) {
  scores_results=list("count"=list("tp"=0,"fp"=0,"tn"=0,"fn"=0),"recall"=0,"precision"=0,"fscore"=0,"dist2opt"=0)
  resultats_comparaison=compare(dag_reference,dag_mmhc)
  
  scores_results$count$tp=resultats_comparaison$tp
  scores_results$count$fp=resultats_comparaison$fp
  scores_results$count$fn=resultats_comparaison$fn
  
  scores_results$count$tn=(nb_variables*(nb_variables-1))-(resultats_comparaison$tp+resultats_comparaison$fp+resultats_comparaison$fn)
  scores_results$recall=scores_results$count$tp/(scores_results$count$tp+scores_results$count$fn)
  scores_results$precision=scores_results$count$tp/(scores_results$count$tp+scores_results$count$fp)
  scores_results$fscore=(2*scores_results$recall*scores_results$precision)/(scores_results$recall+scores_results$precision)
  scores_results$dist2opt=sqrt((1-scores_results$precision)^2+ (1-scores_results$precision)^2)
  return (scores_results)
  }

#main use: compare recall and precision of two algorithms

# load the data.
data(alarm)


# create and plot the true network structure.
modelstring = paste0("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF][LVF]",
                     "[STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA][HRSA|ERCA:HR][ANES]",
                     "[APL][TPR|APL][ECO2|ACO2:VLNG][KINK][MINV|INT:VLNG][FIO2][PVS|FIO2:VALV]",
                     "[SAO2|PVS:SHNT][PAP|PMB][PMB][SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC]",
                     "[MVS][VMCH|MVS][VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG]",
                     "[ACO2|VALV][CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]")
dag_reference = model2network(modelstring)
graphviz.plot(dag_reference)

#on determine les probabilites du dag en fonction du jeu de donnees
donnnes_ajustees=bn.fit(dag_reference,alarm,method="bayes")

#from the true structure, create a dataset of 100 000 lines, and save it
training_set=rbn(dag_reference,100000,alarm)
#potentially, impute the missing rows
training_set_without_missingvalue = impute(object=donnnes_ajustees, data=training_set)

#compare performances of several learning algorithms with the true real structure
#le premier element de la fonction compare correspond au graphe standard, attendu
#mmhc
par(mfrow = c(1, 2))
nb_variables=dim(training_set_without_missingvalue)[2]
dag_mmhc=mmhc(training_set_without_missingvalue)
graphviz.compare(dag_reference,dag_mmhc)
score_perso(dag_reference,dag_mmhc,nb_variables)

#rsamx2
dag_rsmax2=rsmax2(training_set_without_missingvalue)
graphviz.plot(dag_rsmax2)
score_perso(dag_reference,dag_rsmax2,nb_variables)
par(mfrow = c(1, 2))
graphviz.compare(dag_reference,dag_rsmax2)

#inter-iamb
dag_inter_iamb=inter.iamb(training_set_without_missingvalue)
graphviz.plot(dag_inter_iamb)
score_perso(dag_reference,dag_inter_iamb,nb_variables)
par(mfrow = c(1, 2))
graphviz.compare(dag_reference,dag_inter_iamb)


#mmpc
dag_mmpc=inter.iamb(training_set_without_missingvalue)
graphviz.plot(dag_mmpc)
score_perso(dag_reference,dag_mmpc,nb_variables)
par(mfrow = c(1, 3))
graphviz.compare(dag_reference,current=dag_mmpc)





setwd("C:/Users/Bastien/Documents/stage_projet_bayesian_network-20190516T081042Z-001/ancien_programme_r")
require('devtools')
install_github("https://github.com/bastienchassagnol/bnlearn-clone-3.4")
require(bnlearn)



training_set_alarm=read.csv(normalizePath("C:\\Users\\Bastien\\Documents\\stage_projet_bayesian_network-20190516T081042Z-001\\test_pyagrum\\algo_H2PC\\small_alarm.csv"))

for(i in 1:ncol(training_set_alarm)){
  training_set_alarm[,i] <- as.factor(training_set_alarm[,i])
}



hpc(training_set_alarm, test="x2",optimized = FALSE,pc.method="fdr.iamb")
#ci.test('HISTORY','LVFAILURE',data=training_set_alarm,test="x2")












x=c("foo", "bar", "foo", "foo", "bar", "bar")
y=c("one", "two", "one", "two", "one","two")
z1=c("dull", "shiny", "shiny", "dull", "shiny","shiny")
z2=c("hello","hello","hello","hello","mince","mince")
df=data.frame(x,y,z1,z2)
print(df['x'])

sum(is.na(df))
f <- function(x,y ,z=c(), df,test) {
  #convert dataframe as factor df
  for(i in 1:ncol(df)){
    df[,i] <- as.factor(df[,i])
  }
  #two types of tests, according whether z is present or not    
  if (is.null(z)) { resultat=ci.test(x=x,y=y,data=df,test=test) }
  else {resultat=ci.test(x=x,y=y,z=z,data=df,test=test) }                
  return (c(as.numeric(resultat$statistic),as.numeric(resultat$p.value)))                
}
f('x','y',c('z1','z2'),df,test="x2")



par(mfrow = c(2, 1))
dag1 = model2network("[A][B|A][C|A][D|C]")
dag2 = model2network("[A|B:C][B][C|D][D]")
graphviz.plot(dag1)
shd(cpdag(dag1),cpdag(dag2))



require(bnlearn)
dag_alarm_param=read.bif("alarm.bif")
graphviz.plot(dag_alarm_param)

data_random=rbn(dag_alarm_param,50000)
autre_dag=tabu(data_random,tabu=100,max.tabu=20)

dag_alarm=empty.graph(nodes(autre_dag))
amat(dag_alarm)=amat(dag_alarm_param)

hamming(cpdag(dag_alarm),cpdag(autre_dag))

autre_dag_fitted=bn.fit(autre_dag,data_random)
write.bif("random_dag.bif",autre_dag_fitted)

require(bnlearn)
data(learning.test)
configs(learning.test, all = TRUE)

data(learning.test)
res = gs(learning.test)
cpdag(res)
vstructs(res)



require(bnlearn)
bn=model2network("[A][D][B|A:C][C|D]")
bn2=model2network("[A][B|A][C|B][D|B]")
hamming(bn,bn2)
shd(bn,bn2)













