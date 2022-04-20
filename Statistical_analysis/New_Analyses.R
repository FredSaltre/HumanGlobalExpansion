library(maps)
library(ade4)

trajm.dat <- read.csv("AveragedHumanTrajectories(Veghigh)_stats.csv",sep=",",dec=".",header=T)


### fonction qui retourne pour chaque trajectoire sa residuelle normalisee
### par la moyenne et l'ecart-type des trajectoires de meme debut-fin (donne par unique(trajm.dat[,2]))
norme.fct <- function(trajm.dat,v0)
## je rajoute trajm.dat dans les arguments
{
 v <- NULL
 for(i in unique(trajm.dat[,2]))
  {
  sel<- trajm.dat[,2]==i
  v1 <- v0[sel]
####     v <- norme.fct(trajm.dat[,15]*trajm.dat[,4]) rien ? foutre l? : ttention aux coupes colles
  v1<- (v1-mean(v1))/sd(v1)
  v <- c(v,v1)
}
return(v)
}


## test par permutation de l'egalite des moyennes des variables v1 et v2
##
test.fct <- function(v1,v2)
{
   stat.fct <-function(v1,v2)
    { return(mean(v1)-mean(v2))}
   
 obs <- stat.fct(v1,v2)
 sim <- NULL
 for(i in 1:10000)
  {
   w <- sample(c(v1,v2))
   sim <- c(sim,stat.fct(w[1:length(v1)],w[(length(v1)+1):length(w)]))
  }
return(sum(sim>obs)/length(sim))
}



### les variables en boucle
## et on refait pour toutes les variables d'environnement
out<- data.frame(matrix(ncol = 2, nrow = 0))
colnames(out) <- c('Variable', 'Probablilty')

for (i in 8:13) 
{
v <- norme.fct(trajm.dat,trajm.dat[,i]*trajm.dat[,4])
res <- test.fct(v[trajm.dat[,1]!="Y"],v[trajm.dat[,1]=="Y"])
out<-rbind(out,c(names(trajm.dat)[i],abs(res-0.5)))
print(c(names(trajm.dat)[i],"p-valeur de la diff de moyenne",abs(res-0.5)))
}

out$Probablilty<-as.numeric(out$Probablilty)
out$scaleprob<-(out$Probablilty - min(out$Probablilty)) / (max(out$Probablilty) - min(out$Probablilty)) 
colnames(out) <- c('Variable', 'Probablilty', 'Scaled Probablilty')
