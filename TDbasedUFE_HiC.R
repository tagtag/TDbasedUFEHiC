TDbasedUFE_HiC <- function(file_Matrix,file_bed,l_list=c(1),func_list=list(CTCF=T,PLS=T,pELS=T,dELS=T),
                           sd0=0.01,breaks=100,sel=F)
{
  require(Matrix)
  require(irlba)
  require(rtracklayer)
  suppressMessages(library(AnnotationHub))
  th <- function(sd,breaks=100){
    P2<- pchisq(((SVD$u[,k0]-mean(SVD$u[,k0]))/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=breaks,plot=F)
    return(sd(hc$count[1:sum(hc$mid<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
  }
  X <- NULL
  for (i in seq_along(file_Matrix))
  {
    cat(i," ")
    x1 <- read.csv(file_Matrix[i],sep="\t",header=F)
    xx1 <- sparseMatrix(i=x1[,1],j=x1[,2],x=x1[,3],dims=c(max(x1[,1:2]),max(x1[,1:2])))
    xx1 <- xx1/sum(xx1)
    xx1 <- xx1+t(xx1) #ここを忘れていた。
    #dim(xx1) <- c(2806*2806,1)
    X <- cbind(X,xx1)
  }
  SVD <- irlba(X,max(l_list))
  aggregate_scores <- rep(list(NA),4)
  GR <- import(file_bed,format="bed")
  if (func_list[[1]]) #CTCF
  {
    ah <- AnnotationHub()
    query_data <- subset(ah, preparerclass == "CTCF")
    CTCF_hg38_all <- query_data[["AH104727"]]
    overlaps <- findOverlaps(GR, CTCF_hg38_all)
    scores <- score(CTCF_hg38_all)[subjectHits(overlaps)]
    query_indices <- queryHits(overlaps)
    aggregate_scores[[1]] <- tapply(scores, factor(query_indices, levels = seq_along(GR)), sum, default = 0)
  }
  if (func_list[[2]] | func_list[[3]] | func_list[[4]])
  {
    ENCODE <- read.csv("41586_2020_2493_MOESM12_ESM.txt",sep="\t")
    write.table(ENCODE[,1:3],file="ENCODE.bed",col.names=F,row.names=F,quote=F,sep="\t")
    ENCODE2 <- import("ENCODE.bed",format="bed")
    for (i in c(2:4))
    if (func_list[[i]])
    {
      if (i==2)
      {
        index <- grep("PLS",ENCODE[,6])
      } else if (i==3)
      {
        index <- grep("pELS",ENCODE[,6])
      } else if (i==4)
      {
        index <- grep("dELS",ENCODE[,6])
      }
      overlaps <- findOverlaps(GR, ENCODE2[index])
      scores <- rep(1,length(subjectHits(overlaps)))
      query_indices <- queryHits(overlaps)
      aggregate_scores[[i]] <- tapply(scores, factor(query_indices, levels = seq_along(GR)), sum, default = 0)
    }
  }
  func_name <- c("CTCF","PLS","pELS","dELS")
  COR <- list(CTCF=rep(list(list(Pearson=list(NA),Spearman=list(NA))),max(l_list)+1),
              PLS=rep(list(list(Pearson=list(NA),Spearman=list(NA))),max(l_list)+1),
              pELS=rep(list(list(Pearson=list(NA),Spearman=list(NA))),max(l_list)+1),
              dELSF=rep(list(list(Pearson=list(NA),Spearman=list(NA))),max(l_list)+1)
              )
  print(str(COR))
  for (logic in 1:4)
  {
    if (func_list[[logic]])
    {
      COR[[logic]][[1]][[1]] <- cor.test(aggregate_scores[[logic]],rowSums(X))
      COR[[logic]][[1]][[2]] <- cor.test(aggregate_scores[[logic]],rowSums(X),method="spearman")
      for (l in l_list)
      {
        cat(l," ")
        COR[[logic]][[1+l]][[1]] <- cor.test(aggregate_scores[[logic]],SVD$u[,l])
        COR[[logic]][[1+l]][[2]] <- cor.test(aggregate_scores[[logic]],SVD$u[,l],method="spearman")
      }
    }
  }
  MEANSD <- rep(NA,max(l_list))
  if (sel)
  {
    for (l in l_list)
    {
      k0<-l
      P1 <- pchisq(((SVD$u[,k0]-mean(SVD$u[,k0]))/sd0)^2,1,lower.tail=F)
      sd <- optim(sd0,th)$par
      P1<- pchisq(((SVD$u[,k0]-mean(SVD$u[,k0]))/sd)^2,1,lower.tail=F)
      #cat("l=",l," sd=",sd," sum=",sum(p.adjust(P1,"BH")<0.01),"\n")
      if (sum(p.adjust(P1,"BH")<0.01)==0) 
        {
        cat("*** No bins are selected  for",l,"th PC. change sd value *** \n")
        break
        }
      pdf(file=paste("hist_",l,".pdf",sep=""))
      hist(1-P1,breaks=breaks)
      dev.off()
      pdf(file=paste("plot_",l,".pdf",sep=""))
      aa <- seq(0.5*sd,2*sd,by=0.1*sd)
      bb<-apply(matrix(seq(0.5*sd,2*sd,by=0.1*sd),ncol=1),1,th)
      plot(aa,bb,xlab="sigma_l",ylab="sigma_h",type="o")
      arrows(sd,max(bb),sd,min(bb),col=2)
      dev.off()
      pdf(file=paste("image_",l,".pdf",sep=""),width=10,height=5)
      print(image(X[p.adjust(P1,"BH")<0.01,p.adjust(P1,"BH")<0.01]))
      dev.off()
    #
      XX <- X[p.adjust(P1,"BH")<0.01,p.adjust(P1,"BH")<0.01]
      dim <- c(dim(XX)[1],dim(XX)[2]/length(file_Matrix),length(file_Matrix))
      XX <- array(XX,dim)
      MEANSD[l] <- mean(apply(XX,1:2,sd))/mean(XX)
      #
      pdf(file=paste("SVD_",l,".pdf",sep=""),width=10,height=10)
      plot(log10(abs(SVD$u[,k0])),pch=16,cex=0.5,col=c(p.adjust(P1,"BH")<0.01)+1)
      dev.off()
      write.table(file=paste("index_",l,".csv",sep=""),data.frame(which(p.adjust(P1,"BH")<0.01)),sep="\t",col.names=F,row.names=F)
      P_value <- cbind(P1,p.adjust(P1,"BH"))
      colnames(P_value) <- c("P","adjusted P")
      save(file=paste("P_",l,sep=""),P_value)
      }
   }
  return(list(COR=COR,SVD=SVD,MEANSD=MEANSD))
}