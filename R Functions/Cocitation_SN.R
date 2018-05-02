library(tm)
library(igraph)
library(rgexf)   
library(RColorBrewer)
library(ggplot2)
library(dplyr)

joint_cocitationnetwork<-function(df1,df2,df3){
  colnames(s1)<-c("w1.1","w2.1","cor.1")
  colnames(s2)<-c("w1.2","w2.2","cor.2")
  colnames(s3)<-c("w1.3","w2.3","cor.3")
  cw1<-inner_join(s1,s3,by=c("w1.1"="w1.3"))
  cw2<-inner_join(s2,s3,by=c("w1.2"="w1.3"))
  ##1
  cwu1.1<-as.data.frame(unique(cw1[,1]))
  cwu2.1<-as.data.frame(unique(cw1[,2]))
  cwu3.1<-as.data.frame(unique(cw1[,4]))
  colnames(cwu1.1)<-c("c_edges")
  colnames(cwu2.1)<-c("c_edges")
  colnames(cwu3.1)<-c("c_edges")
  cw1<-rbind(cwu1.1,cwu2.1,cwu3.1)
  cw1<-as.data.frame(unique(cw1[,1]))
  colnames(cw1)<-c("c_edges")
  ##2
  cwu1.2<-as.data.frame(unique(cw2[,1]))
  cwu2.2<-as.data.frame(unique(cw2[,2]))
  cwu3.2<-as.data.frame(unique(cw2[,4]))
  colnames(cwu1.2)<-c("c_edges")
  colnames(cwu2.2)<-c("c_edges")
  colnames(cwu3.2)<-c("c_edges")
  cw2<-rbind(cwu1.2,cwu2.2,cwu3.2)
  cw2<-as.data.frame(unique(cw2[,1]))
  colnames(cw2)<-c("c_edges")
  ##
  cw<-inner_join(cw1,cw2,by=c("c_edges"="c_edges"))
  cw$mark<-"CE"
  colnames(cw)<-c("c_edges","CE")
  
  return(cw)
}

all_cocverts<-function(df1,df2,m1,m2,m3){
  colnames(df1)<-c("w1","w2","cor")
  colnames(df2)<-c("w1","w2","cor")
  colnames(df3)<-c("w1","w2","cor")
  df1$mark<-m1
  df2$mark<-m2
  df3$mark<-m3
  df<-rbind(df1,df2,df3)
  df<-df[!duplicated(df[,c(1,2)]), ]
  dfa<-df[,c(1,4)]
  colnames(dfa)<-c("c_edges","mark")
  dfb<-df[,c(2,4)]
  colnames(dfb)<-c("c_edges","mark")
  dfa<-dfa[!duplicated(dfa[,c(1)]), ]
  dfb<-dfb[!duplicated(dfb[,c(1)]), ]
  df<-rbind(dfa,dfb)
  return(df)
}


##Order to runing the cocitations functions
setwd("/Volumes/HD/Archivos/Scripts/Python/TextAn")
s1<-get_semantic_network(".pdf",1:15)
s2<-get_semantic_network(".pdf",1:15)
s3<-get_semantic_network(".pdf",1:15)
j1<-joint_cocitationnetwork(s1,s2,s3)
f1<-all_cocverts(s1,s2,s3,"p88","s80","c07")
f2<-total_nodes(j1,f1)
colnames(s1)<-c("w1","w2","cor")
colnames(s2)<-c("w1","w2","cor")
colnames(s3)<-c("w1","w2","cor")
s<-rbind(s1,s2,s3)
g<-graph.data.frame(s,directed = F,vertices = f2)
g<-simplify(g,remove.multiple = T,remove.loops = F)
setwd("")
saveAsGEXF(g)