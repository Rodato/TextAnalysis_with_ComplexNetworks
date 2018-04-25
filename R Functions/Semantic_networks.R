library(tm)
library(igraph)
library(rgexf)   
library(RColorBrewer)
library(ggplot2)
library(dplyr)

###The function saveasGEFX was not build by me. I get it in internet (don't remember where, maybe 
### in github). saveAsGEXF exports an igraph object and write it in GEXF format
### to deal with the visualization in gephi.

###
get_semantic_network<-function(paper,fqw){

doc<-readPDF(control=list(text="-layout"))(elem=list(uri=paper), language="en")

toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))

text_raw <- doc$content
text_corpus <- Corpus(VectorSource(text_raw))

text_docs <- tm_map(text_corpus, toSpace, "/")
text_docs <- tm_map(text_docs, toSpace, "@")
text_docs <- tm_map(text_docs, toSpace, "\\|")

text_docs<-tm_map(text_docs, content_transformer(function(x) iconv(enc2utf8(x), sub = "byte")))

# Convert text to lower case
text_docs <- tm_map(text_docs, content_transformer(tolower))
# Remove numbers
text_docs <- tm_map(text_docs, removeNumbers)
# Remove english common stopwords
text_docs <- tm_map(text_docs, removeWords, stopwords("english"))
# Remove your own stop word
# specify your stopwords as a character vector
text_docs <- tm_map(text_docs, removeWords, c("journal", "vol", "may", "one",
                                              "effect","two","also","http","jstororg",
                                              "utc","downloaded","however"))
# Remove punctuations
text_docs <- tm_map(text_docs, removePunctuation)
# Eliminate extra white spaces
text_docs <- tm_map(text_docs, stripWhitespace)

dtm <- TermDocumentMatrix(text_docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)

#Constructing the semantic network for the ten most frequent words

semantic_network<-data.frame()
dw<-as.data.frame(d$word)

for(word in dw[fqw,]){
  x<-as.data.frame(findAssocs(dtm, terms=word, corlimit = 0.3))
  x$c<-rownames(x)
  x$pw<-word
  colnames(x)<-c("cor","c2","c1")
  x<-x[,c(3,2,1)]
  semantic_network<-rbind(semantic_network,x)
}
return(semantic_network)
}


###
graph_to.gephi<-function(s,e){
  saveAsGEXF<-function(g, filepath="converted_graph.gexf"){
    require(igraph)
    require(rgexf)
    
    # gexf nodes require two column data frame (id, label)
    # check if the input vertices has label already present
    # if not, just have the ids themselves as the label
    if(is.null(V(g)$label))
      V(g)$label <- as.character(V(g))
    
    # similarily if edges does not have weight, add default 1 weight
    if(is.null(E(g)$weight))
      E(g)$weight <- rep.int(1, ecount(g))
    
    nodes <- data.frame(cbind(V(g), V(g)$label))
    edges <- t(Vectorize(get.edge, vectorize.args='id')(g, 1:ecount(g)))
    
    # combine all node attributes into a matrix (and take care of & for xml)
    vAttrNames <- setdiff(list.vertex.attributes(g), "label") 
    nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&",get.vertex.attribute(g, attr))))
    
    # combine all edge attributes into a matrix (and take care of & for xml)
    eAttrNames <- setdiff(list.edge.attributes(g), "weight") 
    edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&",get.edge.attribute(g, attr))))
    
    # combine all graph attributes into a meta-data
    graphAtt <- sapply(list.graph.attributes(g), function(attr) sub("&", "&",get.graph.attribute(g, attr)))
    
    # generate the gexf object
    output <- write.gexf(nodes, edges, 
                         edgesWeight=E(g)$weight,
                         edgesAtt = edgesAtt,
                         nodesAtt = nodesAtt,
                         meta=c(list(creator="Gopalakrishna Palem", description="igraph -> gexf converted file", keywords="igraph, gexf, R, rgexf"), graphAtt))
    
    print(output, filepath, replace=T)
  }
  g<-graph.data.frame(s,directed = F,vertices = e)
  saveAsGEXF(g)
}

joint_networks<-function(df1,df2){
  colnames(df1)<-c("w1.1","w2.1","cor.1")
  colnames(df2)<-c("w1.2","w2.2","cor.2")
  cw<-inner_join(df1,df2,by=c("w1.1"="w1.2"))
  cwu1<-as.data.frame(unique(cw[,1]))
  cwu2<-as.data.frame(unique(cw[,2]))
  cwu3<-as.data.frame(unique(cw[,4]))
  colnames(cwu1)<-c("c_edges")
  colnames(cwu2)<-c("c_edges")
  colnames(cwu3)<-c("c_edges")
  cw<-rbind(cwu1,cwu2,cwu3)
  cw<-as.data.frame(unique(cw[,1]))
  cw$mark<-"CE"
  colnames(cw)<-c("c_edges","CE")
  return(cw)
}

all_verts<-function(df1,df2){
  colnames(df1)<-c("w1","w2","cor")
  colnames(df2)<-c("w1","w2","cor")
  df1$mark<-"Cook"
  df2$mark<-"Angrist"
  df<-rbind(df1,df2)
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

total_nodes<-function(df1,df2){
  colnames(df1)<-c("c_edges","mark")
  colnames(df2)<-c("c_edges","mark")
  fj<-df2[!df2$c_edges%in%df1$c_edges,]
  f<-rbind(df1,fj)
  f<-f[!duplicated(f[,c(1)]), ]
  return(f)
}


##Order to runing the functions
s1<-get_semantic_network("doc1.pdf",1:15)
s2<-get_semantic_network("doc.pdf",1:15)
j1<-joint_networks(s1,s2)
f1<-all_verts(s1,s2,"doc1","doc2")
f2<-total_nodes(j1,f1)
s<-rbind(s1,s2)
g<-graph.data.frame(s,directed = F,vertices = f2)
g<-simplify(g,remove.multiple = T,remove.loops = F)
saveAsGEXF(g)