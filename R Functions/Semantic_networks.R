library(tm)
library(igraph)
library(RColorBrewer)
library(wordcloud)
library(SnowballC)
library(ggplot2)


setwd("")

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

graph_to.gephi<-function(s){
  g<-graph.data.frame(s,directed = F)
  saveAsGEXF(g)
}
