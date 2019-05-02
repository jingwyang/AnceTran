## ---- eval=FALSE---------------------------------------------------------
#  install.packages('devtools')
#  devtools::install_github("jingwyang/AnceTran")

## ---- eval=FALSE---------------------------------------------------------
#  library('AnceTran')

## ---- message=FALSE, warning=FALSE---------------------------------------
BindingScore.table =read.table(system.file('extdata','HNF4A_meanIntensity_4Mouse.txt',package = 'AnceTran'), header = T)

head(BindingScore.table[,1:5])

## ---- message=FALSE, warning=FALSE---------------------------------------
library('AnceTran')
taxa.objects = tTFConstruct(BSFile=system.file('extdata','HNF4A_meanIntensity_4Mouse.txt',package = 'AnceTran'), taxa="all", tf="all",verbose = TRUE)

## ----message=FALSE, warning=FALSE----------------------------------------
data(TF.objects)

## ----message=FALSE, warning=FALSE----------------------------------------
library('limma')
TF_table = TFtab(objects = TF.objects, taxa = "all", tf = "all",rowindex = NULL, filtering = FALSE, normalize = FALSE, logrithm = FALSE)  
keep<-rowSums((TF_table == 0)) < ncol(TF_table)
TF_table<-TF_table[keep,]
TF_table<-data.frame(log2(normalizeQuantiles(TF_table[,])+1))

## ----message=FALSE, warning=FALSE----------------------------------------
library('ape')
dismat <- TFdist.sou(bsMat = TF_table)
colnames(dismat)=colnames(TF_table)
rownames(dismat)=colnames(dismat)
dismat

## ----message=FALSE, warning=FALSE, fig.height=4, fig.width=6-------------
tf_tree <- NJ(dismat)
tf_tree <- root(tf_tree, outgroup = "CAR_HNF4A", resolve.root = T)
tf_tree <- no0br(tf_tree)

f <- function(xx) {
  
  mat <- TFdist.sou(t(xx))
  # the distance metrics here should be the same as you specified 
  # when you created the TF-binding distance matrix 
  
  colnames(mat) <- rownames(xx)
  rownames(mat) <- colnames(mat)
  
  root(NJ(mat), "CAR_HNF4A", resolve.root = T)
  
}
bs <-  boot.phylo(tf_tree, t(TF_table), f, B = 100) 
tf_tree$node.label = bs
plot(tf_tree, show.node.label = TRUE)

## ------------------------------------------------------------------------
var_mat <- varMatInv(dismat,TF_table,phy = tf_tree)

## ---- warning=FALSE, message=FALSE---------------------------------------
mup20_binding <- TF_table[which(rownames(TF_table) == "ENSMUSG00000078672"),]

## ---- warning=FALSE, message=FALSE---------------------------------------
mup20_anc <- aee(mup20_binding, tf_tree, var_mat, select = "all")

## ---- warning=FALSE, message=FALSE,fig.height=4, fig.width=6-------------
tf_tree$node.label <- sprintf("%.4f",mup20_anc$est)
tf_tree$tip.label <- paste0(tf_tree$tip.label, "  ", sprintf("%.4f", mup20_binding))
plot(tf_tree, edge.color = "grey80", edge.width = 4,show.node.label = T,align.tip.label = T,main="Ancestial HNF4A-Binding Estimation of Gene MUP20")


