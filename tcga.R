performDE <- function(expr, phen) {
	log2expr <- log2(expr+0.5)	
	lowexpr <- (rowMeans(log2expr) < 1)
	log2expr <- log2expr[!lowexpr,]
	sex <- as.numeric(phens[[1]][,'gender']=='male')
	design <- model.matrix(~sex)
	y <- DGEList(expr, remove.zeros=T)
	#y <- calcNormFactors(y)
	v <- voom(y, design)
	#et <- exactTest(y, dispersion=0.4^2)
	fit <- lmFit(log2expr, design)
	fit <- eBayes(fit)
	return(fit)
}

library(edgeR)
library(limma)

exp_dir <- "/home/t.cri.cczysz/tcga/expression"
phen_dir <- "/home/t.cri.cczysz/tcga/phen"

cancer <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

exprs <- list()
phens <- list()

for (i in cancer[1]) {
	files <- list.files(paste(exp_dir, i, sep='/'), pattern='*.txt')
	x <- read.table(paste(exp_dir,i,files[1],sep='/'), header=T, row.names=1, skip=2, sep='\t')
	y <- read.table(paste(exp_dir,i,files[1],sep='/'), header=T, row.names=1, skip=0, sep='\t')
	colnames(x) <- colnames(y)
	exprs[[i]] <- x
	
	files <- list.files(paste(phen_dir, i, sep='/'), pattern='*.txt')
	phen <- read.table(paste(phen_dir,i,files[1],sep='/'), header=F, row.names=1, skip=0, sep="\t")

	#exp_ids <- tolower(apply(matrix(colnames(exprs[[1]]), ncol=1), 1, function(x) {unlist(strsplit(x, split='[.]'))[1:3]})[3,])
	exp_ids <- tolower(apply(matrix(colnames(y), ncol=1), 1, function(x) {unlist(strsplit(x, split='[.]'))[1:3]})[3,])
	phen_ids <- apply(matrix(t(phen)[,1], ncol=1), 1, function(x) {unlist(strsplit(x, split='[-]'))[3]})

	phens[[i]] <- t(phen)[phen_ids%in%exp_ids,]
}

et <- performDE(exprs[[1]], phens[[1]])
