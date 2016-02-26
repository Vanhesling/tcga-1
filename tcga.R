performDE <- function(expr, phen) {
	log2expr <- log2(expr+0.5)	
	lowexpr <- (rowMeans(log2expr) < 1)
	log2expr <- log2expr[!lowexpr,]
	sex <- as.numeric(phen[,'gender']=='male')
	print(as.character(sex))
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

cancer <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

exprs <- list()
phens <- list()

for (i in cancer) {
	files <- list.files(paste(exp_dir, i, sep='/'), pattern='*.data.txt')
	x <- read.table(paste(exp_dir,i,files[1],sep='/'), header=T, row.names=1, skip=2, sep='\t')
	y <- read.table(paste(exp_dir,i,files[1],sep='/'), header=T, row.names=1, skip=0, sep='\t')
	colnames(x) <- colnames(y)
	rm(y)
	
	files <- list.files(paste(phen_dir, i, sep='/'), pattern='*.picked.txt')
	phen <- read.table(paste(phen_dir,i,files[1],sep='/'), header=F, row.names=1, skip=0, sep="\t")

	#exp_ids <- tolower(apply(matrix(colnames(exprs[[1]]), ncol=1), 1, function(x) {unlist(strsplit(x, split='[.]'))[1:3]})[3,])
	exp_ids <- tolower(apply(matrix(colnames(x), ncol=1), 1, function(x) {paste(unlist(strsplit(x, split='[.]'))[1:3], collapse='-')}))
	#phen_ids <- apply(matrix(t(phen)[,1], ncol=1), 1, function(x) {unlist(strsplit(x, split='[-]'))[3]})
	phen_ids <- as.character(t(phen)[,1])

	# Remove entries with duplicates
	colnames(x) <- exp_ids
	unique_exp <- x[, !(exp_ids%in%exp_ids[duplicated(exp_ids)])]

	phens[[i]] <- t(phen)[phen[1,]%in%colnames(unique_exp),]
	exprs[[i]] <- unique_exp
	if (F) {
	if (length(phen_ids) >= length(exp_ids)) {
		phens[[i]] <- t(phen)[phen[1,]%in%exp_ids,]
	} else {
		phens[[i]] <- t(phen)[exp_ids%in%phen_ids,] 
		exprs[[i]] <- x[exp_ids%in%phen_ids] 
	}
	}

}

et_list <- list()
for (i in cancer) {
	et <- performDE(exprs[[i]], phens[[i]])
	et_list[[i]] <- et
}

results <- lapply(et_list, topTable, number=Inf)
save(results, file='/home/t.cri.cczysz/tcga/results.Robj')
