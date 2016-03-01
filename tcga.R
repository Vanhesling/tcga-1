importData <- function(cancer, exp_dir = exp_dir, phen_dir = phen_dir) {

	i <- cancer
	files <- list.files(paste(exp_dir, i, sep='/'), pattern='*.data.txt')

	# Ignore first two header rows and import expression
	x <- read.table(paste(exp_dir, i, files[1],sep='/'), header=T, row.names=1, skip=2, sep='\t')
	# Import header only
	y <- read.table(paste(exp_dir, i, files[1],sep='/'), header=T, row.names=1, skip=0, sep='\t', nrows=2)
	colnames(x) <- colnames(y)
	rm(y)
	
	files <- list.files(paste(phen_dir, i, sep='/'), pattern='*.picked.txt')
	phen <- read.table(paste(phen_dir,i,files[1],sep='/'), header=F, row.names=1, skip=0, sep="\t")

	exp_ids <- tolower(apply(matrix(colnames(x), ncol=1), 1, function(x) {paste(unlist(strsplit(x, split='[.]'))[1:3], collapse='-')}))
	phen_ids <- as.character(t(phen)[,1])

	# Remove entries with duplicates
	colnames(x) <- exp_ids
	unique_exp <- x[, !(exp_ids%in%exp_ids[duplicated(exp_ids)])]

	phens <- t(phen)[phen[1,]%in%colnames(unique_exp),]
	exprs <- unique_exp[, colnames(unique_exp)%in%phens[[i]][,1]]
	return(list(phens, exprs))

}

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
library(biomaRt)
library(ggplot2)
library(goseq)

exp_dir <- "/home/t.cri.cczysz/tcga/expression"
phen_dir <- "/home/t.cri.cczysz/tcga/phen"


cancer <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

full_names <- c('Adrenocortical Carcinoma', 
	'Breast Lobular Carcinoma',
	'Breast Invasive Carcinoma',
	'Cervical Squamous Cell Carcinoma',
	'Cholangiocarcinoma',
	'Colorectal Adenocarcinoma', 
	'Diffuse Large B-cell',
	'Esophageal Carcinoma',
	'Glioblastoma Multiforme',
	'Head/Neck Squamous Cell Carcinoma',
	'Kidney Chromophobe',
	'Kidney Reneal Clear Cell',
	'Kidney Renal Papillary Cell Carcinoma',
	'Acute Myeloid Lukemia',
	'Lower Grade Glioma',
	'Liver Hepatocellular Carcinoma',
	'Lung Adenocarcinoma',
	'Lung Squaous Cell Carcinoma',
	'Ovarian',
	'Pancreatic Adenocarcinoma',
	'Pheochromocytoma and Paraganglioma',
	'Prostate Adenocarcinoma',
	'Rectum Adenocarcinoma',
	'Sarcoma',
	'Stomach Adenocarcinoma',
	'Testicular Germ Cell',
	'Thyroid Cancer',
	'Thymoma',
	'Uterine Corpus Endometrial Carcinoma',
	'Uterine Carcinosarcoma',
	'Uveal Melanoma'
)
if (!file.exists('/home/t.cri.cczysz/tcga/results.Robj')) {
	exprs <- list()
	phens <- list()

	for (i in cancer) {
		out <- importData(i)
		exprs[[i]] <- out[[2]]
		phens[[i]] <- out[[1]]
	}

	et_list <- list()
	for (i in cancer) {
		et <- performDE(exprs[[i]], phens[[i]])
		et_list[[i]] <- et
	}

	#results <- lapply(et_list, topTable, number=Inf)
	save(et_list, file='/home/t.cri.cczysz/tcga/results.Robj') 
} else {load('/home/t.cri.cczysz/tcga/results.Robj')}

# et_list[[cancer]] gives lmFit objects from lm of gene ~ sex test

# Perform GO enrichment analysis
go.data <- list()
if (!file.exists('go.Robj')) {
	for (i in names(et_list)) {
		tmp <- et_list[[i]]
		ttable <- topTable(tmp, number=Inf)
		sig <- topTable(tmp, number=Inf, p.value=0.05)
		all_genes <- rownames(topTable(tmp, number=Inf))
		go_in <- as.numeric(all_genes%in%rownames(sig))
		all_genes <- apply(matrix(all_genes, ncol=1), 1, function (x) {unlist(strsplit(x, '[|]'))[2]})
		names(go_in) <- all_genes
		pwf <- nullp(go_in, "hg19", "refGene")
		GO.wall <- goseq(pwf, "hg19", "refGene", method='Sampling', repcnt=10000)
		GO.wall$q.value <- p.adjust(GO.wall$over_represented_pvalue,method="fdr")
		go.data[[i]] <- GO.wall
		#GO.enriched <- subset(GO.wall, q.value < 0.1)
	}
	save(go.data, file='go.Robj') 
	} else {load(file='go.Robj')}
q()

ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
meta_list <- c('cancer', 'name', 'sample_size', 'nmale', 'nfemale', 'ngenes', 'nsiggenes', 'malebiased', 'percentmale', 'femalebiased', 'percentfemale')
if (!file.exists('/home/t.cri.cczysz/tcga/out.Robj')) {
	out_list <- list()
	for (i in seq(length(et_list))) {
		cancer_type <- names(et_list)[i]
		nmale <- sum(et_list[[i]]$design[,2])
		nfemale <- sum(!et_list[[i]]$design[,2])

		all_genes <- topTable(et_list[[i]], number=Inf)
		sig_genes <- topTable(et_list[[i]], number=Inf, p.value=0.05)

		ngenes <- nrow(all_genes)	
		nsig <- nrow(sig_genes)
		male_bias <- sum(sig_genes$logFC > 0)
		female_bias <- sum(sig_genes$logFC < 0)

		meta_list <- rbind(meta_list, c(cancer_type, full_names[i], nmale+nfemale, nmale, nfemale, ngenes, nsig, male_bias, signif(100 * male_bias/nsig, 2), female_bias, signif(100*female_bias/nsig, 2)))

		sig_symbols <- apply(matrix(rownames(sig_genes), ncol=1), 1, function(x) {unlist(strsplit(x, split='[|]'))[1]})
		sig_entrez <- apply(matrix(rownames(sig_genes), ncol=1), 1, function(x) {unlist(strsplit(x, split='[|]'))[2]})

		if (length(sig_symbols)==0) {
			#out_list[[cancer_type]] <- NA
			next
		}

		biomart.results=getBM(ensembl, attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene", "chromosome_name"),filters="entrezgene",values=sig_entrez)

		# Reorder biomart results by entrez id
		x <- biomart.results[match(sig_entrez, biomart.results$entrezgene),]
		out <- cbind(x, logFC = signif(sig_genes$logFC,2), p.val = signif(sig_genes$P.Value,2), q.value = signif(sig_genes$adj.P.Val,2))
		rownames(out) <- rownames(sig_genes)
		out_list[[cancer_type]] <- out
		save(out_list, file='/home/t.cri.cczysz/tcga/out.Robj')
	}
} else {load(file='/home/t.cri.cczysz/tcga/out.Robj')}
	
write.table(meta_list, file='meta.csv', sep=',', row.names=F, col.names=F, quote=F)

plot.data <- data.frame(meta_list[-1,-1])
colnames(plot.data) <- meta_list[1,][-1]
rownames(plot.data) <- meta_list[,1][-1]

pdf(file='sig_genes.pdf', width=11, height=8)
	g <- ggplot(data=plot.data, aes(x=rownames(plot.data), y=as.numeric(nsiggenes)))
	g + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	g <- ggplot(data=plot.data, aes(x=rownames(plot.data), y=log10(as.numeric(nsiggenes))))
	g + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

if (F) {
for (cancer_type in names(out_list)) {
	f.out <- paste(cancer_type, 'results.csv', sep='.')
	write.csv(out_list[[cancer_type]], file=f.out, sep=',', row.names=T, col.names=T, quote=F)
}
}
