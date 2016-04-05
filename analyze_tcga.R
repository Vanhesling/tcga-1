library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")

getBioMart <- function(ensembl, genes) {
	sig_entrez <- apply(matrix(genes, ncol=1), 1, function(x) {unlist(strsplit(x, split='[|]'))[2]})
	biomart.results=getBM(ensembl, attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene", "chromosome_name"),filters="entrezgene",values=sig_entrez)

	# Reorder biomart results by entrez id
	x <- biomart.results[match(sig_entrez, biomart.results$entrezgene),]
	return(x)
}

setwd('/home/t.cri.cczysz/tcga')

load('de_results_sva.Robj')

dat <- list()
for (i in seq(length(et_list))) {
        coefs <- et_list[[i]]$coefficients
        pvals <- et_list[[i]]$p.value
        if (ncol(coefs)==1) {
                q.vals <- p.adjust(pvals, method='fdr')
                df <- data.frame(coefs = coefs, pvals = pvals, qvals = q.vals)
		dat[[names(et_list)[i]]] <- df
        }
}

sig_list <- list()

for (i in seq(length(dat))) {
	cancer_type <- names(dat)[i]
	sig_genes <- subset(dat[[i]], (qvals <= 0.05 & (coefs < -0.5 | coefs > 0.5)))
	sig_list[[cancer_type]] <- sig_genes
	if (nrow(sig_genes) > 0) {
		bio_out <- getBioMart(ensembl, rownames(sig_genes))
		f.out <- paste(cancer_type, 'sig_results.csv', sep='_')
		write.table(cbind(rownames(sig_genes), sig_genes, bio_out), file=f.out, quote=F, row.names=F, col.names=F, sep=',')
		write.table(na.omit(bio_out$entrezgene), file=paste('go_data', f.out, sep='/'), row.names=F, col.names=F, quote=F)
		
		bio_out <- getBioMart(ensembl, rownames(dat[[i]]))
		f.out <- paste(cancer_type, 'all.csv', sep='_')
		write.table(na.omit(bio_out$entrezgene), file=paste('go_data', f.out, sep='/'), row.names=F, col.names=F, quote=F)
	}
}
