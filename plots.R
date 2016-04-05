setwd('/home/t.cri.cczysz/tcga')

load('de_results_sva.Robj') # et_list

pdf('foldChange.pdf', width=8, height=8)
for (i in seq(length(et_list))) {
	coefs <- et_list[[i]]$coefficients
	if (ncol(coefs)==1) {
		plot(density(coefs), main=names(et_list)[i])
	}
}
dev.off()

pdf('volcanoPlots.pdf', width=8, height=8)
for (i in seq(length(et_list))) {
	coefs <- et_list[[i]]$coefficients
	pvals <- et_list[[i]]$p.value
	if (ncol(coefs)==1) {
		q.vals <- p.adjust(pvals, method='fdr')
		df <- data.frame(coefs = coefs, pvals = pvals, qvals = q.vals)
		try(sig_p <- -log10(df[order(df$qvals),][max(which(df[order(df$qvals),'qvals']<=0.05)),'pvals']))
		plot(coefs, -log10(pvals), main=names(et_list)[i])
		if (!is.na(sig_p)) { abline(sig_p, 0) }
	}
}
dev.off()

