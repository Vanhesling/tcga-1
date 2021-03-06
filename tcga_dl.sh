#!/bin/sh

# Expression link
# http://gdac.broadinstitute.org/runs/stddata__2015_11_01/data/BRCA/20151101/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0.tar.gz

# Clinical phenotypes link
# http://gdac.broadinstitute.org/runs/stddata__2015_11_01/data/BRCA/20151101/gdac.broadinstitute.org_BRCA.Clinical_Pick_Tier1.Level_4.2015110100.0.0.tar.gz

diseasetypes='ACC BLCA BRCA CESC CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM'

for i in $(echo $diseasetypes); do
	explink=http://gdac.broadinstitute.org/runs/stddata__2015_11_01/data/$i/20151101/gdac.broadinstitute.org_$i.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0.tar.gz;
	outf=$i.tar.gz;
	curl -L $explink -o ~/tmp/$i.tar.gz;
	mkdir -p tcga/expression/$i;
	tar -xvf ~/tmp/$i.tar.gz -C tcga/expression/$i/ --strip-components=1;
	rm ~/tmp/*.tar.gz;

	phenlink=http://gdac.broadinstitute.org/runs/stddata__2015_11_01/data/$i/20151101/gdac.broadinstitute.org_$i.Clinical_Pick_Tier1.Level_4.2015110100.0.0.tar.gz
	curl -L $phenlink -o ~/tmp/$i.tar.gz;
	mkdir -p tcga/phen/$i;
	tar -xvf ~/tmp/$i.tar.gz -C tcga/phen/$i/ --strip-components=1;
done;
