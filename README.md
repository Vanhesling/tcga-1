# Sex-biased expression in normalized TCGA data

These scripts provide a cursory overview of the sex-specific expression of genes across different cancer types:

* ACC
* BLCA
* BRCA
* CESC
* CHOL
* COAD
* DLBC
* ESCA
* GBM
* HNSC
* KICH
* KIRC
* KIRP
* LAML
* LGG
* LIHC
* LUAD
* LUSC
* OV
* PAAD
* PCPG
* PRAD
* READ
* SARC
* STAD
* TGCT
* THCA
* THYM
* UCEC
* UCS
* UVM

Note: Some cancer types removed from analysis (MESO and SKCM) due to failure to download

# Workflow

Expression and phenotype files were downloaded and sorted in the working directory using the `dl_tcga.sh` script, to have the following directory structure:

	├── expression
	│   ├── ACC
	│   │   ├── ACC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── BLCA
	│   │   ├── BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── BRCA
	│   │   ├── BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── CESC
	│   │   ├── CESC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── CHOL
	│   │   ├── CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── COAD
	│   │   ├── COAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── DLBC
	│   │   ├── DLBC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── ESCA
	│   │   ├── ESCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── GBM
	│   │   ├── GBM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── HNSC
	│   │   ├── HNSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── KICH
	│   │   ├── KICH.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── KIRC
	│   │   ├── KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── KIRP
	│   │   ├── KIRP.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── LAML
	│   │   ├── LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── LGG
	│   │   ├── LGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── LIHC
	│   │   ├── LIHC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── LUAD
	│   │   ├── LUAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── LUSC
	│   │   ├── LUSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   │   └── MANIFEST.txt
	│   ├── MESO
	│   ├── OV
	│   │   ├── MANIFEST.txt
	│   │   └── OV.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── PAAD
	│   │   ├── MANIFEST.txt
	│   │   └── PAAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── PCPG
	│   │   ├── MANIFEST.txt
	│   │   └── PCPG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── PRAD
	│   │   ├── MANIFEST.txt
	│   │   └── PRAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── READ
	│   │   ├── MANIFEST.txt
	│   │   └── READ.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── SARC
	│   │   ├── MANIFEST.txt
	│   │   └── SARC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── SKCM
	│   ├── STAD
	│   │   ├── MANIFEST.txt
	│   │   └── STAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── TGCT
	│   │   ├── MANIFEST.txt
	│   │   └── TGCT.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── THCA
	│   │   ├── MANIFEST.txt
	│   │   └── THCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── THYM
	│   │   ├── MANIFEST.txt
	│   │   └── THYM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── UCEC
	│   │   ├── MANIFEST.txt
	│   │   └── UCEC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   ├── UCS
	│   │   ├── MANIFEST.txt
	│   │   └── UCS.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	│   └── UVM
	│       ├── MANIFEST.txt
	│       └── UVM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
	├── phen
	│   ├── ACC
	│   │   ├── ACC.clin.merged.picked.txt
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_ACC.tsv
	│   ├── BLCA
	│   │   ├── All_CDEs.txt
	│   │   ├── BLCA.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_BLCA.tsv
	│   ├── BRCA
	│   │   ├── All_CDEs.txt
	│   │   ├── BRCA.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_BRCA.tsv
	│   ├── CESC
	│   │   ├── All_CDEs.txt
	│   │   ├── CESC.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_CESC.tsv
	│   ├── CHOL
	│   │   ├── All_CDEs.txt
	│   │   ├── CHOL.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_CHOL.tsv
	│   ├── COAD
	│   │   ├── All_CDEs.txt
	│   │   ├── COAD.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_COAD.tsv
	│   ├── DLBC
	│   │   ├── All_CDEs.txt
	│   │   ├── DLBC.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_DLBC.tsv
	│   ├── ESCA
	│   │   ├── All_CDEs.txt
	│   │   ├── ESCA.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_ESCA.tsv
	│   ├── GBM
	│   │   ├── All_CDEs.txt
	│   │   ├── GBM.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_GBM.tsv
	│   ├── HNSC
	│   │   ├── All_CDEs.txt
	│   │   ├── HNSC.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_HNSC.tsv
	│   ├── KICH
	│   │   ├── All_CDEs.txt
	│   │   ├── KICH.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_KICH.tsv
	│   ├── KIRC
	│   │   ├── All_CDEs.txt
	│   │   ├── KIRC.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_KIRC.tsv
	│   ├── KIRP
	│   │   ├── All_CDEs.txt
	│   │   ├── KIRP.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_KIRP.tsv
	│   ├── LAML
	│   │   ├── All_CDEs.txt
	│   │   ├── LAML.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_LAML.tsv
	│   ├── LGG
	│   │   ├── All_CDEs.txt
	│   │   ├── LGG.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_LGG.tsv
	│   ├── LIHC
	│   │   ├── All_CDEs.txt
	│   │   ├── LIHC.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_LIHC.tsv
	│   ├── LUAD
	│   │   ├── All_CDEs.txt
	│   │   ├── LUAD.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_LUAD.tsv
	│   ├── LUSC
	│   │   ├── All_CDEs.txt
	│   │   ├── LUSC.clin.merged.picked.txt
	│   │   ├── MANIFEST.txt
	│   │   └── stage3_params_clin_selection_LUSC.tsv
	│   ├── MESO
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── MESO.clin.merged.picked.txt
	│   │   └── stage3_params_clin_selection_MESO.tsv
	│   ├── OV
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── OV.clin.merged.picked.txt
	│   │   └── stage3_params_clin_selection_OV.tsv
	│   ├── PAAD
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── PAAD.clin.merged.picked.txt
	│   │   └── stage3_params_clin_selection_PAAD.tsv
	│   ├── PCPG
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── PCPG.clin.merged.picked.txt
	│   │   └── stage3_params_clin_selection_PCPG.tsv
	│   ├── PRAD
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── PRAD.clin.merged.picked.txt
	│   │   └── stage3_params_clin_selection_PRAD.tsv
	│   ├── READ
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── READ.clin.merged.picked.txt
	│   │   └── stage3_params_clin_selection_READ.tsv
	│   ├── SARC
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── SARC.clin.merged.picked.txt
	│   │   └── stage3_params_clin_selection_SARC.tsv
	│   ├── SKCM
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── SKCM.clin.merged.picked.txt
	│   │   └── stage3_params_clin_selection_SKCM.tsv
	│   ├── STAD
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── STAD.clin.merged.picked.txt
	│   │   └── stage3_params_clin_selection_STAD.tsv
	│   ├── TGCT
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── stage3_params_clin_selection_TGCT.tsv
	│   │   └── TGCT.clin.merged.picked.txt
	│   ├── THCA
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── stage3_params_clin_selection_THCA.tsv
	│   │   └── THCA.clin.merged.picked.txt
	│   ├── THYM
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── stage3_params_clin_selection_THYM.tsv
	│   │   └── THYM.clin.merged.picked.txt
	│   ├── UCEC
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── stage3_params_clin_selection_UCEC.tsv
	│   │   └── UCEC.clin.merged.picked.txt
	│   ├── UCS
	│   │   ├── All_CDEs.txt
	│   │   ├── MANIFEST.txt
	│   │   ├── stage3_params_clin_selection_UCS.tsv
	│   │   └── UCS.clin.merged.picked.txt
	│   └── UVM
	│       ├── All_CDEs.txt
	│       ├── MANIFEST.txt
	│       ├── stage3_params_clin_selection_UVM.tsv
	│       └── UVM.clin.merged.picked.txt
	├── README.md
	├── results.Robj
	├── Rplots.pdf
	├── sub_tcga.sh
	├── tcga_dl.sh
	└── tcga.R

The expression txt files contain two headers, followed by normalized RPKM values for genes (rows) by individuals (columns)
The clinical txt files contain various phenotype info (rows) by individuals (columns). The first row corresponds to truncated individual ids of the expression data.

The `tcga.R` script processes the expression and phenotype files, convertng the expression individual ids
