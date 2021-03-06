setwd("~/PhD/fat-NAT")
set.seed(123)

#############################################################
# Load libraries ############################################

# Libraries
library(tidyverse)
library(magrittr)
library(foreach)
source("~/third_party_sofware/utilities_normalization.R")
library(sva)
library(edgeR)

# RUV
library(ruv)
source("~/third_party_sofware//RUV4.R")

# Heatmap
#library(ggdendro)

library(RColorBrewer)

# as_matrix function
source("https://gist.githubusercontent.com/stemangiola/39b08d529157ce59a5ff5dc1653951c5/raw/5ceb71b95d33254016b15fdc28c3ec40d7dbe137/as_matrix.R")

# EGSEA
library(EGSEA)
library(EGSEAdata)

# Safe detach
# detach("package:AnnotationDbi", unload=TRUE, force = T)
# detach("package:hgu95av2.db", unload=TRUE, force = T)

# Deconvolution
library(devtools)
#install_github("stemangiola/ARMET", args = "--preclean", build_vignettes = FALSE, auth_token = "37c5c6238136a6804d336d9a7078eece993ce870", password="x-oauth-basic")  
library(ARMET)

# Cibersort
source("~/PhD/deconvolution/ARMET_BK_Apr2017/comparison_methods/cibersort/CIBERSORT_annotated.R")


# Set theme plots
my_theme = 	
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=8),
		legend.position="bottom",
		aspect.ratio=1,
		axis.text.x = element_text(angle = 90, hjust = 1),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

#############################################################
# Build data set ############################################

if(0){
	bam.files <- dir("alignment_hg19", pattern=".sam", full.names=T, recursive=T)
	bam.counts <- featureCounts(bam.files,countMultiMappingReads=TRUE,fraction=TRUE, annot.inbuilt="hg19", nthreads=40)
	d <- DGEList(counts=bam.counts$counts, genes=bam.counts$annotation[,c("GeneID","Length")])
}
load("counts_paired_end.RData")

#library(hgu95av2.db)

d$genes$symbol <- 
	AnnotationDbi:::mapIds(
		org.Hs.eg.db::org.Hs.eg.db,
		keys=as.character(d$genes$GeneID),
		column="SYMBOL",
		keytype="ENTREZID",
		multiVals="first"
	)
colnames(d$counts) = 
	sapply(
		colnames(d$counts), 
		function(cc) strsplit(cc, ".", fixed=T)[[1]][2]
	) 

# Build annotation df
annot = 
	read_csv("info.csv") %>%
	filter(Label %in% c("Neoadjuvant", "High")) %>%
	mutate(Recurrence =
				 	ifelse(
				 		Sample %in% c("11104PP","11204PP","11218PP","11086PP"),
				 		1,
				 		ifelse(
				 			Sample %in% c("11182PP","11184PP","11160PP","11165PP"),
				 			0,
				 			NA
				 		)
				 	)
	) %>%
	mutate_if(is.character, as.factor) %>%
	arrange(Sample)

d = d[,annot %>% pull(Sample) %>% as.character()]

# print sequencing output
writeLines("Sequencing output")
d$counts %>% colSums() %>% summary() 

#############################################################
# Normalise #################################################

# Function that takes annotated tibble and save MDS plot for the first 8 PC
tbl_to_MDS_plot = function(my_df_mds, annot, read_count = "Read count", file_name, limits = c(-1.5, 1.5)){
	
	foreach(
		components = list(c(1, 2), c(3, 4), c(5, 6), c(7, 8)), 
		#my_df_mds = (.),
		.combine = bind_rows
	) %do% {
		my_df_mds %>%
			dplyr::select(symbol, Sample, !!read_count) %>%
			spread( Sample, !!read_count) %>%
			dplyr::select(-symbol) %>%
			as.matrix() %>%
			plotMDS(dim.plot = components) %>% 
			{
				tibble(Sample = names((.)$x), x = (.)$x, y = (.)$y, PCx = components[1], PCy = components[2])
			}
	} %>%
		
		# Annotate
		left_join(annot) %>%
		mutate(Sample = gsub("PP", "", as.character(Sample))) %>%
		mutate_if(is.character, as.factor) %>%
		
		gather(Annotation, Value, c("Batch", "Label", "kit", "Recurrence")) %>%
		{
			
			ggplot(data=(.), aes(x = x, y = y, label = Sample, PCx = PCx, PCy = PCy)) + 
				geom_point(aes(fill = Value), size=3, shape=21, color="grey20") +
				ggrepel::geom_text_repel(
					size = 1, 
					point.padding = 0.3, 
					#fontface = 'bold', 
					# label.padding = 0.1, 
					# label.size = 0,
					segment.size = 0.2,
					seed = 123
				) +
				xlim(limits) +
				ylim(limits) +
				#geom_text(color = "grey20", size = 2 ) +
				scale_fill_brewer(palette = "Set1") +
				facet_grid(sprintf("PC %s", interaction(PCx, PCy))~Annotation) +
				theme_bw() +
				theme(
					panel.border = element_blank(), 
					axis.line = element_line(),
					panel.grid.major = element_line(size = 0.2),
					panel.grid.minor = element_line(size = 0.1),
					text = element_text(size=12),
					legend.position="bottom",
					aspect.ratio=1,
					strip.background = element_blank(),
					axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
					axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
				) +
				xlab("Principal component x") + ylab("Principal component y") 
		} %>%
		ggsave(plot = .,
					 file_name,
					 useDingbats=FALSE,
					 units = c("mm"),
					 width = 183 ,
					 height = 183 
		)
}


# collect expression data
d_adj =
	d$counts %>% 
	as_tibble %>%
	mutate(symbol = d$genes$symbol) %>%
	gather(Sample, value, -symbol) %>%
	drop_na() %>%
	
	# Normalise
	norm_RNAseq( 
		sample_column = "Sample", 
		gene_column = "symbol", 
		value_column = "value"
	) %>%
	filter(!filt_for_calc) %>%
	dplyr::select(symbol, Sample, `value normalised`, value) %>%
	mutate(`value normalised log` = log(`value normalised` + 1)) %>%
	
	# Plot densities
	{
		
		getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
		
		(.) %>%
			gather(is_normalised, value, c("value", "value normalised")) %>%
			{
				ggplot((.), aes(value + 1, color=Sample)) +
					geom_density() +
					facet_grid(~ is_normalised) +
					scale_color_manual(values = getPalette( (.) %>% distinct(Sample) %>% nrow )) +
					#scale_color_brewer(palette = "Set1") +
					scale_x_log10() +
					theme_bw() +
					theme(
						panel.border = element_blank(), 
						axis.line = element_line(),
						panel.grid.major = element_line(size = 0.2),
						panel.grid.minor = element_line(size = 0.1),
						text = element_text(size=12),
						legend.position="bottom",
						aspect.ratio=1,
						strip.background = element_blank(),
						axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
						axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
					)
			} %>%
			ggsave(plot = .,
						 "out_treatment_vs_high/nat_vs_high_density.pdf",
						 useDingbats=FALSE,
						 units = c("mm"),
						 width = 183 ,
						 height = 183 /2 + 30 
			)
		
		(.)
	} %>%
	
	# Adjustment + MDS
	{
		
		my_df_mds = (.)
		
		# 1. Landscape
		my_df_mds %>%
			tbl_to_MDS_plot(
				annot, 
				read_count = "value normalised log",
				file_name = "out_treatment_vs_high/mds_plot_landscape.pdf"
			)
		
		# 2. MDS plot after batch correction
		hk_600 = 
			as.character(read.table("~/PhD/deconvolution/hk_600.txt", quote="\"", comment.char="")[,1]) %>% 
			intersect(my_df_mds %>% pull(symbol))
		
		my_df_mds.RUV = 
			my_df_mds %>%
			dplyr::select(symbol, Sample, `value normalised log`) %>%
			spread( Sample, `value normalised log`) %>%
			arrange(desc(symbol %in% hk_600)) %>%
			as.data.frame() %>%
			magrittr::set_rownames(.$symbol) %>%
			dplyr::select(-symbol) %>%
			as.matrix() %>%
			
			# Execute RUV
			{
				Y = (.)
				
				RUV4(
					Y %>% t(),
					annot %>% 
						arrange(match(Sample, colnames(Y))) %>% 
						mutate(is_neoadjuvant = Label == "Neoadjuvant") %>% 
						pull(is_neoadjuvant) %>% 
						as.numeric() %>% 
						as.matrix(),
					rownames(Y)%in%hk_600, 
					k=1, 
					Z = 1, 
					hkg = hk_600
				) %>%
				{
					RUVres = (.)
					
					# Plot unwanted covariates
					
					# Recalculated corrected expresison
					matrix(
						pmax(
							Y %>% t() - 
								RUVres$W %*% RUVres$alpha_all,
							0
						),
						nr=nrow(Y %>% t()),
						ncol=ncol(Y %>% t())
					) %>% 
						
						magrittr::set_colnames(Y %>% rownames()) %>%
						magrittr::set_rownames(Y %>% colnames()) %>%
						as_tibble(rownames="Sample") %>%
						bind_cols(W = RUVres$W) %>%
						gather(symbol,`value RUV log`, -Sample, -W ) 
				} 
			} 
		
		# Produce MDS
		my_df_mds.RUV %>%
			dplyr::select(Sample, symbol, `value RUV log`) %>%
			tbl_to_MDS_plot(
				annot, 
				read_count = "value RUV log",
				file_name = "mds_plot_landscape_RUV.pdf"
			)
		
		# 3. Combat correction of KIT
		my_df_mds.Combat = 
			my_df_mds %>%
			dplyr::select(symbol, Sample, `value normalised log`) %>%
			spread( Sample, `value normalised log`) %>%
			arrange(desc(symbol %in% hk_600)) %>%
			as.data.frame() %>%
			magrittr::set_rownames(.$symbol) %>%
			dplyr::select(-symbol) %>%
			as.matrix() %>%
			ComBat(
				batch=annot %>% arrange(match(Sample, colnames((.)))) %>% pull(kit), 
				mod=model.matrix(	~ annot %>% arrange(match(Sample, colnames((.)))) %>% pull(Label)	)
			)  %>%
			as_tibble(rownames="symbol") %>%
			gather(Sample,`value Combat log`, -symbol ) 
		#left_join(my_df_mds) %>%
		
		# Produce MDS
		my_df_mds.Combat %>%
			tbl_to_MDS_plot(
				annot, 
				read_count = "value Combat log",
				file_name = "out_treatment_vs_high/mds_plot_landscape_Combat.pdf"
			)
		
		# Return the normalisation data sets
		my_df_mds %>%
			left_join(my_df_mds.RUV) %>%
			left_join(my_df_mds.Combat)
	} %>%
	
	# Attach annotation
	left_join(annot) %>%
	mutate_if(is.character, as.factor)

#############################################################
# DE ########################################################

design = 
	model.matrix(
		~ 
			d_adj %>% distinct(Sample, Label, W) %>% arrange(Sample) %>% pull(Label) +
			d_adj %>% distinct(Sample, Label, W) %>% arrange(Sample) %>% pull(W)
	) %>%
	magrittr::set_colnames(c("(Intercept)", "neoadjuvant", "W"))

# contrasts = 
# 	makeContrasts(
# 		MUvsWT=neoadjuvant-high,
# 		levels=design
# 	)

DE.obj <-
	d %>%
	
	# Filter f
	{
		(.)[
			(.)$genes$symbol %in% 
				(
					d_adj %>% 
						distinct(symbol) %>% 
						pull(symbol) %>% 
						as.character
				),
			]
	} %>%
	
	# Add group info to d
	{
		d_with_group = (.)
		d_with_group$samples$group = as.factor(design[,2])
		d_with_group
	} %>%
	
	# Calculate statistics
	calcNormFactors(method="TMM") %>%
	estimateGLMCommonDisp(design) %>%
	estimateGLMTagwiseDisp(design) %>%
	
	# Downstream analyses
	{
		y = (.)
		
		# Output list
		list(
			# EdgeR
			top = {
				y %>%
					glmFit(design) %>%
					glmLRT(coef = 2) %>%
					topTags(n=999999) %$%
					table %>%
					as_tibble() %>%
					
					# Mark DE genes
					mutate(is_de = FDR < 0.05 & abs(logFC) > 1) %>%
					# Write 
					{
						(.) %>% write_csv("out_treatment_vs_high/DE_table_RUV.csv")
						(.)
					} %>%
					
					# Plot smear
					{
						top = (.)
						
						plotSmear(
							y, 
							de.tags=
								top %>% 
								filter(is_de) %>% 
								pull(GeneID)
						) %>%
							do.call(bind_cols, .)	%>%
							mutate(symbol = y$genes$symbol) %>%
							left_join(top) %>%
							mutate(symbol = ifelse(is_de, symbol, "")) %>%
							{
								ggplot((.), aes(x = A, y = M,  label=symbol)) +
									geom_point(aes(color = is_de, alpha = is_de, size = is_de)) +
									ggrepel::geom_text_repel(
										size = 1, 
										point.padding = 0.3, 
										segment.size = 0.2,
										seed = 123
									) +
									scale_color_manual(values = c("FALSE" = "grey20", "TRUE" = "#db2523")) +
									scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
									scale_size_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
									theme_bw() +
									theme(
										panel.border = element_blank(), 
										axis.line = element_line(),
										panel.grid.major = element_line(size = 0.2),
										panel.grid.minor = element_line(size = 0.1),
										text = element_text(size=12),
										legend.position='none',
										aspect.ratio=1,
										strip.background = element_blank(),
										axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
										axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
									) +
									xlab("Average log count per million") + ylab("Log fold change") 
							}	 %>%
							ggsave(plot = .,
										 "out_treatment_vs_high/plot_smear_RUV.pdf",
										 useDingbats=FALSE,
										 units = c("mm"),
										 width = 89 ,
										 height = 89 
							)
						
						(.)
					} %>%
					
					# Heat map
					{
						top = (.)
						
						# # Top genes from MDS - NOT USED ANYMORE BUT HELPFUL
						# source("~/third_party_sofware/myMDS.R")
						# top_genes = my_df_mds %>%
						# 	dplyr::select(symbol, Sample, !!read_count) %>%
						# 	spread( Sample, !!read_count) %>%
						# 	as_matrix(rownames = "symbol") %>%
						# 	my.mds(gene.selection = "common") %$% 
						# 	top_genes 
						
						library(superheat)
						
						pdf("out_treatment_vs_high/heatmap_RUV.pdf", width = 183 * 0.0393701, height = 183*0.0393701, useDingbats=F)
						d_adj %>%
							left_join(annot) %>%
							unite(Sample_label, c("Sample", "Label"), remove = F) %>%
							dplyr::select(symbol, Sample_label, `value RUV log`) %>%
							filter(
								symbol %in% 
									(
										top %>% 
											filter(is_de) %>% 
											pull(symbol)
									)
							) %>%
							spread( Sample_label, `value RUV log`) %>%
							as_matrix(rownames = "symbol") %>%
							t() %>%
							scale() %>%
							t() %>%
							
							superheat(
								row.dendrogram = T,
								col.dendrogram = T,
								bottom.label.text.angle = 90,
								left.label.text.size = 1,
								bottom.label.text.size = 3,
								left.label.size = 0.1,
								grid.hline.size	= 0.1,
								grid.vline.size = 0.1
							) 
						dev.off()
						
						detach("package:superheat", unload=TRUE, force = T)
						
						(.)
					} %>%
					
					# GSEA
					{
						
						# source("https://bioconductor.org/biocLite.R")
						# biocLite("GSEABase")
						# biocLite("GSVA")
						obesity_signature <- as.character(unlist(read.table("obesity_signature.csv", quote="\"", comment.char="")))
						
						(.) %>%
							arrange(desc(abs(logFC))) %>%
							dplyr::select(symbol, logFC) %>%
							drop_na() %>%
							write_delim("out_treatment_vs_high/DE_path_obesity_GSEA.rnk", col_names = F, delim = "\t")
						
						system(sprintf("java -cp gsea2-2.2.2.jar -Xmx8g xtools.gsea.GseaPreranked -gmx %s -rnk %s -collapse false -mode Max_probe -norm meandiv -nperm 5000 -scoring_scheme weighted -include_only_symbols true -make_sets true -plot_top_x 100 -rnd_seed timestamp -set_max 500 -set_min 1 -zip_report false -out gsea -gui false", "obesity.gmt", "out_treatment_vs_high/DE_path_obesity_GSEA.rnk"))
						
						writeLines("Top genes present in the signature")
						
						(.) %>%
							arrange(desc(abs(logFC))) %>%
							left_join(
								readLines("obesity.gmt") %>% 
									strsplit("\t") %>% 
									as.data.frame() %>% 
									as_tibble() %>% 
									setNames("symbol") %>%
									mutate(is_in_obesity = TRUE)
							) %>%
							filter(is_in_obesity == TRUE) %>%
							head(n=10) %>%
							pull(symbol) %>%
							paste(collapse=" ") %>%
							writeLines()
						
						# Return original data frame
						(.)
					}
			},
			
			# EGSEA
			egsea.res = {
				library(AnnotationDbi)
				
				file_name = "out_treatment_vs_high/EGSEA_high_vs_treated.RData"
				
				er = switch(
					
					# Condition
					(!file.exists(file_name)) + 1,
					
					# If file exist
					{
						load(file_name)
						egsea.results
					},
					
					# If file does NOT exist
					y %>%
						voom(design, plot=FALSE) %>%
						
						# Run gene enrichment
						{
							v = (.)
							colnames(v$genes) = c("ENTREZID", "length", "SYMBOL")
							v$genes = v$genes[,c(1,3,2)]
							
							idx = buildIdx(entrezIDs=rownames(v), species="human")
							
							v %>%
								egsea(
									contrasts=2, 
									gs.annots=idx, 
									symbolsMap=
										v %$% 
										genes %>% 
										dplyr::select(1:2) %>%
										setNames(c("FeatureID", "Symbols")),
									baseGSEAs = egsea.base()[-c(6, 7, 8, 9)],
									sort.by="med.rank",
									num.threads = 1
								)
						} %>%
						{
							egsea.results = (.)
							save(egsea.results, file=file_name)
							(.)
						}
				) 
				
				detach("package:AnnotationDbi", unload=TRUE, force = T)
				
				# return after detach trick 
				er
			}
		)
		
		
	}


#############################################################
# Plot DE ###################################################

# Summary EGSEA
t = topSets(DE.obj$egsea.res, names.only=FALSE, number = Inf, verbose = FALSE)
t[grep("LIM_", rownames(t)), c("p.adj", "Rank", "med.rank", "vote.rank")]


d_adj %>% 
	
{
	top_de = DE.obj %$%
		top %>% 
		filter(is_de) 
	
	# Add annotation
	(.) %>% 
		inner_join(top_de) %>%
		
		# Shape raw and RUV counts
		mutate(`Counts RUV` = exp(`value RUV log`)) %>%
		dplyr::rename(`Counts raw` = value) %>%
		gather(is_normalised, `Read count`, c("Counts raw", "Counts RUV")) %>%
		
		# Set order factors
		mutate(symbol = factor(symbol, levels = top_de %>% pull(symbol)))
	
} %>%
	
	# Plot
{
	ggplot((.), aes(x=is_normalised, y = `Read count`, fill = Label)) +
		geom_boxplot(outlier.size = 0, lwd=0.2, position = position_dodge(width=0.8)) +
		geom_point(position=position_jitterdodge(dodge.width=0.8), size = 0.2, shape = 21 ) +
		#geom_jitter(height = 0, width = 0.2) +
		facet_wrap(~ symbol, scale="free") +
		scale_y_log10() +
		scale_fill_manual(values = c("Neoadjuvant" = "#e31e1e", "High" = "#999999")) +
		theme_bw() +
		theme(
			panel.border = element_blank(), 
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=9),
			legend.position="bottom",
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(t = 30, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.text.x = element_text(angle = 30, vjust = 0.5)
		)
} %>%
	ggsave(plot = .,
				 "out_treatment_vs_high/boxplot_DE_genes.pdf",
				 useDingbats=FALSE,
				 units = c("mm"),
				 width = 183 ,
				 height = 183/2*3 
	)

# top %>% filter(FDR<0.05) %>% mutate(`Fold change` = exp(abs(logFC))) %>%  summarise(median(`Fold change`), max(`Fold change`))

#############################################################
# Manual annotation #########################################

DE.obj %$%
	top %>%
	filter(is_de) %>%
	
	# Annotate
	left_join(
		biomaRt::getBM(
			filters=c("hgnc_symbol"),
			values = (.)$symbol , 
			attributes = c(
				'ensembl_gene_id', 
				'entrezgene','hgnc_symbol',
				'description',
				'name_1006',    'namespace_1003', 'definition_1006'
			),
			mart = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
		) %>%
			as_tibble() %>%
			dplyr::rename(symbol=hgnc_symbol),
		by = "symbol"
	) %>%
	mutate(`Recurrent patterns` = NA) %>%
	mutate(
		`Recurrent patterns` = 
			ifelse(
				grepl("lipid|steroid|hormone|fat|fatty|adipose|sperm|estradiol|insulin|sex|corticoid|cortisol|spermatogenesis|thyroid", name_1006),
				"Hormone/fat homeostasis",
				`Recurrent patterns`
			)
	) %>%
	mutate(
		`Recurrent patterns` = 
			ifelse(
				grepl("axon|neuron|neuro|synapse|brain", name_1006),
				"Neural",
				`Recurrent patterns`
			)
	) %>%
	mutate(
		`Recurrent patterns` = 
			ifelse(
				grepl("immune|complement|leukocyte|inflammation|infection|neutrophil|myeloid", name_1006),
				"Inflammation",
				`Recurrent patterns`
			)
	) %>%
	
	{
		(.) %>% write_csv("Supplementary_table_1_anotation_DE.csv")
		(.)
	} %>%
	
	filter(!is.na(`Recurrent patterns`)) %>%
	distinct(symbol, logFC, logCPM, `Recurrent patterns`) %>%
	group_by(symbol) %>%
	summarise(`Recurrent patterns` = paste(`Recurrent patterns`, collapse=" | "), logFC = unique(logFC), logCPM = unique(logCPM)) %>%
	mutate_if(is.character, as.factor) %>%
	{
		df = (.)
		
		df %>% 
			filter(`Recurrent patterns` == "Hormone/fat homeostasis") %>%
			mutate(x = runif(n(), -1, 0.5), y=runif(n(), -1, 0.5)) %>%
			bind_rows(
				df %>% 
					filter(`Recurrent patterns` == "Hormone/fat homeostasis | Inflammation") %>%
					mutate(x = 1.2, y=runif(n(), -1, 0.5)) 
			) %>%
			bind_rows(
				df %>% 
					filter(`Recurrent patterns` == "Inflammation") %>%
					mutate(x = runif(n(), 2.5, 4), y=runif(n(), -1, 0.5)) 
			) %>%
			bind_rows(
				df %>% 
					filter(grepl("Neural", `Recurrent patterns`)) %>%
					mutate(x = runif(n(), 0.2, 1.7), y=runif(n(), -5, -3.5)) 
			)
	} %>%
	{
		set.seed(123)
		my_max = (.) %>% pull(logFC) %>% max
		ggplot((.), aes(x = x, y = y, label = symbol)) + 
			geom_point(aes(fill=logFC, size = logCPM), shape=21) +
			ggrepel::geom_text_repel(
				size = 3, 
				point.padding = 0.3, 
				segment.size = 0.2,
				seed = 123
			) +
			scale_fill_distiller(
				palette = "Spectral",
				na.value = 'white',
				direction = 1,
				limits=c(	-my_max,my_max)
			) +
			#scale_fill_brewer(palette = "Set1") +
			theme_bw() +
			theme(
				panel.border = element_blank(), 
				axis.line = element_line(),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.title =  element_blank(),
				axis.text = element_blank(),
				axis.ticks = element_blank(),
				text = element_text(size=12),
				legend.position="bottom",
				aspect.ratio=2/4.5,
				strip.background = element_blank()
			)
	} %>%
	ggsave(plot = .,
				 "out_treatment_vs_high/plot_DE_recurrent_pathways.pdf",
				 useDingbats=FALSE,
				 units = c("mm"),
				 width = 283 ,
				 height = 283 + 20 
	)


#############################################################
# Tissue composition analyses ###############################

tissue_composition = 
	d_adj %>%
	mutate( `Count RUV` = exp(`value RUV log`)) %>%
	dplyr::select(symbol, Sample, `Count RUV`) %>%
	spread(Sample, `Count RUV`) %>%
	group_by(symbol) %>% 
	filter(row_number() == 1) %>%
	ungroup() %>%
	{
		df = (.)
		
		list(
			
			# Cibersort
			cibersort = df %>%
			{
				file_name = sprintf("out_treatment_vs_high/read_count_for_cibersort.tab")
				write_delim((.), path = file_name, delim = "\t")
				file_name
			} %>%
				
				# Inference
			{
				CIBERSORT("~/PhD/deconvolution/ARMET_BK_Apr2017/comparison_methods/cibersort/LM22.txt",	(.))$proportions %>%
					as_tibble(rownames = "Sample") %>% 
					dplyr::select(-`P-value`, -Correlation, -RMSE)
			} %>%
				
				# test
			{
				proportions = (.)
				dd = (.) %>% left_join(d_adj %>% distinct(Sample, Label, W)) %>% dplyr::select(-Sample) 
				AL = (.) %>% dplyr::select(-Sample) %>% DirichletReg::DR_data()
				DirichletReg::DirichReg(AL ~ Label + W, data=dd) %>%
					summary %>%
					{
						
						res = (.)
						res %$% 
							coef.mat %>% 
							as_tibble(rownames="cov") %>% 
							mutate(
								`Cell type` = 
									res %$%	
									residuals %>% 
									colnames() %>% 
									sapply(., rep, 3) %>% 
									c
							) %>%
							filter(cov == "LabelNeoadjuvant") %>%
							mutate(FDR = p.adjust(`Pr(>|z|)`, "bonferroni")) %>%
							mutate(Sig = ifelse(FDR < 0.05, "*", "")) 
					} %>%
					
					# Unite proportions with test
					left_join(
						proportions %>%
							gather(`Cell type`, Proportion, -Sample) %>%
							left_join(d_adj %>% distinct(Sample, Label))
					) 
			},
			
			# ARMET-tc
			armet = 
			{
			browser()
				file_name = "out_treatment_vs_high/armet_high_vs_treated.RData"
				
				switch(
					
					# Condition
					(!file.exists(file_name)) + 1,
					
					# If file exist
					{
						load(file_name)
						armet.res
					},
					
					# If file does NOT exist
					{
						
						d$counts %>% as_tibble(rownames="GeneID") %>% left_join(d$genes %>% as_tibble %>% mutate(GeneID = GeneID %>% as.character)) %>% dplyr::select(-GeneID, -Length) %>% gather(sample, `read count`, -symbol) %>% drop_na %>% spread(sample, `read count`)%>%
							dplyr::rename(gene=symbol) %>%
							ARMET_tc(
								my_design = 
									design %>% 
									as_tibble() %>% 
									dplyr::mutate(
										sample = 
											df %>%
											dplyr::select(-symbol) %>% 
											colnames()
									),
								cov_to_test = "neoadjuvant"
							) %>%
							{
								armet.res = (.)
								save(armet.res, file=file_name)
								armet.res
							} %>%
							{
								# Plot the immune cell components for validation
								armet.res %$%
									tree %$%
									estimate_prop_with_uncertanties  %>%
									left_join(	annot %>% distinct(Sample, Label) %>% rename(sample = Sample)	) %>%
									ggplot(aes(x=sample, y=estimate, shape = ct, color = Label)) + 
									geom_errorbar(aes(ymin=.lower, ymax=.upper), color="grey") +
									geom_point() +  
									geom_tile("Immune cells for validation") +
									my_theme
								
								(.)
							}
					}
				)
			}
		)
	} %>%
	
	# plot Comparison
	{

		getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

		bind_rows(
			
			# Cibersort
			(.) %$% 
				cibersort %>%
				
				# Convert cell types to match
				left_join(
					read_csv("cell_type_conversion.csv") %>% dplyr::select(1,2) %>% dplyr::rename(`Cell type` = Cibersort)
				) %>%
				rowwise() %>%
				mutate(`Cell type` = ifelse(!is.na(ARMET), ARMET, `Cell type`)) %>%
				ungroup() %>%
				dplyr::select(-ARMET) %>%
				mutate(`Cell type` = gsub("_", " ", `Cell type`)) %>%
				mutate(`Cell type` = Hmisc::capitalize(`Cell type`)) %>%
				
				# Add significance to cell name
				unite(`Cell type sig`, c("Cell type", "Sig"), sep = " ", remove = F) %>%
				mutate(`Cell type sig` = str_trim(`Cell type sig`)) %>%
				dplyr::rename(sample = Sample) %>%
				mutate(Algorithm = "Cibersort") ,
			
			# ARMET
			(.) %$% 
				armet %>%
			{
				res = (.)
			
				res %$% 
					proportions %>% 
					rename(`Cell type` = ct, Proportion = absolute_proportion) %>%
					dplyr::rename(Label = neoadjuvant) %>%
					mutate(Label = as.factor(Label)) %>%
					left_join(
						res %$% 
							stats %>% 
							data.tree::ToDataFrameTable("name", "Direction", "Sig",   "Driver", "CI.low", "CI.high") %>% 
							as_tibble() %>%
							dplyr::rename(`Cell type` = name)
					) %>%
					rowwise() %>%
					mutate(Driver = paste(rep(Driver, 2), collapse="")) %>%
					ungroup() %>%
					unite(Sig, c("Sig", "Driver"), sep = " ") %>%
					mutate(Sig = str_trim(Sig)) %>%
					drop_na() %>%
					
					# Beautify labels
					mutate(`Cell type` = gsub("_", " ", `Cell type`)) %>%
					mutate(`Cell type` = Hmisc::capitalize(`Cell type`)) %>%
					
					# Add significance to cell name
					unite(`Cell type sig`, c("Cell type", "Sig"), sep = " ", remove = F) %>%
					mutate(`Cell type` = str_trim(`Cell type`)) %>%
					
					# Mutate Label to match Cibersort
					rowwise() %>%
					mutate(Label = ifelse(Label == 0, "High", "Neoadjuvant")) %>%
					ungroup()
			} %>%
				mutate(Algorithm = "ARMET") 
		) %>%
		
		#select best sample for valudation	
		{
				
				(left_join(
					(.) %>% select(Algorithm, sample, `Cell type`, Proportion, Label) %>%
						spread(Algorithm, Proportion),
					(.) %>% filter(Algorithm == "ARMET") %>%
						select(sample, `Cell type`, Label, .lower_absolute, .upper_absolute)
				)	 %>%
				#drop_na() %>%
				# group_by(sample) %>%
				# mutate(ARMET = ARMET / sum(ARMET), Cibersort = Cibersort / sum(Cibersort)) %>%
				# ungroup() %>%
				filter(grepl("Macrophage", `Cell type`) | `Cell type` == "Monocyte") %>%
				separate(`Cell type`,c("Cell type"), sep= " ") %>%
				group_by(`Cell type`, sample, Label) %>%
				summarise(
					ARMET = ARMET %>% sum,
					Cibersort = Cibersort %>% sum,
					.lower_absolute = .lower_absolute %>% sum,
					.upper_absolute = .upper_absolute %>% sum
				) %>%
				gather(Algorithm,Proportion, c("ARMET", "Cibersort")) %>%
				mutate(
					.lower_absolute = ifelse(Algorithm == "Cibersort", NA, .lower_absolute),
					.upper_absolute = ifelse(Algorithm == "Cibersort", NA, .upper_absolute)
				) %>%
				ggplot(aes(x=sample, y=Proportion, shape = `Cell type`, color = Label)) + 
				geom_errorbar(aes(ymin=.lower_absolute, ymax=.upper_absolute), color="grey") +
				geom_point() + 
				facet_grid(~ Algorithm) + 
				my_theme) %>%
				plot()
			
			(.)
		} %>%
			# plot
		{
			res = 
				(.) %>%
				mutate(
					`Hypothesis test label` = ifelse(
						Algorithm == "Cibersort",
						sprintf("FDR = %s", FDR %>% formatC(format = "e", digits = 2)),
						sprintf("CI = %s/%s", CI.low, CI.high)
					)
				)
			
			browser()
			res %>%
			{
				ggplot((.), aes(x=Label, y = Proportion, fill = `Cell type`, label = Sig)) +
					geom_boxplot(outlier.size = 0) +
					geom_jitter(size = 0.2, height = 0) +
					facet_wrap( 
						`Cell type` ~ Algorithm + `Hypothesis test label`, 
						scales = "free", 
						ncol = 8, 
						labeller = label_wrap_gen(width=10)
					) +
					#scale_fill_manual(values = c("Neoadjuvant" = "#e31e1e", "High" = "#999999")) +
					scale_fill_manual(values = getPalette( res %>% distinct(`Cell type`) %>% nrow )) +
					scale_y_continuous(labels = scales::scientific_format(digits = 1)) +
					theme_bw() +
					theme(
						panel.border = element_blank(), 
						axis.line = element_line(),
						panel.grid.major = element_line(size = 0.2),
						panel.grid.minor = element_line(size = 0.1),
						text = element_text(size=9),
						legend.position="bottom",
						strip.background = element_blank(),
						strip.text = element_text(margin = margin(0,0,0,0, "cm")),
						axis.title.x  = element_text(margin = margin(t = 30, r = 10, b = 10, l = 10)),
						axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
						axis.text.x = element_text(angle = 30, vjust = 0.5)
					)
			} %>%
				ggsave(plot = .,
							 "out_treatment_vs_high/tissue_composition_cibersort_vs_ARMET.pdf",
							 useDingbats=FALSE,
							 units = c("mm"),
							 width = 256 ,
							 height = 333
				)
		}
		
		# Return all results
		(.)
		
	} %>%
	
	# Plot ARMET polar plot
	{
		
		(.) %$% 
			armet %>% 
			ARMET_plotPolar(
				size_geom_text = 2.5,
				my_breaks=c(0, 0.01, 0.1,0.5,1),
				prop_filter = 0.01,
				barwidth = 0.5, barheight = 3,
				legend_justification = 0.76
			) %>%
			ggsave(plot = .,
						 "out_treatment_vs_high/tissue_composition_ARMET_polar.pdf",
						 useDingbats=FALSE,
						 units = c("mm"),
						 width = 98 ,
						 height = 98 + 25
			)
		
		# Return all results
		(.)
		
	}






#############################################################
# qRT-PCR validation ########################################

# Plot of other potential markers
c(
	"CYP1A1", 
	"IGKV1D-39", 
	"IGKV1-39", 
	"MMP8", 
	"SLC16A12", 
	"ART3", 
	"DIO2", 
	"OR51E2", 
	"MUC16", 
	"TPRG1", 
	"IL4I1"
) 


	
	
DE.obj$top %>% 
	filter(abs(logFC)>1.5 ) %>% 
	head(n=40) %>% 
	pull(symbol)%>%

{
	
	gene_to_validate = (.)
	
	# Annotate
	biomaRt::getBM(
		filters=c("hgnc_symbol"),
		values = gene_to_validate, 
		attributes = c(
			'ensembl_gene_id', 
			'entrezgene','hgnc_symbol',
			'description',
			'name_1006',    'namespace_1003', 'definition_1006'
		),
		mart = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
	) %>%
		as_tibble() %>%
		dplyr::rename(symbol=hgnc_symbol) %>%
		filter(namespace_1003 %in% c("biological_process", "molecular_function")) %>%
		select(symbol, description, name_1006) %>%
		{
			(.) %>% write_csv("out_treatment_vs_high/annotation_for_gene_validation.csv")
			(.)
		}
	
	d_adj %>% 
		filter(symbol %in% gene_to_validate) %>%
		
		# Shape raw and RUV counts
		mutate(`Counts RUV` = exp(`value RUV log`)) %>%
		dplyr::rename(`Counts raw` = value) %>%
		gather(is_normalised, `Read count`, c("Counts raw", "Counts RUV")) %>%
		
		# Set order factors
		mutate(symbol = factor(symbol, levels = gene_to_validate)) %>%
		
		# Plot
		ggplot(aes(x=is_normalised, y = `Read count`, fill = Label)) +
		geom_boxplot(outlier.size = 0, lwd=0.2, position = position_dodge(width=0.8)) +
		geom_point(position=position_jitterdodge(dodge.width=0.8), size = 0.2, shape = 21 ) +
		#geom_jitter(height = 0, width = 0.2) +
		facet_wrap(~ symbol, scale="free") +
		scale_y_log10() +
		scale_fill_manual(values = c("Neoadjuvant" = "#e31e1e", "High" = "#999999")) +
		my_theme
}
	
	
# Validated genes	
readxl::read_excel("PCR_validation/Ryan qPCR fat analysis Oct18.xlsx", sheet = "raw") %>%
	dplyr::rename(CSDC2 = CSCD5) %>%
	bind_rows(
		(.) %>% head(n=3) %>% mutate_if(is.character, function(.) "HG02")  %>% mutate_if(is.numeric, function(.) NA)
	) %>%
	gather(symbol, `qRT-PCR CT value`, -`Sample Name`) %>%
		
	# Plot raw values
	{
		(
			(.) %>%
				ggplot(., aes(x=symbol, y=`qRT-PCR CT value`, color=`Sample Name`)) + 
				geom_point() +
				my_theme 
		) %>%
		plot()
	
		(.)
	} %>%
		
	# Filter outliers
	mutate(is_outlier = `qRT-PCR CT value` < 20 | is.na(`qRT-PCR CT value`)) %>%
	group_by(`Sample Name`, symbol) %>%
	mutate(how_many_outliers = is_outlier %>% sum ) %>%
	ungroup() %>%
	filter(!is_outlier | how_many_outliers == 3) %>%
	
	# Calculate normalised means
	group_by(`Sample Name`, symbol) %>%
	summarise(`qRT-PCR CT value mean` = mean(`qRT-PCR CT value`, na.rm=T)) %>%
	ungroup() %>%
		
	# Add control value
	left_join(
		readxl::read_excel("PCR_validation/Ryan qPCR fat analysis Oct18.xlsx", sheet = "average") %>%
			select(`Sample Name`, Control)
	) %>%
		
	# Normalise
	mutate(`qRT-PCR CT value mean normalised log` = -(`qRT-PCR CT value mean` - Control)) %>%
	mutate(`qRT-PCR CT value mean normalised` = 2^`qRT-PCR CT value mean normalised log`) %>%
		
	# Integrate the gene LYZ
	bind_rows(
		readxl::read_excel("PCR_validation/Ryan qPCR fat analysis Oct18.xlsx", sheet = "previousPatRun") %>%
			mutate(`qRT-PCR CT value mean normalised log` = `qRT-PCR CT value mean normalised` %>% log2)
	) %>%
	
	# Plot normalised values and return
	{
		(
			(.) %>%
			ggplot(aes(x=symbol, y=`qRT-PCR CT value mean normalised`, color=`Sample Name`)) + 
			geom_point() +
			scale_y_log10() +
			my_theme
		) %>%
		print()
		
		(.)
	} %>%
		
	# Integrate RNA seq data
	filter(!symbol %in% c("GAPDH", "TBP", "POL2")) %>%
	mutate(`Sample Name` = gsub("HG0", "HG", `Sample Name`)) %>%
	left_join(
		read_csv("PCR_to_RNAseq_sampleID_conversion.csv") %>%
			mutate(prefix="PP") %>%
			unite(Sample, c(Sample, prefix), sep	="")
	) %>%
	select(symbol, `qRT-PCR CT value mean normalised log`, Sample) %>%
	
	# Add annotation
	left_join( d_adj %>% distinct(Sample, Label) ) %>%
		
	# Plot vs. RNA sequencing
	{
		(.) %>%
			left_join(
				d_adj %>% 
					dplyr::rename(`Counts RUV log` = `value RUV log`) %>%
					mutate(`Counts raw log` = value %>% `+` (1) %>% log)
			) %>%
			
			# Standardise gene abundance
			gather(Source, `Gene transcript abundance log`, c("Counts raw log", "Counts RUV log", "qRT-PCR CT value mean normalised log")) %>%
			group_by(symbol, Source) %>%
			mutate(`Gene transcript abundance log standardised` =  `Gene transcript abundance log` %>% scale %>% as.numeric) %>%
			ungroup() %>%

			# Plot
			{
				ggplot((.), aes(x=Source, y =`Gene transcript abundance log standardised`, fill = Label)) +
				geom_boxplot(outlier.size = 0, lwd=0.2, position = position_dodge(width=0.8)) +
				geom_point(position=position_jitterdodge(dodge.width=0.8), size = 0.2, shape = 21 ) +
				facet_wrap(~ symbol, scale="free") +
				scale_fill_manual(values = c("Neoadjuvant" = "#e31e1e", "High" = "#999999")) +
				my_theme
			} %>% 
			print()
		
		(.)
	} %>%
	
	# Plot t-test
	{
		(.) %>%
		mutate(Label = ifelse(Label == "High", "Naive", "Treated")) %>%
		{
			ggplot((.), aes(x=Label, y =`qRT-PCR CT value mean normalised log`, fill = Label)) +
			geom_boxplot(outlier.size = 0, lwd=0.2, position = position_dodge(width=0.8)) +
			geom_point(position=position_jitterdodge(dodge.width=0.8), size = 0.05, shape = 21 ) +
			# ggpubr::stat_compare_means(	method = "t.test", size=0.8) +
			facet_wrap(~ symbol, scales = "free") +
			scale_fill_manual(values = c("Treated" = "#e31e1e", "Naive" = "#999999")) +
			my_theme +
			theme(
				text = element_text(size=4),
				legend.key.size = unit(3, "mm"),
				axis.line = element_line(size = 0.2),
				axis.ticks = element_line(size = 0.2)
			) 
		}  %>%
		ggsave(plot = .,
		 "out_treatment_vs_high/PCR_validation_boxplot.pdf",
		 useDingbats=FALSE,
		 units = c("mm"),
		 width = 89 ,
		 height = 89
		)
		
		(.)
	} %>%

	# T-test
	group_by(symbol) %>%
	do(

		(.) %>% 
			mutate(
				`p-value` = t.test(`qRT-PCR CT value mean normalised log` ~ Label, .)$p.value / 2
			) %>%
			distinct(symbol, `p-value`)
	) %>%
	ungroup() %>%
	mutate(`p-value adjusted` = `p-value` %>% p.adjust("bonferroni")) %>%
  mutate_if(is.numeric, function(x) x %>% formatC(format = "e", digits = 2))
	
	
	


#############################################################
# BMI check #################################################

# Load BMI
read_csv("Patient height_weight_BMI_NMC_treated_vs_naive.csv") %>%
	dplyr::rename(`Sample Name` = `Patient ID` ) %>%
	mutate(`Sample Name` = gsub("HG0", "HG", `Sample Name`)) %>%
	left_join(
		read_csv("PCR_to_RNAseq_sampleID_conversion.csv") %>%
			mutate(prefix="PP") %>%
			unite(Sample, c(Sample, prefix), sep	="")
	) %>%
	
	# Load CAPRA
	left_join( read_csv("treated_vs_high_CAPRA.csv") ) %>%
	
	# Add annotation
	left_join( d_adj %>% distinct(Sample, Label) ) %>%
	mutate(Label = ifelse(Label == "High", "Naive", "Treated")) %>%
	mutate(Labels = as.factor(Label)) %>%
	
	# Gather for facetting
	dplyr::rename(`CAPRA` = CAPRA_HR, `CAPRA-S` = CAPRA_S_HR) %>%
	gather(`Clinical variable`, Value, c("BMI", "CAPRA", "CAPRA-S")) %>%
	
	# Plot
	{
		(.) %>%
		ggplot(aes(x = Label, y = Value, fill = Label)) +
		geom_boxplot(outlier.size = 0, lwd=0.2, position = position_dodge(width=0.8)) +
		geom_point(position=position_jitterdodge(dodge.width=0.8), size = 0.05, shape = 21 ) +
		ggpubr::stat_compare_means(	method = "t.test") +
		facet_wrap(~ `Clinical variable`, scales = "free") +
		scale_fill_manual(values = c("Treated" = "#e31e1e", "Naive" = "#999999")) +
		my_theme
	} %>%
	ggsave(
		plot = .,
		"out_treatment_vs_high/Clinical_values_boxplot.pdf",
		useDingbats=FALSE,
		units = c("mm"),
		width = 189 ,
		height = 189/3+50
	)
