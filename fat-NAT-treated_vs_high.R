setwd("~/PhD/fat-NAT")
set.seed(123)

# Libraries
library(tidyverse)
library(magrittr)
library(foreach)
source("~/third_party_sofware/utilities_normalization.R")
library(sva)

# RUV
library(ruv)
source("~/third_party_sofware//RUV4.R")

# Heatmap
library(ggdendro)

library(RColorBrewer)

# to_matrix function
source("https://gist.githubusercontent.com/stemangiola/39b08d529157ce59a5ff5dc1653951c5/raw/4538f8f5665810d3625fba8e06cbedce6486e887/to_matrix.R")

# EGSEA
library(EGSEA)
library(EGSEAdata)

# Safe detach
detach("package:AnnotationDbi", unload=TRUE, force = T)
detach("package:hgu95av2.db", unload=TRUE, force = T)

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
		org.Hs.eg.db,
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

#############################################################
# MDS plots #################################################

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
		
		getPalette = colorRampPalette(brewer.pal(9, "Set1"))
		
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
	
# Plot and save

design = 
	model.matrix(
		~ 
			0 +
			d_adj %>% distinct(Sample, Label, W) %>% arrange(Sample) %>% pull(Label) +
			d_adj %>% distinct(Sample, Label, W) %>% arrange(Sample) %>% pull(W)
	) %>%
	magrittr::set_colnames(c("high", "neoadjuvant", "W"))

contrasts = 
	makeContrasts(
		MUvsWT=neoadjuvant-high,
		levels=design
	)

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
					glmLRT(contrast = contrasts) %>%
					topTags(n=999999) %$%
					table %>%
					as_tibble() %>%
					
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
								filter(FDR < 0.05) %>% 
								pull(GeneID)
						) %>%
							do.call(bind_cols, .)	%>%
							mutate(symbol = y$genes$symbol) %>%
							left_join(top) %>%
							mutate(symbol = ifelse(FDR<0.05, symbol, "")) %>%
							{
								ggplot((.), aes(x = A, y = M,  label=symbol)) +
									geom_point(aes(color = FDR<0.05, alpha = FDR<0.05, size = FDR<0.05)) +
									ggrepel::geom_text_repel(
										size = 1.6, 
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
										legend.position="bottom",
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
										 width = 183 ,
										 height = 183 
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
						# 	to_matrix(rownames = "symbol") %>%
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
											filter(FDR < 0.05) %>% 
											pull(symbol)
									)
							) %>%
							spread( Sample_label, `value RUV log`) %>%
							to_matrix(rownames = "symbol") %>%
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
					}
			},
	
			# EGSEA
			egsea.res = {
				library(AnnotationDbi)
				
				y %>%
				voom(design, plot=FALSE) %>%
					
				# Run gene enrichment
				{
					v = (.)
					colnames(v$genes) = c("ENTREZID", "length", "SYMBOL")
					v$genes = v$genes[,c(1,3,2)]
					browser()
					idx = buildIdx(entrezIDs=rownames(v), species="human")
					v %>%
						egsea(
							contrasts=contrasts, 
							gs.annots=idx, 
							symbolsMap=
								v %$% 
								genes %>% 
								dplyr::select(1:2) %>%
								setNames(c("FeatureID", "Symbols")),
							baseGSEAs = egsea.base()[-c(4)],
							sort.by="med.rank"
						)
				} %>%
				{
					save((.), file="out_treatment_vs_high/EGSEA_high_vs_treated.RData")
					(.)
				}
				
				detach("package:AnnotationDbi", unload=TRUE, force = T)
	
			}
		)
	}


# top %>% filter(FDR<0.05) %>% mutate(`Fold change` = exp(abs(logFC))) %>%  summarise(median(`Fold change`), max(`Fold change`))
		

# Manual annotation
top %>%
	# Annotate
	filter(FDR<0.05) %>%
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
				grepl("immune|complement|leukocyte|inflammation|infection|neutrophil|myeloid", name_1006),
				"Inflammation",
				`Recurrent patterns`
			)
	) %>%
	filter(!is.na(`Recurrent patterns`)) %>%
	distinct(symbol, logFC, `Recurrent patterns`) %>%
	group_by(symbol) %>%
	mutate(`Recurrent patterns` = ifelse(n() > 1, "Both", `Recurrent patterns`)) %>%
	distinct() %>%
	mutate_if(is.character, as.factor) %>%
	{
		df = (.)

		df %>% 
		filter(`Recurrent patterns` == "Hormone/fat homeostasis") %>%
		mutate(x = runif(n(), -1, 1), y=runif(n(), -1, 1)) %>%
		bind_rows(
			df %>% 
				filter(`Recurrent patterns` == "Both") %>%
				mutate(x = 1.2, y=runif(n(), -1, 1)) 
		) %>%
		bind_rows(
			df %>% 
				filter(`Recurrent patterns` == "Inflammation") %>%
				mutate(x = runif(n(), 1.5, 3.5), y=runif(n(), -1, 1)) 
		)
	} %>%
	ggplot(aes(x = x, y = y, label = symbol)) + 
	geom_point(aes(fill=`Recurrent patterns`), shape=21, size = 3) +
	ggrepel::geom_text_repel(
		size = 3, 
		point.padding = 0.3, 
		segment.size = 0.2,
		seed = 123
	) +
	scale_fill_brewer(palette = "Set1") +
	theme_bw() +
	theme(
		panel.border = element_blank(), 
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		aspect.ratio=2/4.5,
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

	