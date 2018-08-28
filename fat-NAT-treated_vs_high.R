setwd("~/PhD/fat-NAT")
set.seed(123)

# Libraries
library(hgu95av2.db)
library(tidyverse)
source("~/third_party_sofware/utilities_normalization.R")
library(sva)

#############################################################
# Build data set ############################################

if(0){
	bam.files <- dir("alignment_hg19", pattern=".sam", full.names=T, recursive=T)
	bam.counts <- featureCounts(bam.files,countMultiMappingReads=TRUE,fraction=TRUE, annot.inbuilt="hg19", nthreads=40)
	d <- DGEList(counts=bam.counts$counts, genes=bam.counts$annotation[,c("GeneID","Length")])
}
load("counts_paired_end.RData")

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
	mutate_if(is.character, as.factor) 

d = d[,annot %>% pull(Sample) %>% as.character()]

#############################################################
# MDS plots #################################################

# Function that takes annotated tibble and save MDS plot for the first 8 PC
tbl_to_MDS_plot = function(my_df_mds, annot, read_count = "Read count", file_name){
	
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
	dplyr::select(symbol, Sample, `value normalised`) %>%
	mutate(`value normalised log` = log(`value normalised` + 1)) %>%
	
	# MDS calculation
	{
		
		tbl_to_MDS_plot((.), annot, read_count = "value normalised log", "mds_plot_landscape.pdf")
		
		# Landscape
		my_df_mds = (.)
		
		foreach(
			components = list(c(1, 2), c(3, 4), c(5, 6), c(7, 8)), 
			#my_df_mds = (.),
			.combine = bind_rows
		) %do% {
			my_df_mds %>%
				dplyr::select(symbol, Sample, `value normalised log`) %>%
				spread( Sample, `value normalised log`) %>%
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
			
		# Plot
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
					 sprintf("mds_plot_landscape.pdf"),
					 useDingbats=FALSE,
					 units = c("mm"),
					 width = 183 ,
					 height = 183 
		)
		

		# MDS plot after batch correction
		library(ruv)
		source("~/third_party_sofware//RUV4.R")
		
		hk_600 = 
			as.character(read.table("~/PhD/deconvolution/hk_600.txt", quote="\"", comment.char="")[,1]) %>% 
			intersect(my_df_mds %>% pull(symbol))
browser()
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
					# Recalculated corrected expresison
					matrix(
						pmax(
							Y %>% t() - 
							(.)$W %*% (.)$alpha_all,
							0
						),
						nr=nrow(Y %>% t()),
						ncol=ncol(Y %>% t())
					)
				} %>%
				t() %>%
				magrittr::set_rownames(Y %>% rownames()) %>%
				magrittr::set_colnames(Y %>% colnames())
			} %>%
		
			# Annotate
			as_tibble(rownames="symbol") %>%
			gather(Sample,`value RUV log`, -symbol ) %>%
			left_join(my_df_mds) %>%
			left_join(annot) %>%
			mutate(Sample = gsub("PP", "", as.character(Sample))) %>%
			mutate_if(is.character, as.factor) %>%


		
	}
	




	# Batch correction
	{
		my_df = (.)
		
		my_df %>% 
			dplyr::select(symbol, sample, `value normalised log`) %>%
			spread(sample, `value normalised log`) %>%
			{
				mat = (.) %>% dplyr::select(-symbol) %>% as.matrix
				rownames(mat) = (.) %>% pull(symbol)
				mat
			} %>%
			ComBat(
				batch=my_df %>% distinct(sample, batch) %>% pull(batch), 
				mod=model.matrix(~ my_df %>% distinct(sample, CAPRA_groups) %>% pull(CAPRA_groups))
			) %>%
			{
				
				if(ct=="F")
					(.) %>%
					ComBat(
						batch=my_df %>% distinct(sample, ClinStageT) %>% pull(ClinStageT), 
						mod=model.matrix(~ my_df %>% distinct(sample, CAPRA_groups) %>% pull(CAPRA_groups))
					)
				else (.)
			} %>%
			
			# MDS calculation
			{
				
				my_df_mds = (.)
				
				foreach(
					components = list(c(1, 2), c(3, 4), c(5, 6), c(7, 8)), 
					#my_df_mds = (.),
					.combine = bind_rows
				) %do% {
					my_df_mds %>%
						plotMDS(dim.plot = components) %>% 
						{
							tibble(sample = names((.)$x), x = (.)$x, y = (.)$y, PCx = components[1], PCy = components[2])
						}
				} %>%
					left_join(my_df %>% distinct(sample, ClinStageT, batch, CAPRA_groups, CAPRA_TOTAL, sample_label)) %>%
					mutate(ct = ct)
			}
	} %>%
	
	# Plot and save





