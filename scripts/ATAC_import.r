
# Import ATAC-seq data from 10x Genomics

import_all_from_10x <- function(){

	atac_base_dir = paste0(tg_atac_dir, "/mm_gastru/%s/outs/filtered_feature_bc_matrix/")
	
	import_from_10x(sprintf(atac_base_dir, "DNMT_TKO_KO"), genome = "mm10", id = "DNMT_TKO_KO")
    import_from_10x(sprintf(atac_base_dir, "DNMT_DKO_KO"), genome = "mm10", id = "DNMT_DKO_KO")
	import_from_10x(sprintf(atac_base_dir, "TET_TKO_KO"), genome = "mm10", id = "TET_TKO_KO")

    import_from_10x(sprintf(atac_base_dir, "DNMT_TKO_host"), genome = "mm10", id = "DNMT_TKO_host")
    import_from_10x(sprintf(atac_base_dir, "DNMT_DKO_host"), genome = "mm10", id = "DNMT_DKO_host")
	import_from_10x(sprintf(atac_base_dir, "TET_TKO_host"), genome = "mm10", id = "TET_TKO_host")

	import_from_10x(sprintf(atac_base_dir, "WT_sm"), genome = "mm10", id = "WT_sm")
	import_from_10x(sprintf(atac_base_dir, "WT_lg"), genome = "mm10", id = "WT_lg")
	
	import_from_10x(sprintf(atac_base_dir, "WT_E65"), genome = "mm10", id = "WT_E65")
	import_from_10x(sprintf(atac_base_dir, "WT_E70"), genome = "mm10", id = "WT_E70")
	import_from_10x(sprintf(atac_base_dir, "WT_E75"), genome = "mm10", id = "WT_E75")
								  
}		

# Read the sparse count matrix from 10x Genomics into mcATAC format SCC

read_SCCs <- function(raw_data){

	atac_base_dir = paste0(tg_atac_dir, "/mm_gastru/%s/outs/reads/")
	genotypes = c(#"WT_E65", "WT_E70", "WT_E75", 
				  "DNMT_TKO_KO", "DNMT_DKO_KO", "TET_TKO_KO", "WT_lg", "WT_sm",
				  "DNMT_TKO_host", "DNMT_DKO_host", "TET_TKO_host")
	
	for (genotype in genotypes){
		sc_counts <- scc_read(sprintf(atac_base_dir, genotype))
		raw_data[[paste(genotype, "scc", sep = "_")]] = sc_counts
	}	
		
	return(raw_data)
}

organize_md <- function(raw_data){
	
	models = c("host", "DNMT_TKO_KO", "DNMT_DKO_KO", "TET_TKO_KO", "WT")
	mds = c()
	for (model in models){
	
		# Load tables
		# MC2 outputs
		cells = read_h5ad(sprintf("./data/MC_models/multiome_%s_cells.h5ad", model))
        cli_alert_info(text = sprintf("%s: %s cells in h5ad.", model, nrow(cells$obs)))
                       
		mcells = read_h5ad(sprintf("./data/MC_models/multiome_%s_metacells_projected.h5ad", model))
		cell_md = read.csv("./data/ref/cell_md.csv") %>% select(id, cell_name, cell_name_mod) # id to cell name mapping

		# Unify all ~
			
		if ("metacell_name" %in% colnames(mcells$obs)) {
            mcells$obs = mcells$obs %>% select(-metacell_name)
        }
		mc_md = mcells$obs %>% rownames_to_column("metacell_name") %>% select(metacell_name, projected_type) %>% dplyr::rename(cell_type = projected_type) %>% 
			mutate_all(as.character)

		cell_md = cells$obs %>% rownames_to_column("id") %>% 
                mutate(origin_KO_host = paste(origin, KO_host, sep="_")) %>% 
                mutate(origin_KO_host = gsub("_WT$", "", origin_KO_host)) %>%
            select(id, origin_KO_host, metacell_name) %>% filter(id %in% rownames(cells$obs)) %>% 
			left_join(cell_md) %>% mutate(metacell_num = gsub("M|\\..*", "", as.character(metacell_name))) %>% left_join(mc_md) %>%
            select(id, cell_name, cell_name_mod, origin_KO_host, metacell_name, metacell_num, cell_type)
							 
		mds[[model]] = cell_md
	}
	
	raw_data$cell_md <- do.call(rbind, mds)
    rownames(raw_data$cell_md) = NULL

    cell_type_key = read.csv("./data/ref/cell_type_key.csv")
	raw_data$cell_type_key <- cell_type_key
	
	return(raw_data)
	
}

write_MCCs <- function(raw_data, group_by = "shorthand"){

    atac_base_dir = paste0(tg_atac_dir, "/mm_gastru/%s/outs/reads/")    

	cell_md <- raw_data$cell_md
	cell_type_key <- raw_data$cell_type_key

	genotypes = c("DNMT_TKO_KO", "DNMT_DKO_KO", "TET_TKO_KO", 
                  "DNMT_TKO_host", "DNMT_DKO_host", "TET_TKO_host",
                  "WT_sm", "WT_lg")
	
	for (genotype in genotypes){
        sc_counts <- raw_data[[paste(genotype, "scc", sep = "_")]]
        genotype_cells = sc_counts@cell_names[sc_counts@cell_names %in% 
                                              cell_md[cell_md$origin_KO_host == genotype, "cell_name"]]
		
	cli_alert_info(text = sprintf("%s: %s duplicate cells.", genotype, sum(duplicated(genotype_cells))))
						
    # Associate each cell with a group (either metacell or cell type AKA "shorthand")
    cell_to_group = cell_md %>% left_join(cell_type_key) %>%
        filter(origin_KO_host == genotype, cell_name %in% genotype_cells) %>% 
        dplyr::rename(group = {{group_by}}) %>% filter(!is.na(group)) %>% 
        select(cell_name, group) %>% 
        dplyr::rename(cell_id = cell_name) %>%
        mutate(group = ifelse(is.na(group), "Missing", group))
			
        mc_counts <- scc_to_mcc(sc_counts = sc_counts, 
                                cell_to_metacell = cell_to_group, 
                                ignore_metacells = c("Missing", "Excluded", "Outliers"))
							
		raw_data[[paste(genotype, "MCC", sep = "_")]] = mc_counts

    }						
	
	return(raw_data)
}

write_tracks <- function(raw_data, manual_choice_of_genotypes = NULL){

    if (!is.null(manual_choice_of_genotypes)){
        genotypes = manual_choice_of_genotypes
    } else {
	genotypes = c("DNMT_TKO_KO", "DNMT_DKO_KO", "TET_TKO_KO", 
                  "DNMT_TKO_host", "DNMT_DKO_host", "TET_TKO_host",
                  "WT_sm", "WT_lg")
    }
	for (genotype in genotypes){
		mc_counts <- raw_data[[paste(genotype, "MCC", sep = "_")]]
		mcc_to_tracks(mc_counts, glue::glue("mut_multi2.{genotype}"), overwrite = TRUE, create_marginal_track = TRUE, 
					  window_size = NULL, resolution = 20, normalize = FALSE)		
	}

	# Remove "^mc" from track names (function default)
	base_track_dir = file.path(mm10_trackdb, "tracks/mut_multi2/")
	for (genotype in genotypes){
		subdir = paste0(base_track_dir, genotype)
		setwd(subdir)
		
		files <- list.files()
		mc_files <- grep("^mc", files, value = TRUE)
		new_names <- sub("^mc", "", mc_files)
        for (dir_to_remove in new_names){
            if (file.exists(dir_to_remove)) {
                unlink(dir_to_remove, recursive = TRUE)
            }
        }
        
		file.rename(mc_files, new_names)
	}
}	

