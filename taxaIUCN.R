UK_only <- TRUE
my_model <- '2021_Francesca'
my_model <- '2021_Ellcur'
group <- file.path('/data-s3/occmods/Ladybirds')
group <- file.path('/data-s3/occmods/Ephemeroptera')
source('plotOccRmdFunctions.R')
create_group_data <- function(group, my_model, UK_only = TRUE){
  group_file <- paste0(basename(group),'_',my_model,'_filenames.rds')
  
  if(!file.exists(group_file)){
    input <- list.files(file.path(group,'input_data',my_model),
                        full.names = TRUE)
    if(length(input)>1){
      input <- input[!grepl('visitData',input)]
    }
    output_df <-
      data.frame(output_file = list.files(file.path(group,
                                                    'occmod_outputs',
                                                    my_model),
                                          full.names = TRUE),
                 stringsAsFactors = FALSE)
    output_df$species <- file_path_sans_ext(basename(output_df$output_file))
    
    # Extract daisy number
    output_df$daisy <- output_df$species %>%
      str_extract('(?<=[a-z]_)[0-9]+(?=_[0-9])') %>% as.numeric()
    if(all(is.na(output_df$daisy))){
      output_df <- output_df %>% dplyr::select(-daisy)
      output_df$metadata_file <- output_df$output_file
    } else {
      maxit <- max(output_df$daisy[!is.na(output_df$daisy)])
      keepers_max <- (is.na(output_df$daisy))|(output_df$daisy == maxit)
      output_df_max <- output_df %>% filter(keepers_max)
      output_df_max$species[!is.na(output_df_max$daisy)] <-
        str_extract(output_df_max$species[!is.na(output_df_max$daisy)], '^.*(?=_[0-9]+_[0-9])')
      output_df_max <- output_df_max %>% dplyr::select(-daisy)
      
      minit <- min(output_df$daisy[!is.na(output_df$daisy)])
      keepers_min <- (is.na(output_df$daisy))|(output_df$daisy == minit)
      output_df_min <- output_df %>% filter(keepers_min)
      output_df_min$species[!is.na(output_df_min$daisy)] <-
        str_extract(output_df_min$species[!is.na(output_df_min$daisy)], '^.*(?=_[0-9]+_[0-9])')
      output_df_min <- output_df_min %>% dplyr::select(-daisy)
      names(output_df_min)[1] <- c('metadata_file')
      
      output_df_minmax <- merge(output_df_max, output_df_min,
                                by.x = 'species', by.y = 'species',
                                all.y = TRUE, all.x = TRUE)
      
      output_df_unchained <- output_df %>% filter(is.na(output_df$daisy))
      output_df_unchained <- output_df_unchained %>% dplyr::select(-daisy)
      output_df_unchained$metadata_file <- output_df_unchained$output_file
      output_df <- bind_rows(output_df_minmax, output_df_unchained)
    }
    
    input_data <- loadRfile(input)
    # Subset to UK only
    if(any(names(input_data)=='SQ_1KM')){
      names(input_data)[names(input_data)=='SQ_1KM'] <- 'TO_GRIDREF'
    }
    if(UK_only){
      input_data <- input_data[grepl('^[A-Z]{2}[0-9]{4}$',
                                     as.character(input_data$TO_GRIDREF)),]
    }
    
    input_data$years <- lubridate::year(input_data$TO_STARTDATE)
    tet_quad <- BRCmap::det_tet_quad(input_data$TO_GRIDREF)
    input_data$tetrad <- tet_quad$TETRAD_GR
    input_data$square <- substr(input_data$tetrad, 1, nchar(input_data$tetrad)-1)
    all_input_species <- unique(input_data$CONCEPT)
    missing_species <- all_input_species[!(all_input_species) %in% output_df$species]
    if(length(missing_species)>0){
      missing_df <- data.frame(species = missing_species,
                               output_file = NA,
                               metadata_file = NA,
                               stringsAsFactors = FALSE)
    } else {
      missing_df <- NULL
    }
    output_df <- bind_rows(output_df, missing_df)
    output_df <- output_df %>% arrange(species)
    
    filenames <- pblapply(1:nrow(output_df), FUN = function(i){
      input_results <-
        spec_input_function(output_df$species[i], input_data)
      if(is.na(output_df$output_file[i])){
        output_results <- NA
      } else {
        output_results <- occ_outputs_function(i, output_df)
      }
      list(i = i,
           species = output_df$species[i],
           input_results = input_results,
           output_results = output_results)
    })
    group_data <- list(missing_species = missing_species,
                       filenames = filenames)
    saveRDS(group_data, group_file)
  } else {
    group_data <- readRDS(group_file)
  }
  return(group_data)
}