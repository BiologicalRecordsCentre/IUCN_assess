#devtools::install_github('BiologicalRecordsCentre/wrappeR')
#unzip("BRCmap_0.10.3.8.zip")
#install.packages('BRCmap', repos = NULL)
library(pbapply)
library(lubridate)
library(ggplot2)
library(tools)
library(stringr)
library(dplyr)
suppressWarnings({library(BRCmap)})
library(gridExtra)
library(grid)
library(png)
library(adehabitatHR)
library(raster)
#source('loadRData.R')

tetrad_counts <- function(year, df, species){
  tmp <- df %>% filter(years == year)
  data.frame(Year = year,
             Tetrad_Count = tmp$tetrad %>% unique() %>% length(),
             Species = species,
             stringsAsFactors = FALSE)
}

# Function to set convergence colours.  Converged -> blue, not converged -> red
#   Arguments:
#     vec = converged vector of TRUE/FALSE for converged yes or no
#   Returns:
#     colours = vector of colours for plotting
set_colours <- function(vec){
  if(all(vec %>% as.logical())){
    return('Blue')
  } else if(all(!(vec %>% as.logical()))){
    return('Red')
  } else {
    return(c('Blue','Red'))
  }
}

create_occ_plot <- function(df, species, all_years){
  occupancy <- df[grepl('^psi.fs\\[[0-9]+\\]$',rownames(df)),]
  occupancy <- occupancy %>% data.frame() %>%
    dplyr::select(c('mean','X2.5.','X97.5.','Rhat'))
  occupancy$converged <- (occupancy$Rhat < 1.1)
  rownames(occupancy) <- str_extract(rownames(occupancy), '[0-9]+')
  occupancy$Year <-
    min(all_years) + (row.names(occupancy) %>% as.numeric()) - 1
  occupancy$converged <- factor(occupancy$converged, levels = c(TRUE, FALSE))
  
  colours <- set_colours(occupancy$converged)
  
  # Create smooth line through occupancy using spline
  spline_25  <- as.data.frame(spline(occupancy$Year, occupancy$X2.5.))
  spline_975 <- as.data.frame(spline(occupancy$Year, occupancy$X97.5.))
  
  # Build plot
  occ_plot <- ggplot() +
    geom_point(data = occupancy, aes(y = mean, x = Year, colour = converged)) +
    geom_vline(aes(xintercept = as.numeric(format(Sys.Date(), "%Y"))-10),
               alpha = .2) +
    geom_vline(aes(xintercept = as.numeric(format(Sys.Date(), "%Y"))-1),
               alpha = .2) +
    ylab('') + xlab('Year') +
    geom_line(data = as.data.frame(spline(occupancy$Year, occupancy$mean)),
              aes(y = y, x = x)) +
    geom_ribbon(aes(x = spline_25$x,
                    ymin = spline_25$y,
                    ymax = spline_975$y),
                alpha = 0.1) +
    scale_colour_manual(values = colours)
  occ_plot_zerone <- occ_plot +
    geom_ribbon(aes(x = c(seq(as.numeric(format(Sys.Date(), "%Y"))-10,
                              as.numeric(format(Sys.Date(), "%Y"))-1,
                              length.out = nrow(occupancy))),
                    ymin = 0, ymax = 1),
                fill = '#000080', alpha = .1) +
    coord_cartesian(ylim = c(0,1))
  occ_plot <- occ_plot +
    geom_ribbon(aes(x = c(seq(as.numeric(format(Sys.Date(), "%Y"))-10,
                              as.numeric(format(Sys.Date(), "%Y"))-1,
                              length.out = nrow(occupancy))),
                    ymin = min(occupancy$X2.5.), ymax = max(occupancy$X97.5.)),
                fill = '#000080', alpha = .1)
  return(list(occ_plot = occ_plot, occ_plot_zerone = occ_plot_zerone))
}

occ_outputs_function <- function(i, output_df){
  output_data_species_meta <-
    suppressWarnings({loadRfile(output_df$metadata_file[i])})
  all_years <- output_data_species_meta$min_year:output_data_species_meta$max_year
  species <- output_df$species[i]
  output_data_species <-
    suppressWarnings({loadRfile(output_df$output_file[i])})
  if('out' %in% names(output_data_species)){
    output_data_species <- output_data_species$out
  }
  occupancy_species <- output_data_species$BUGSoutput$sims.list$psi.fs %>%
    data.frame()
  converged <- output_data_species$BUGSoutput$summary

  occ_plot <- create_occ_plot(converged, species, all_years)
  occ_plot_fn <- file.path('tmp',paste0('occ_plot_',species,'.png'))
  ggsave(plot = grid.arrange(occ_plot$occ_plot,
                             occ_plot$occ_plot_zerone,
                             ncol=2),
         occ_plot_fn, width = 40, height = 10, units = 'cm')

  this_year <- format(Sys.Date(), "%Y")
  colyearnums <- as.numeric(str_extract(colnames(occupancy_species),'[0-9]+'))
  if(output_data_species_meta$max_year == this_year){
    occupancy_species <- occupancy_species[,!(colyearnums==max(colyearnums))]
  }
  if(ncol(occupancy_species)>=10){
    last10 <- occupancy_species[,(ncol(occupancy_species)-9):ncol(occupancy_species)]
  }
  last10$change_occ <- (last10[,ncol(last10)])/(last10[,1])
  last10_noinf <- last10[!(last10$change_occ==Inf),]
  occ_densplot <- ggplot(last10_noinf) +
    stat_density(aes(change_occ), geom = "line") +
    xlab('10 year change in occupancy') +
    ylab('Density') +
    scale_x_continuous(breaks = c(0,.1,.2,.35,.5,.6,.7,.85,1,1.5),
                       labels=c("0" = "0",".1" = "CR",
                                ".2" = "-80%",".35" = "EN",
                                ".5" = "-50%",".6" = "VU",
                                ".7" = "-30%",".85" = "LC",
                                "1" = "0%", "1.5" = "+50%")) +
    geom_vline(aes(xintercept = .7), alpha = .2) +
    geom_vline(aes(xintercept = .5), alpha = .2) +
    geom_vline(aes(xintercept = .2), alpha = .2) +
    geom_ribbon(aes(ymin = 0, ymax = Inf,
                    x = c(seq(0,.2,length.out = nrow(last10_noinf)))),
                fill = '#d12f1d', alpha = .2) +
    geom_ribbon(aes(ymin = 0, ymax = Inf,
                    x = c(seq(.2,.5,length.out = nrow(last10_noinf)))),
                fill = '#d19526', alpha = .2) +
    geom_ribbon(aes(ymin = 0, ymax = Inf,
                    x = c(seq(.5,.7,length.out = nrow(last10_noinf)))),
                fill = '#ebdb34', alpha = .2) +
    geom_ribbon(aes(ymin = 0, ymax = Inf,
                    x = c(seq(.7,2,length.out = nrow(last10_noinf)))),
                fill = '#67d126', alpha = .2) +
    coord_cartesian(xlim = c(0,2))
  ecdf10 <- ecdf(last10$change_occ)
  median <- quantile(x = last10$change_occ, probs = .5, na.rm = TRUE)
  
  occ_densplot_fn <- file.path('tmp',paste0('occ_densplot_',species,'.png'))
  
  ggsave(filename = occ_densplot_fn, plot = occ_densplot,
         width = 30, height = 10, units = 'cm')
  
  return(list(occ_densplot = occ_densplot_fn,
              occ_plot = occ_plot_fn,
              ecdf10 = ecdf10,
              median = median))
}

spec_input_function <- function(species, input_data){
  species_input <- input_data %>% filter(CONCEPT == species)
  if(!dir.exists('tmp')) dir.create('tmp')
  
  full_years <- unique(input_data$years) %>% sort()
  tetrad_input <-
    lapply(full_years,
           tetrad_counts,
           species_input,
           species) %>%
    bind_rows()
  tetrad_last_10 <- species_input
  all_but_this_year <- full_years[!(full_years %in% format(Sys.Date(), "%Y"))]
  last_10_years <-
    all_but_this_year[(length(all_but_this_year)-9):length(all_but_this_year)]
  obs_10_years <- species_input %>% data.frame() %>%
    filter(years %in% last_10_years)
  tetrad_10_years <- obs_10_years %>% pull(tetrad) %>% unique()
  num_tetrad_10_years <- length(tetrad_10_years)
  area_10_years <- num_tetrad_10_years * 4
  
  area_sq <- 0
  spdf <- BRCmap::gr2sp_points(tetrad_10_years)
  BR_TF <- grepl('^[A-Z]{2}[0-9]{2}[A-Z]{1}$',
                 as.character(spdf$GRIDREF))
  spdf_BR <- spdf[BR_TF,]
  spdf_EI <- spdf[!BR_TF,]
  if(sum(BR_TF)<5){
    BR <- FALSE
    area_sq_BR <- 0
  } else {
    BR <- TRUE
    spdf_BR$GRIDREF <- factor(x = rep(species,sum(BR_TF)),
                              levels = species)
    mcp_BR <- mcp(spdf_BR, percent = 100)
    inter_poly_BR <- suppressWarnings({
      raster::intersect(mcp_BR,
        UK_countries_lowres[as.character(UK_countries_lowres$REGION)!='Ireland',])
    })
    area_sq_BR <- sum(unlist(lapply(1:nrow(inter_poly_BR), FUN = function(j){
      area(inter_poly_BR[j,])
    })))/1000000
  }
  if(sum(!BR_TF)<5){
    EI <- FALSE
    area_sq_EI <- 0
  } else {
    EI <- TRUE
    spdf_EI$GRIDREF <- factor(x = rep(species,sum(!BR_TF)),
                              levels = species)
    mcp_EI <- mcp(spdf_EI, percent = 100)
    inter_poly_EI <- suppressWarnings({
      raster::intersect(mcp_EI,
        UK_countries_lowres[as.character(UK_countries_lowres$REGION)=='Ireland',])
    })
    area_sq_EI <- sum(unlist(lapply(1:nrow(inter_poly_EI), FUN = function(j){
      area(inter_poly_EI[j,])
    })))/1000000
  }
  area_sq <- sum(c(area_sq_BR,area_sq_EI))
  
  my_file <- file.path('tmp',paste0('spec_extent_',species,'.png'))
  png(filename = my_file, width = 1500, height = 1000)
  plot_GIS(UK_countries_lowres,
           xlim = c(-200000,600000),
           ylim = c(0,1050000))
  plotUK_gr_symbols(unique(obs_10_years$square), col = 'blue')
  while(!is.null(dev.list())) dev.off()
  
  my_file_mcp <- file.path('tmp',paste0('mcp_spec_extent_',species,'.png'))
  png(filename = my_file_mcp, width = 1500, height = 1000)
  plot_GIS(UK_countries_lowres,
           xlim = c(-200000,600000),
           ylim = c(0,1050000))
  if(BR){
    suppressWarnings({
      plot(inter_poly_BR,
           col = rgb(0,0,1,alpha=0.3),
           add = TRUE)
    })
  }
  if(EI){
    suppressWarnings({
      plot(inter_poly_EI,
           col = rgb(0,0,1,alpha=0.3),
           add = TRUE)
    })
  }
  plotUK_gr_symbols(unique(obs_10_years$square), col = 'blue')
  while(!is.null(dev.list())) dev.off()

  tet_plot <- file.path('tmp',paste0('tet_plot_',species,'.png'))
  tet_plot_p <- ggplot(tetrad_input) +
    geom_ribbon(aes(x = c(seq(as.numeric(format(Sys.Date(), "%Y"))-10.5,
                              as.numeric(format(Sys.Date(), "%Y"))-.5,
                              length.out = nrow(tetrad_input))),
                    ymin = 0, ymax = max(tetrad_input$Tetrad_Count)),
                fill = '#000080', alpha = .1) +
    geom_col(aes(x = Year, y = Tetrad_Count), fill = 'blue', alpha = .5) +
    ylab('Tetrad Count') + xlab('Year') +
    geom_vline(aes(xintercept = as.numeric(format(Sys.Date(), "%Y"))-10.5),
               alpha = .2) +
    geom_vline(aes(xintercept = as.numeric(format(Sys.Date(), "%Y"))-.5),
               alpha = .2)
  ggsave(plot = tet_plot_p,
         tet_plot, width = 15, height = 10, units = 'cm')
  
  myTet <- png::readPNG(tet_plot)
  myTet <- grid::rasterGrob(myTet, interpolate = TRUE)
  myImage <- png::readPNG(my_file)
  myImage <- grid::rasterGrob(myImage, interpolate = TRUE)
  myImage_mcp <- png::readPNG(my_file_mcp)
  myImage_mcp <- grid::rasterGrob(myImage_mcp, interpolate = TRUE)
  
  twocol <- paste0('tmp/',species,'_extent.png')
  ggsave(twocol, plot = grid.arrange(myImage, myTet, ncol = 2),
         width = 15, height = 5, units = 'cm')
  twocolmcp <- paste0('tmp/',species,'_mcp_extent.png')
  ggsave(twocolmcp, plot = grid.arrange(myImage_mcp, myTet, ncol = 2),
         width = 15, height = 5, units = 'cm')
  
  return(list(spec_extent = my_file,
              mcp_spec_extent = my_file_mcp,
              tet_plot = tet_plot,
              twocol = twocol,
              twocolmcp = twocolmcp,
              areas = list(area_10_years = area_10_years,
                           area_sq = area_sq)))
}