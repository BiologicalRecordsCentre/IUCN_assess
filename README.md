# IUCN_assess
This repository contains code to perform IUCN assessments based on occupancy data

taxaAssessIUCN.Rmd
  - this is the Rmarkdown file which creates a pdf, with a header page and a 1 page IUCN classification summary for each species within a taxonomic group. The user needs to load it up and manually change lines 3 and 18-20.

code/taxaIUCN.R
  - the top level script which the Rmd file calls. It takes in taxa name and model name and saves a lot of local files, most of which are images of graphs and some metadata used in compiling the Rmd. If these files have already been created locally, it just passes back the metadata and the locations of all the images of the graphs. All the files are stored in a local folder called 'tmp', so if new data is loaded, or data is updated, you will need to delete the folder 'tmp' to force it to re-run all species.

code/plotOccRmdFunctions.R
  - lower level functions. These create the plot files which the Rmd inserts as images. There are multiple functions in this file.

