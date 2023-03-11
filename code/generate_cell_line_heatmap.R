library(readr)
library(pheatmap)
library(tidyverse)

#Reading Data and wrangling
data_cell_line <- read_csv("data/processed_data/drug_response.csv")
data_cell_line<-na.omit(data_cell_line)

data_cell_line_wide<- data_cell_line[,-1] %>% pivot_wider(names_from = CELL_LINE_NAME,
                                                     values_from = SYNERGY_OBS_EMAX)
data_matrix<-t(as.matrix(data_cell_line_wide[,5:23]))

# Create column names
names <- paste(data_cell_line_wide$ANCHOR_NAME,
                   data_cell_line_wide$LIBRARY_NAME,
                   data_cell_line_wide$ANCHOR_CONC,
                   data_cell_line_wide$LIBRARY_CONC,
                   sep = "  ")

colnames(data_matrix)<-names

#Color pallete

rocket_r_colors <- c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7")

png(file='figures/drug_response_by_cell_line_heatmap.png',width = 1500, height = 800, units = "px",)
# Create a heatmap with row and column annotations
pheatmap(data_matrix,
         cluster_rows =F,
         cluster_cols=T,
              fontsize_row = 12,
         fontsize_col = 9.5,
         color = rocket_r_colors)

dev.off()