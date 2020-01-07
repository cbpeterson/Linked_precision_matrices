pdf('./Case_Study/Output/Figure2.pdf')
source("./Case_Study/Figures/create_histogram.R")
dev.off()

pdf('./Case_Study/Output/Figure3.pdf')
source("./Case_Study/Figures/create_heatmap.R")
dev.off()

pdf('./Case_Study/Output/Figure4.pdf')
source("./Case_Study/Figures/create_full_networks.R")
dev.off()

pdf('./Case_Study/Output/Figure5.pdf')
source("./Case_Study/Figures/create_networks_by_lobe.R")
dev.off()
