##### Run phylogenetic signal #####
library(tidyverse)
library(phytools)
### Load wrapper functions
source("https://raw.githubusercontent.com/jesusNPL/SpecEvolution/main/R/wrappers.R")

### hyperspectral data 
NEON_HS <- read_csv("data/NEON_airborne_spectra_CLEAN_NAS.csv")

### Load Phylogeny 
phyloNEON <- read.nexus("data/phylo_ALL_NEON_S3.nex")

### Match phylogeny and HS data
phyloNEON <- drop.tip(phyloNEON, 
                      setdiff(phyloNEON$tip.label, NEON_HS$species))

### Check info
setdiff(phyloNEON$tip.label, NEON_HS$species)

### Run Blomberg's K 
# run only 10 bands
NEON_HS_PhyloSig <- demon_K(spectra = data.frame(NEON_HS), # spectra in data.frame
                            tree = phyloNEON, # phylogeny
                            nSIM = 1000, # number of simulations
                            nBands = 10) # number of bands

### Read results
NEON_HS_PhyloSig <- read_csv("output/NEON_K_phylosig_CLEAN.csv")

### Transform wide to long format
NEON_HS_long <- NEON_HS %>% 
  pivot_longer(!species, 
               names_to = "wavelength", 
               values_to = "reflectance") 

### Remove water bands
wb1 <- c(1340, 1445)
wb2 <- c(1790, 1955)

NEON_HS_long_masked <- NEON_HS_long

# Mask out all values within each of the two atmospheric absorbtion bands
NEON_HS_long_masked[NEON_HS_long_masked$wavelength > 
                           wb1[1] & NEON_HS_long_masked$wavelength < wb1[2], ]$reflectance <- NA 

NEON_HS_long_masked[NEON_HS_long_masked$wavelength > 
               wb2[1] & NEON_HS_long_masked$wavelength < wb2[2], ]$reflectance <- NA 

### Plot mean species and Blomberg's K values
library(ggnewscale)

p_ps_K <- NEON_HS_long %>% 
  ggplot(aes(x = as.numeric(wavelength), 
             y = reflectance, 
             colour = factor(species)
  )) +
  geom_line(linewidth = 0.5, alpha = 0.5, show.legend = FALSE) + 
  scale_colour_viridis_d(option = "turbo", direction = -1) + 
  # add a new scale color
  new_scale_color() + 
  # plot phylogenetic signal 
  geom_point(aes(x = as.numeric(WL_num), 
                 y = extra_y, 
                 color = species, 
                 size = `Blomberg's K`), # Pval_cat = species
             data = NEON_HS_PhyloSig , 
             inherit.aes = FALSE) + 
  scale_color_manual(name = "P-value", 
                     labels = c("< 0.05", "> 0.05"), 
                     values = c("no_diff" = "black", "diff" = "red")) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  xlab("Wavelength (nm)") + 
  scale_y_continuous(
    # Features of the first axis
    name = "Reflectance",
    # Add a second axis and specify its features
    # sec.axis = sec_axis(trans = ~.*1, name = "Blomberg's K")
  ) + 
  annotate("rect", 
           xmin = 400, xmax = 700, # Visible
           ymin = -Inf, ymax = Inf,  fill = "grey2", alpha = 0.1) + 
  annotate("text", x = 559:560, y = 0.2, label = "VIS", size  = 7) + 
  annotate("rect",
           xmin = 800, xmax = 1340, # NIR
           ymin = -Inf, ymax = Inf,  fill = "grey2", alpha = 0.1) + 
  annotate("text", x = 1070:1071, y = 0.2, label = "NIR", size  = 7) + 
  annotate("rect",
           xmin = 1445, xmax = 1790, # SWIR 1
           ymin = -Inf, ymax = Inf,  fill = "grey2", alpha = 0.1) + 
  annotate("text", x = 1620:1621, y = 0.2, label = "SWIR1", size  = 7) + 
  annotate("rect",
           xmin = 1955, xmax = 2400, # SWIR 2 
           ymin = -Inf, ymax = Inf,  fill = "grey2", alpha = 0.1) +  
  annotate("text", x = 2170:2171, y = 0.2, label = "SWIR2", size  = 7) + 
  theme_minimal() + 
  theme(
    legend.position = "right",
    legend.text = element_text(size = 18), 
    legend.title = element_text(size = 20), 
    axis.text = element_text(size = 17, colour = "black"),
    axis.title = element_text(size = 20), 
    plot.title = element_text(size = 20, vjust = -12, hjust = 0.05),
    plot.subtitle = element_text(size = 20, vjust = -10, hjust = 0.05)
  )

p_ps_K

### Save results
pdf(file = "output/Fig_NEON_HS_phyloSignal.pdf", 
    height = 7, width = 10)

p_ps_K

dev.off()
