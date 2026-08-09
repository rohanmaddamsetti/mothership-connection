## ARG-MGE-PCN-analysis.R by Rohan Maddamsetti.
## analyse the distribution of antibiotic resistance genes (ARGs)
## on chromosomes versus plasmids in fully-sequenced genomes and plasmids
## in the NCBI RefSeq database first analyzed in Maddamsetti et al. 2025.

library(tidyverse)
library(cowplot)
library(data.table)
library(ggExtra)
library(scico) ## Nice scientific palettes: https://github.com/thomasp85/scico

## size threshold for calling small plasmids
SIZE_THRESHOLD <- 10000

################################################################################
## Regular expressions used in this analysis.

## "Tra" actually matches genes involved in conjugative transfer, which is better matched for plasmid functions rather than transposon functions!
## This is an oversight in the regular expressions used in the Zeevi et al. 2019 paper.

## match MGE genes using the following keywords in the "product" annotation
transposon.keywords <- "IS|transpos\\S*|insertion|Transpos\\S*|Tn[0-9]|tranposase|Tnp|Ins|ins"
plasmid.keywords <- "relax\\S*|conjug\\S*|Tra[A-Z]|Tra[0-9]|tra[A-Z]|mob\\S*|plasmid|chromosome partitioning|chromosome segregation|Mob\\S*|Plasmid|Rep|Conjug\\S*"
phage.keywords <- "capsid|phage|Tail|tail|head|tape measure|antiterminatio|Phage|virus|Baseplate|baseplate|coat|entry exclusion|Integrase|integrase"
other.HGT.keywords <- "excision\\S*|exonuclease|recomb|toxin|restrict\\S*|resolv\\S*|topoisomerase|reverse transcrip|intron|antitoxin|toxin|Toxin|Reverse transcriptase|hok|Hok|competence|addiction|type IV|conjugate transposon|post-segregation killing"

MGE.keywords <- paste(transposon.keywords, plasmid.keywords, phage.keywords, other.HGT.keywords, sep="|")

## antibiotic-specific keywords.
chloramphenicol.keywords <- "chloramphenicol|Chloramphenicol"
tetracycline.keywords <- "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
lincosamide.keywords <- "lincosamide|lincomycin|clindamycin|pirlimycin"
multidrug.keywords <- "Multidrug resistance|multidrug resistance"
beta.lactam.keywords <- "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\\S*"
glycopeptide.keywords <- "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
polypeptide.keywords <- "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
## DHFR inhibitors include sulfonamides and dihydropyrimidines like trimethoprim.
DHFR.inhibitor.keywords <- "trimethoprim|dihydrofolate reductase|dihydropteroate synthase|sulfonamide|Sul1|sul1|sulphonamide"
aminoglycoside.and.quinolone.keywords <- "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt|quinolone|Quinolone|oxacin|qnr|Qnr"
macrolide.keywords <- "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA|streptogramin"
antimicrobial.keywords <- "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\\S*"


antibiotic.keywords <- paste(
  chloramphenicol.keywords,
  tetracycline.keywords,
  lincosamide.keywords,
  multidrug.keywords,
  beta.lactam.keywords,
  glycopeptide.keywords,
  polypeptide.keywords,
  DHFR.inhibitor.keywords,
  aminoglycoside.and.quinolone.keywords,
  macrolide.keywords,
  antimicrobial.keywords, sep="|")

antibiotic.or.MGE.keywords <- paste(MGE.keywords,antibiotic.keywords,sep="|")


make_PCN_base_plot <- function(my.PCN.data) {
  ## PCN plot colored by PredictedMobility
  my.PCN.data |> 
    ggplot(aes(
      x = log10_replicon_length,
      y = log10_PIRACopyNumber,
      color = PredictedMobility)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    scale_color_manual(values=c("#fc8d62","#66c2a5","#8da0cb"), name="plasmid mobility") +
    ## make the points in the legend larger.
    guides(color = guide_legend(override.aes = list(size = 5))) +
    xlab("log10(length)")  +
    ylab("log10(copy number)") +
    theme(legend.position = "bottom") +
    theme(strip.background = element_blank()) +
    theme(
      axis.title.x = element_text(size=11),
      axis.title.y = element_text(size=11),
      axis.text.x  = element_text(size=11),
      axis.text.y  = element_text(size=11))
}


make_ARG_PCN_base_plot <- function(ARG.annotated.PCN.data) {
  
  ## Make the basic plot for Figure 1,
  ## before adding the marginal histograms or facetting
  ARG.annotated.PCN.data |>
    ggplot() +
    geom_point(
      data = subset(ARG.annotated.PCN.data, has_ARG == FALSE),
      aes(x = log10_replicon_length, y = log10_PIRACopyNumber),
      color = "grey80",
      size = 0.5, alpha=0.8) +
    geom_point(
      data = subset(ARG.annotated.PCN.data, has_ARG == TRUE),
      aes(x = log10_replicon_length, y = log10_PIRACopyNumber, color = ARG_count),
      size = 0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    scale_color_scico(palette = "berlin", name="ARG count") +
    ## make the points in the legend larger.
    xlab("log10(length)")  +
    ylab("log10(copy number)") +
    theme(legend.position = "bottom") +
    theme(strip.background = element_blank()) +
    theme(
      axis.title.x = element_text(size=11),
      axis.title.y = element_text(size=11),
      axis.text.x  = element_text(size=11),
      axis.text.y  = element_text(size=11))
}


make_transposon_PCN_base_plot <- function(transposon.annotated.PCN.data) {
  ## Make the basic plot for Figure 2,
  ## before adding the marginal histograms or facetting
  transposon.annotated.PCN.data |>
    ggplot() +
    geom_point(
      data = subset(transposon.annotated.PCN.data, has_transposon == FALSE),
      aes(x = log10_replicon_length, y = log10_PIRACopyNumber),
      color = "grey80",
      size = 0.5, alpha=0.8) +
    geom_point(
      data = subset(transposon.annotated.PCN.data, has_transposon == TRUE),
      aes(x = log10_replicon_length, y = log10_PIRACopyNumber, color = transposon_gene_count),
      size = 0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    scale_color_scico(palette = "vanimo", name="transposon count") +
    ## make the points in the legend larger.
    xlab("log10(length)")  +
    ylab("log10(copy number)") +
    theme(legend.position = "bottom") +
    theme(strip.background = element_blank()) +
    theme(
      axis.title.x = element_text(size=11),
      axis.title.y = element_text(size=11),
      axis.text.x  = element_text(size=11),
      axis.text.y  = element_text(size=11))
}


################################################################################
## Set up the key data structures for the analysis,
## gbk.annotation, in particular.

## import plasmid PCN data.
PCN.data <- read.csv(
  "../data/Maddamsetti2025-results/Source-Data/Fig1BC-Source-Data.csv") |>
  ## drop this useless column.
  select(-AnnotationAccession_right)

## annotate source sequences as plasmid or chromosome.
episome.database <- read.csv("../data/Maddamsetti2025-results/chromosome-plasmid-table.csv") |>
    as_tibble()

gbk.annotation <- read.csv(
  "../data/Maddamsetti2025-results/computationally-annotated-gbk-annotation-table.csv") |>
  as_tibble() |>
  ## refer to NA annotations as "Unannotated".
  mutate(Annotation = replace_na(Annotation,"Unannotated")) |>
  ## collapse Annotations into a smaller number of categories as follows:
  ## Marine, Freshwater --> Water
  ## Sediment, Soil, Terrestrial --> Earth
  ## Plants, Agriculture, Animals --> Plants & Animals
  mutate(Annotation = replace(Annotation, Annotation == "Marine", "Water")) |>
  mutate(Annotation = replace(Annotation, Annotation == "Freshwater", "Water")) |>
  mutate(Annotation = replace(Annotation, Annotation == "Sediment", "Earth")) |>
  mutate(Annotation = replace(Annotation, Annotation == "Soil", "Earth")) |>
  mutate(Annotation = replace(Annotation, Annotation == "Terrestrial", "Earth")) |>
  mutate(Annotation = replace(Annotation, Annotation == "Plants", "Plants & Animals")) |>
  mutate(Annotation = replace(Annotation, Annotation == "Agriculture", "Plants & Animals")) |>
  mutate(Annotation = replace(Annotation, Annotation == "Animals", "Plants & Animals")) |>
  ## get species name annotation from episome.database.
  left_join(episome.database) |>
  ## Annotate the genera.
  mutate(Genus = stringr::word(Organism, 1))

## assert that episome.database and gbk.annotation are consistent with each other.
stopifnot(nrow(filter(episome.database, !(AnnotationAccession %in% gbk.annotation$AnnotationAccession))) == 0)

## import the file containing plasmid proteins.
## Save memory by using the data.table package for import.
## import the sequence column to cross-check with duplicated.ARGs
## for Figure 5 analysis.
plasmid.proteins <- data.table::fread("../results/filtered-plasmid-proteins.csv") |>
  left_join(gbk.annotation)

plasmid.MGE.and.ARG.proteins <- plasmid.proteins |>
  filter(str_detect(product, antibiotic.or.MGE.keywords))

plasmid.ARGs <- plasmid.MGE.and.ARG.proteins |> 
  filter(str_detect(product, antibiotic.keywords))

plasmid.MGEs <- plasmid.MGE.and.ARG.proteins |> 
  filter(str_detect(product, MGE.keywords))

## count all ARGs.
plasmid.ARG.totals.df <- plasmid.ARGs |>
  summarize(ARG_count = n(),
            .by = c(AnnotationAccession, SeqID, SeqType)) |>
  arrange(AnnotationAccession, SeqID, SeqType)


################################################################################
## Figure 1. Experiments have shown that small multicopy plasmids can
## promote antibiotic resistance in several ways.
## This is a conceptual figure / illustration that I made separately by hand.

################################################################################
## Figure 2. ARGs are concentrated in larger plasmids, and rarer in small plasmids.

ARG.annotated.PCN.data <- PCN.data |>
  left_join(plasmid.ARG.totals.df) |>
  ## set NAs to zero in the ARG_count column
  mutate(ARG_count = replace_na(ARG_count, 0)) |>
  ## Note: this is the minus sign character "−" U+2212 in Unicode.
  mutate(has_ARG = ifelse(ARG_count > 0, TRUE, FALSE))

## 12,006 plasmids.
nrow(ARG.annotated.PCN.data)

## 3364 small plasmids.
nrow(filter(ARG.annotated.PCN.data, replicon_length < SIZE_THRESHOLD))

## 3870 plasmids have ARGs.
nrow(filter(ARG.annotated.PCN.data, ARG_count > 0))

##223 small plasmids (< 10 kB have ARGs)
nrow(filter(filter(ARG.annotated.PCN.data, ARG_count > 0), replicon_length < SIZE_THRESHOLD))

## scatterplot of log10(Normalized plasmid copy number) vs. log10(plasmid length).
## using a somewhat unusual plotting strategy to avoid overplotting points with ARGs
Fig2_base <- make_ARG_PCN_base_plot(ARG.annotated.PCN.data)

## add the marginal histogram to Figure A.
Fig2A <- ggExtra::ggMarginal(Fig2_base, margins="both")

## Figure 2B: facet by ecological annotation.
## By eye, looks like there may be some enrichment of ARGs on small plasmids
## in humans but unclear-- need statistics to be rigorous.
Fig2B <- Fig2_base + guides(color = "none") + facet_wrap(.~Annotation)


## Figure 2CD: show differences in mobility type.
## IMPORTANT: report in figure legend that NA datapoints were removed.
Fig2C_base <- ARG.annotated.PCN.data |>
  filter(!is.na(PredictedMobility)) |>
  filter(has_ARG==FALSE) |>
  make_PCN_base_plot() + ggtitle("ARG- plasmids")

## Get the legend.
Fig2CD_legend <- get_legend(Fig2C_base)
## now remove the legend from base figure.
Fig2C_base <- Fig2C_base + guides(color="none")
Fig2C <- ggExtra::ggMarginal(Fig2C_base, groupColour = TRUE, groupFill = TRUE, margins="both")

Fig2D_base <- ARG.annotated.PCN.data |>
  filter(!is.na(PredictedMobility)) |>
  filter(has_ARG==TRUE) |>
  make_PCN_base_plot() + ggtitle("ARG+ plasmids") + guides(color="none")
Fig2D <- ggExtra::ggMarginal(Fig2D_base, groupColour = TRUE, groupFill = TRUE, margins="both")

Fig2CD <- plot_grid(plot_grid(Fig2C, Fig2D, nrow=1, labels=c("C","D")),Fig2CD_legend, nrow=2,rel_heights=c(1,0.05))

Fig2 <- plot_grid(plot_grid(Fig2A, Fig2B, nrow = 1, labels=c("A","B")), Fig2CD, nrow = 2)

## Draft figure 2, showing that ARGs are largely on large conjugative plasmids,
## and rarely on small plasmids (but this is observed).
ggsave("../results/Fig2.pdf", Fig2, width=7.5,height=8.5)


################################################################################
## Supplementary Figures S1 - S12:
## analyse the data with annotated ARG classes.
## analyze distribution across different ARG classes using
## plasmid.ARG.class.totals.df.

## count ARGs per class.

plasmid.chloramphenicol.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, chloramphenicol.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "chloramphenicol",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.tetracycline.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, tetracycline.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "tetracyclines",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.lincosamide.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, lincosamide.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "lincosamides",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.multidrug.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, multidrug.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "multidrug",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.beta.lactam.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, beta.lactam.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "beta-lactams",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.glycopeptide.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, glycopeptide.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "glycopeptides",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.polypeptide.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, polypeptide.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "polypeptides",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.DHFR.inhibitor.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, DHFR.inhibitor.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "DHFR-inhibitors",
            .by = c(AnnotationAccession, SeqID, SeqType))

## combined because AAC(6')-Ib-cr5 and AAC(6')-Ib-cr7 are bifunctional acetyltransferases that confer resistance to both aminoglycosides and fluoroquinolones.
plasmid.aminoglycoside.and.quinolone.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, aminoglycoside.and.quinolone.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "aminoglycosides and quinolones",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.macrolide.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, macrolide.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "macrolides",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.antimicrobial.totals.df <- plasmid.ARGs |>
  filter(str_detect(product, antimicrobial.keywords)) |> 
  summarize(ARG_count = n(), ARG_type = "other antimicrobials",
            .by = c(AnnotationAccession, SeqID, SeqType))

plasmid.ARG.class.totals.df <- rbind(
  plasmid.chloramphenicol.totals.df,
  plasmid.tetracycline.totals.df,
  plasmid.lincosamide.totals.df,
  plasmid.multidrug.totals.df,
  plasmid.beta.lactam.totals.df,
  plasmid.glycopeptide.totals.df,
  plasmid.polypeptide.totals.df,
  plasmid.DHFR.inhibitor.totals.df,
  plasmid.aminoglycoside.and.quinolone.totals.df,
  plasmid.macrolide.totals.df,
  plasmid.antimicrobial.totals.df)

S1Fig <- PCN.data |>
  inner_join(plasmid.ARG.class.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~ARG_type)

## Facet by ecological annotation for each ARG class.
S2Fig <- PCN.data |>
  inner_join(plasmid.chloramphenicol.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Chloramphenicol resistance")

S3Fig <- PCN.data |>
  inner_join(plasmid.tetracycline.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Tetracycline resistance")

S4Fig <- PCN.data |>
  inner_join(plasmid.lincosamide.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Lincosamide resistance")

S5Fig <- PCN.data |>
  inner_join(plasmid.multidrug.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Multidrug resistance")

S6Fig <- PCN.data |>
  inner_join(plasmid.beta.lactam.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Beta-lactam resistance")

S7Fig <- PCN.data |>
  inner_join(plasmid.glycopeptide.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Glycopeptide resistance")

S8Fig <- PCN.data |>
  inner_join(plasmid.polypeptide.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Polypeptide resistance")

S9Fig <- PCN.data |>
  inner_join(plasmid.DHFR.inhibitor.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("DHFR inhibitor resistance")

S10Fig <- PCN.data |>
  inner_join(plasmid.aminoglycoside.and.quinolone.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Aminoglycoside and quinolone resistance")

S11Fig <- PCN.data |>
  inner_join(plasmid.macrolide.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Macrolide resistance")

S12Fig <- PCN.data |>
  inner_join(plasmid.antimicrobial.totals.df) |>
  mutate(has_ARG = TRUE) |> ## needed for make_ARG_PCN_base_plot()
  make_ARG_PCN_base_plot() +
  facet_wrap(.~Annotation) +
  ggtitle("Other antimicrobial resistance")

ggsave("../results/S1Fig.pdf", S1Fig, width=8)
ggsave("../results/S2Fig.pdf", S2Fig, width=8)
ggsave("../results/S3Fig.pdf", S3Fig, width=8)
ggsave("../results/S4Fig.pdf", S4Fig, width=8)
ggsave("../results/S5Fig.pdf", S5Fig, width=8)
ggsave("../results/S6Fig.pdf", S6Fig, width=8)
ggsave("../results/S7Fig.pdf", S7Fig, width=8)
ggsave("../results/S8Fig.pdf", S8Fig, width=8)
ggsave("../results/S9Fig.pdf", S9Fig, width=8)
ggsave("../results/S10Fig.pdf", S10Fig, width=8)
ggsave("../results/S11Fig.pdf", S11Fig, width=8)
ggsave("../results/S12Fig.pdf", S12Fig, width=8)


##################################################################
## This data frame is just to check internal consistency, not for plotting.
## check counts summed over individual ARG classes
## against counts summed over all ARGs to check for
## internal consistency.
plasmid.all.ARG.classes.totals.df <- plasmid.ARG.class.totals.df |>
  summarize(ARG_count = sum(ARG_count),
            .by = c(AnnotationAccession, SeqID, SeqType)) |>
  arrange(AnnotationAccession, SeqID, SeqType)
  
## assertion for internal consistency.
stopifnot(all.equal(plasmid.all.ARG.classes.totals.df, plasmid.ARG.totals.df))


##################################################################
## Figure 3. analyze the distribution of transposon genes on plasmids.
## Don't run on all MGE genes to show that this is not an artifact of
## counting plasmid function genes!
## For the sake of this paper, we mainly care about transposons.

plasmid.transposon.totals.df <- plasmid.MGEs |>
  filter(str_detect(product, transposon.keywords)) |>
  summarize(transposon_gene_count = n(),
            .by = c(AnnotationAccession, SeqID, SeqType))

transposon.annotated.PCN.data <- PCN.data |>
  left_join(plasmid.transposon.totals.df) |>
  ## set NAs to zero in the transposon_count column
  mutate(transposon_gene_count = replace_na(transposon_gene_count, 0)) |>
  mutate(has_transposon = ifelse(transposon_gene_count > 0, TRUE, FALSE))

## scatterplot of log10(Normalized plasmid copy number) vs. log10(plasmid length).
## using a somewhat unusual plotting strategy to avoid overplotting points with transposons
Fig3_base <- make_transposon_PCN_base_plot(transposon.annotated.PCN.data)

## add the marginal histogram to Figure A.
Fig3A <- ggExtra::ggMarginal(Fig3_base, margins="both")

## Figure 3B: facet by ecological annotation.
## By eye, looks like there may be some enrichment of transposons on small plasmids
## in humans but unclear-- need statistics to be rigorous.
Fig3B <- Fig3_base + guides(color = "none") + facet_wrap(.~Annotation)


## Figure 3CD: show differences in mobility type.
## IMPORTANT: report in figure legend that NA datapoints were removed.
Fig3C_base <- transposon.annotated.PCN.data |>
  filter(!is.na(PredictedMobility)) |>
  filter(has_transposon==FALSE) |>
  make_PCN_base_plot() + ggtitle("transposon- plasmids")

## Get the legend.
Fig3CD_legend <- get_legend(Fig3C_base)
## now remove the legend from base figure.
Fig3C_base <- Fig3C_base + guides(color="none")
Fig3C <- ggExtra::ggMarginal(Fig3C_base, groupColour = TRUE, groupFill = TRUE, margins="both")

Fig3D_base <- transposon.annotated.PCN.data |>
  filter(!is.na(PredictedMobility)) |>
  filter(has_transposon==TRUE) |>
  make_PCN_base_plot() + ggtitle("transposon+ plasmids") + guides(color="none")
Fig3D <- ggExtra::ggMarginal(Fig3D_base, groupColour = TRUE, groupFill = TRUE, margins="both")

Fig3CD <- plot_grid(plot_grid(Fig3C, Fig3D, labels = c("C","D"), nrow=1), Fig3CD_legend, nrow=2,rel_heights=c(1,0.05))

Fig3 <- plot_grid(plot_grid(Fig3A, Fig3B, labels = c("A", "B"), nrow = 1), Fig3CD, nrow = 2)

## Draft figure 3, showing that transposons are largely on large conjugative plasmids,
## and rarely on small plasmids (but this is observed).
ggsave("../results/Fig3.pdf", Fig3, width=7.5, height=8.5)

###########################################################################
## Given these findings, what functions ARE found on smaller plasmids <10kB in size?
## Let's take a look.
## Finding: significant linkage between ARGs and transposases on small plasmids.

small.plasmid.proteins <- plasmid.proteins |> filter(replicon_length < SIZE_THRESHOLD)

small.plasmid.function.summary <- small.plasmid.proteins |>
  summarize(function_count = n(), .by = c(product)) |>
  arrange(desc(function_count))

## there are 1336 ARGs found on small plasmids.
small.plasmid.ARGs <- small.plasmid.proteins |> 
  filter(str_detect(product, antibiotic.keywords))
nrow(small.plasmid.ARGs)

## what are these small plasmids with ARGs?
small.plasmids.with.ARGs <- PCN.data |>
  filter(SeqID %in% small.plasmid.ARGs$SeqID)
## There are 223 of these small plasmids with ARGs.
nrow(small.plasmids.with.ARGs)

## there are 3,364 small plasmids.
small.plasmids <- PCN.data |>
  filter(replicon_length < SIZE_THRESHOLD) |>
  as_tibble()
nrow(small.plasmids)

## what about transposons on these guys?
transposons.on.small.plasmids.with.ARGs <- plasmid.MGEs |>
  filter(str_detect(product, transposon.keywords)) |> 
  filter(SeqID %in% small.plasmid.ARGs$SeqID)
## THIS IS A NICE RESULT!
## there are 197 cases of transposons on these 223 small plasmids with ARGs.
nrow(transposons.on.small.plasmids.with.ARGs)

## now compare to baseline.
## There are 974 cases of transposases on small plasmids.
transposons.on.small.plasmids <- plasmid.MGEs |>
  filter(str_detect(product, transposon.keywords)) |> 
  filter(replicon_length < SIZE_THRESHOLD)
nrow(transposons.on.small.plasmids)

## 768 small plasmids have transposons, out of 3364 small plasmids.
length(unique(transposons.on.small.plasmids$SeqID))
length(unique(small.plasmids$SeqID))

## 161 out of 223 small plasmids with ARGs have transposons.
length(unique(transposons.on.small.plasmids.with.ARGs$SeqID))
length(unique(small.plasmids.with.ARGs$SeqID))

## This looks like a super significant result, showing an association
## between ARGs and transposons on small plasmids!!!

## Make a contingency table to test.
# ARG+/Tn+ ARG+/Tn-
# ARG-/Tn+ ARG-/Tn-
ARG.Tn.contingency.table <- matrix(
  c(161, 223-161,
    768-161, (3364 - 161 - (223-161) - (768 - 161))),
  nrow = 2,
  byrow = TRUE
)

## Odds ratio: 10.8
fisher.test(ARG.Tn.contingency.table)
## p < 1e-55
fisher.test(ARG.Tn.contingency.table)$p.value

###################################################################################
## Hypothesis: genomes containing small plasmids with ARGs are associated with bloodstream infections.
## Finding: Yes, highly significant association.

## 36 genomes with small plasmids with ARGs isolated from blood.
length(unique(filter(small.plasmids.with.ARGs, str_detect(isolation_source, "blood"))$AnnotationAccession))

##210 genomes have small plasmids with ARGs.
genomes.with.small.ARG.plasmids <- unique(small.plasmids.with.ARGs$AnnotationAccession)
length(genomes.with.small.ARG.plasmids)

## 298 genomes isolated from blood.
length(unique(filter(PCN.data, str_detect(isolation_source, "blood"))$AnnotationAccession))

## 4644 genomes total.
length(unique(PCN.data$AnnotationAccession))

## Make a contingency table to test.
# Blood+/small_ARG_plasmid+ Blood-/small_ARG_plasmid+
# Blood+/small_ARG_plasmid- Blood-/small_ARG_plasmid-
blood.small.ARG.plasmid.contingency.table <- matrix(
  c(36, 210-36,
    298-36, (4464 - 36 - (210-36) - (298-36))), 
  nrow = 2,
  byrow = TRUE
)

## Odds ratio: 3.15
fisher.test(blood.small.ARG.plasmid.contingency.table)
## p-value = 7.68e-08
fisher.test(blood.small.ARG.plasmid.contingency.table)$p.value

################################################################################
## Analyze duplicate pairs that are found on small plasmids and on big plasmids.
## import duplicate proteins, filter for transposons or transposons, and plasmid_count >= 2.
## then search in the filtered_plasmid_proteins to see if these occur on separate plasmids.

duplicated.proteins <- read.csv("../results/duplicate-proteins.csv")

duplicated.ARGs.only.on.plasmids <- duplicated.proteins |>
  filter(str_detect(product, antibiotic.keywords))

## how many of the small plasmid ARGs are duplicated?
## to estimate, filter duplicated.ARGs based on the (AnnotationAccession, product, sequence)
## rows in small.plasmid.ARGs.
small.plasmids.with.ARGs.filter.df <- small.plasmid.ARGs |>
  select(AnnotationAccession, product, sequence) |> distinct()

## there are 108 cases of duplicated ARGs with a copy on a small plasmid.
duplicated.ARGs.with.copy.on.small.plasmid <- duplicated.ARGs.only.on.plasmids |>
  semi_join(small.plasmids.with.ARGs.filter.df)
nrow(duplicated.ARGs.with.copy.on.small.plasmid)

duplicated.ARG.with.copy.on.small.plasmid.filter <- duplicated.ARGs.with.copy.on.small.plasmid |>
  select(AnnotationAccession, product, sequence) |> distinct()
  
## now, check to see if these are duplicated on the same plasmid,
## or duplicated across plasmids (say between large and small).

distinct.plasmid.filter.for.mothership.ARGs <- plasmid.ARGs |>
  semi_join(duplicated.ARG.with.copy.on.small.plasmid.filter) |>
  arrange(AnnotationAccession,product,replicon_length) |>
  ## get rid of ARGs that are duplicated on the same plasmid
  select(-SeqIndex) |>
  distinct() |>  ## get rid of duplicate ARGs on the same plasmid
  summarize(distinct_plasmids_with_ARG = n_distinct(SeqID),
            .by=c(AnnotationAccession, product, sequence, SeqType,
                  host, isolation_source, Annotation, Organism, Strain,
                  TaxonomicGroup, TaxonomicSubgroup, Genus)) |>
  ## filter for genomes with the ARG on multiple plasmids
  filter(distinct_plasmids_with_ARG > 1)

## 184 mothership ARGs
candidate.mothership.plasmid.ARGs <- plasmid.ARGs |>
  semi_join(distinct.plasmid.filter.for.mothership.ARGs) |>
  arrange(AnnotationAccession,product,replicon_length) |>
  ## manually remove a misannotated plasmid
  ## (same plasmid annotated as two separate plasmids)
  filter(AnnotationAccession != "GCF_001558295.2_ASM155829v2")
nrow(candidate.mothership.plasmid.ARGs)

## write to file to examine by eye.
write.csv(candidate.mothership.plasmid.ARGs,
          "../results/candidate-mothership-plasmid-ARGs.csv", quote=F, row.names=F)
## These look quite promising!
## This will be the basis for Figure 5.

## 63 genomes with candidate mothership ARGs
length(unique(candidate.mothership.plasmid.ARGs$AnnotationAccession))

## let's look at their data.
candidate.mothership.gbk.annotation <- gbk.annotation |>
  filter(AnnotationAccession %in% candidate.mothership.plasmid.ARGs$AnnotationAccession) |>
  arrange(AnnotationAccession, replicon_length)

## write to file to examine by eye.
write.csv(candidate.mothership.gbk.annotation,
          "../results/candidate-mothership-genomes.csv", quote=F, row.names=F)
## These look quite promising!

## Let's look at duplicated proteins in these candidate genomes.

duplicated.plasmid.proteins.in.candidate.mothership.genomes <- duplicated.proteins |>
  filter(plasmid_count >= 1) |> 
  filter(AnnotationAccession %in% candidate.mothership.gbk.annotation$AnnotationAccession) |>
  arrange(SeqIndex)

## write to file to examine by eye.
write.csv(duplicated.plasmid.proteins.in.candidate.mothership.genomes,
          "../results/duplicated-plasmid-proteins-in-candidate-mothership-genomes.csv", quote=F, row.names=F)

## look at the isolation sources for these genomes.
unique(candidate.mothership.gbk.annotation$isolation_source)

################################################################################
## Figure 4: Mathematical modeling figure using Pluto notebook in Julia.

################################################################################
## Figure 5: Bioinformatic evidence for mothership hypothesis:
## duplicated ARG-transposons found on small multicopy plasmids and on another plasmid as well.

## Design: an arc diagram. Each row is a genome; each point is a distinct
## plasmid (positioned by log10 length); a shallow arc connects two plasmids
## within a genome whenever they carry the same ARG (identical product
## annotation AND identical protein sequence). Multiple ARGs shared by the
## same plasmid pair are drawn as nested arcs.

## classify each ARG into a drug class, reusing the antibiotic keyword
## regexes already defined above so the classification stays consistent
## with the rest of the analysis.
classify_ARG_class <- function(product) {
  case_when(
    str_detect(product, chloramphenicol.keywords) ~ "Chloramphenicol",
    str_detect(product, tetracycline.keywords) ~ "Tetracycline",
    str_detect(product, lincosamide.keywords) ~ "Lincosamide",
    str_detect(product, beta.lactam.keywords) ~ "Beta-lactam",
    str_detect(product, glycopeptide.keywords) ~ "Glycopeptide",
    str_detect(product, polypeptide.keywords) ~ "Polypeptide",
    str_detect(product, DHFR.inhibitor.keywords) ~ "Sulfonamide/DHFR",
    str_detect(product, aminoglycoside.and.quinolone.keywords) ~ "Aminoglycoside/Quinolone",
    str_detect(product, macrolide.keywords) ~ "Macrolide",
    str_detect(product, multidrug.keywords) ~ "Multidrug",
    str_detect(product, antimicrobial.keywords) ~ "Other antimicrobial",
    TRUE ~ "Other"
  )
}


ARG_CLASS_COLORS <- c(
  "Beta-lactam" = "#d55e00",
  "Aminoglycoside/Quinolone" = "#0072b2",
  "Sulfonamide/DHFR" = "#009e73",
  "Chloramphenicol" = "#cc79a7",
  "Other antimicrobial" = "#999999")

## vertical offset between parallel lines connecting the same plasmid pair
## (subway-map style), so multiple shared ARGs stay visually distinct.
LINE_OFFSET_STEP <- 0.10

################################################################################
## get myLLannotator annotations for clinical genomes.
clinical.gbk.annotation <- read.csv("../results/labeled-gbk-annotation-table.csv") |>
  rename(ClinicalAnnotation = Annotation)

Fig5.data <- candidate.mothership.plasmid.ARGs |>
  ## Add genome labels: Organism + Strain is not always unique (some records have
  ## blank/NA Strain), so append the short RefSeq accession to guarantee a
  ## unique, traceable label for each genome.
  mutate(short_accession = str_extract(AnnotationAccession, "^GCF_[0-9.]+")) |>
  mutate(Sample = paste0(Organism, " (", short_accession, ")")) |>
  left_join(clinical.gbk.annotation) |>
  filter(ClinicalAnnotation == "Clinical")


## order genomes by the position of their smallest plasmid (the leftmost
## node), so genomes with the largest "smallest plasmid" are shown at the
## top of the figure and genomes with the smallest "smallest plasmid" are
## shown at the bottom.
genome.smallest.plasmid <- Fig5.data |>
  distinct(Sample, AnnotationAccession, SeqID, replicon_length) |>
  summarize(min_replicon_length = min(replicon_length), .by = Sample) |>
  arrange(min_replicon_length)

## order by smallest plasmid
sample.order <- genome.smallest.plasmid$Sample

ARG.plasmid.nodes <- Fig5.data |>
  distinct(Sample, AnnotationAccession, product, sequence, SeqID, replicon_length) |>
  mutate(Sample = factor(Sample, levels = sample.order)) |>
  mutate(log10_length = log10(replicon_length)) |>
  mutate(PlasmidSize = ifelse(replicon_length < SIZE_THRESHOLD, "<10 kb", ">=10 kb"))


## every plasmid-plasmid pair, within a genome, that shares an ARG
## (identical product annotation AND identical protein sequence).
Fig5.edges <- ARG.plasmid.nodes |>
  inner_join(ARG.plasmid.nodes,
             by = c("AnnotationAccession", "product", "sequence", "Sample"),
             suffix = c("1", "2"),
             relationship = "many-to-many") |>
  filter(SeqID1 < SeqID2) |> ## keep each unordered pair once
  mutate(ARG_class = classify_ARG_class(product)) |>
  ## offset multiple ARGs shared by the exact same plasmid pair into
  ## parallel lines, centered on the node's y-position (subway-map
  ## style), so the offsets stay
  ## symmetric and legible when several edges share a pair.
  mutate(stack_rank = row_number(), n_stack = n(),
         .by = c(AnnotationAccession, SeqID1, SeqID2)) |>
  mutate(y = as.numeric(Sample) + (stack_rank - (n_stack + 1) / 2) * LINE_OFFSET_STEP)


Fig5.nodes <- ARG.plasmid.nodes |>
  distinct(AnnotationAccession, SeqID, Sample, log10_length, PlasmidSize, replicon_length)


## keep y numeric throughout (nodes + arcs) so both layers share one
## continuous scale, then relabel the axis with genome names by hand.
sample.levels <- levels(Fig5.nodes$Sample)

Fig5_base <- ggplot() +
  geom_segment(data = Fig5.edges,
               aes(x = log10_length1, xend = log10_length2, y = y, yend = y,
                   color = ARG_class),
               linewidth = 0.6, alpha = 0.8, lineend = "round") +
  geom_point(data = Fig5.nodes,
             aes(x = log10_length, y = as.numeric(Sample), shape = PlasmidSize),
             size = 2.6, color = "grey20", fill = "white", stroke = 0.9) +
  scale_shape_manual(values = c(21, 19), name = "plasmid size") +
  scale_color_manual(values = ARG_CLASS_COLORS, name = "shared ARG class", na.value = "grey50") +
  scale_y_continuous(breaks = seq_along(sample.levels), labels = sample.levels,
                     expand = expansion(add = 0.6)) +
  theme_classic() +
  xlab("log10(plasmid length)") +
  ylab(NULL) +
  ggtitle("The mothership connection: amplification of ARGs\non small plasmids in clinical isolates") +
  theme(
    axis.text.y = element_text(size = 7),
    plot.title = element_text(size = 11),
    legend.position = "right") +
  guides(shape = "none")

## Get the legend.
Fig5_legend <- get_legend(Fig5_base)
## now remove the legend from base figure.
Fig5_base <- Fig5_base + guides(color="none")

Fig5 <- plot_grid(plot_grid(Fig5_legend, NULL,nrow=1), Fig5_base, nrow=2, rel_heights=c(0.2,0.8))
ggsave("../results/Fig5_base.pdf", Fig5, width = 7, height = 7)
