## ARG-MGE-PCN-analysis.R by Rohan Maddamsetti.
## analyse the distribution of antibiotic resistance genes (ARGs)
## on chromosomes versus plasmids in fully-sequenced genomes and plasmids
## in the NCBI RefSeq database first analyzed in Maddamsetti et al. 2025.

library(tidyverse)
library(cowplot)
library(data.table)
library(ggExtra)

library(viridis) ## for viridis palettes
library(scico) ## in case these palettes useful: https://github.com/thomasp85/scico

################################################################################
## Regular expressions used in this analysis.

## We build on the regular expressions used by Zeevi et al. (2019).
## Transposon: ‘transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|IS[0-9]|conjugate transposon’
## plasmid: ‘relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation’
## phage: ‘capsid|phage|tail|head|tape measure|antiterminatio’
## other HGT mechanisms: ‘integrase|excision\S*|exo- nuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip’
## antibiotic resistance: ‘azole resistance|antibiotic resistance|TetR|tetracycline resistance|VanZ|betalactam\S*|beta-lactam|antimicrob\S*|lantibio\S*’.


## unknown protein keywords.
unknown.protein.keywords <- "unknown|Unknown|hypothetical|Hypothetical|Uncharacterized|Uncharacterised|uncharacterized|uncharacterised|DUF|unknow|putative protein in bacteria|Unassigned|unassigned"

## NOTE: some hypothetical proteins are "ISXX family insertion sequence hypothetical protein"
## so filter out those cases, when counting unknown proteins.

## "Tra" actually matches genes involved in conjugative transfer, which is better matched for plasmid functions rather than transposon functions!
## This is an oversight in the regular expressions used in the Zeevi et al. 2019 paper.

## match MGE genes using the following keywords in the "product" annotation
transposon.keywords <- "IS|transpos\\S*|insertion|Transpos\\S*|Tn[0-9]|tranposase|Tnp|Ins|ins"
plasmid.keywords <- "relax\\S*|conjug\\S*|Tra[A-Z]|Tra[0-9]|tra[A-Z]|mob\\S*|plasmid|chromosome partitioning|chromosome segregation|Mob\\S*|Plasmid|Rep|Conjug\\S*"
phage.keywords <- "capsid|phage|Tail|tail|head|tape measure|antiterminatio|Phage|virus|Baseplate|baseplate|coat|entry exclusion|Integrase|integrase"
other.HGT.keywords <- "excision\\S*|exonuclease|recomb|toxin|restrict\\S*|resolv\\S*|topoisomerase|reverse transcrip|intron|antitoxin|toxin|Toxin|Reverse transcriptase|hok|Hok|competence|addiction|type IV|conjugate transposon|post-segregation killing"


MGE.keywords <- paste(transposon.keywords, plasmid.keywords, phage.keywords, other.HGT.keywords, sep="|")
MGE.or.unknown.protein.keywords <- paste(MGE.keywords,unknown.protein.keywords,sep="|")


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


categorize.as.MGE.ARG.or.other <- function(product) {
  if (is.na(product))
    "Other function"
  else if (str_detect(product, antibiotic.keywords))
    "ARG"
  else if (str_detect(product, MGE.keywords))
    "MGE"
  else
    "Other function"
}


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
    scale_color_scico(palette = "batlow", name="MGE count") +
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


make_phage_PCN_base_plot <- function(phage.annotated.PCN.data) {
  ## Make the basic plot for Figure 3,
  ## before adding the marginal histograms or facetting
  phage.annotated.PCN.data |>
    ggplot() +
    geom_point(
      data = subset(phage.annotated.PCN.data, has_phage == FALSE),
      aes(x = log10_replicon_length, y = log10_PIRACopyNumber),
      color = "grey80",
      size = 0.5, alpha=0.8) +
    geom_point(
      data = subset(phage.annotated.PCN.data, has_phage == TRUE),
      aes(x = log10_replicon_length, y = log10_PIRACopyNumber, color = phage_gene_count),
      size = 0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    scale_color_scico(palette = "batlow", name="MGE count") +
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
## I can save a ton of memory if I don't import the sequence column,
## and by using the data.table package for import.
plasmid.proteins <- data.table::fread(
  "../results/filtered-plasmid-proteins.csv",
  drop="sequence") |>
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
## Figure 1. ARGs are concentrated in larger plasmids, and rarer in small plasmids.

ARG.annotated.PCN.data <- PCN.data |>
  left_join(plasmid.ARG.totals.df) |>
  ## set NAs to zero in the ARG_count column
  mutate(ARG_count = replace_na(ARG_count, 0)) |>
  ## Note: this is the minus sign character "−" U+2212 in Unicode.
  mutate(has_ARG = ifelse(ARG_count > 0, TRUE, FALSE))

## scatterplot of log10(Normalized plasmid copy number) vs. log10(plasmid length).
## using a somewhat unusual plotting strategy to avoid overplotting points with ARGs
Fig1_base <- make_ARG_PCN_base_plot(ARG.annotated.PCN.data)

## add the marginal histogram to Figure A.
Fig1A <- ggExtra::ggMarginal(Fig1_base, margins="both")

## Figure 1B: facet by ecological annotation.
## By eye, looks like there may be some enrichment of ARGs on small plasmids
## in humans but unclear-- need statistics to be rigorous.
Fig1B <- Fig1_base + guides(color = "none") + facet_wrap(.~Annotation)


## Figure 1CD: show differences in mobility type.
## IMPORTANT: report in figure legend that NA datapoints were removed.
Fig1C_base <- ARG.annotated.PCN.data |>
  filter(!is.na(PredictedMobility)) |>
  filter(has_ARG==FALSE) |>
  make_PCN_base_plot() + ggtitle("ARG− plasmids")

## Get the legend.
Fig1CD_legend <- get_legend(Fig1C_base)
## now remove the legend from base figure.
Fig1C_base <- Fig1C_base + guides(color="none")
Fig1C <- ggExtra::ggMarginal(Fig1C_base, groupColour = TRUE, groupFill = TRUE, margins="both")

Fig1D_base <- ARG.annotated.PCN.data |>
  filter(!is.na(PredictedMobility)) |>
  filter(has_ARG==TRUE) |>
  make_PCN_base_plot() + ggtitle("ARG+ plasmids") + guides(color="none")
Fig1D <- ggExtra::ggMarginal(Fig1D_base, groupColour = TRUE, groupFill = TRUE, margins="both")

Fig1CD <- plot_grid(plot_grid(Fig1C, Fig1D, nrow=1),Fig1CD_legend, nrow=2,rel_heights=c(1,0.05))

Fig1CD

Fig1 <- plot_grid(plot_grid(Fig1A, Fig1B, nrow = 1), Fig1CD, nrow = 2)

## Draft figure 1, showing that ARGs are largely on large conjugative plasmids,
## and rarely on small plasmids (but this is observed).
Fig1


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
## Figure 2. analyze the distribution of transposon genes on plasmids.
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
Fig2_base <- make_transposon_PCN_base_plot(transposon.annotated.PCN.data)

## add the marginal histogram to Figure A.
Fig2A <- ggExtra::ggMarginal(Fig2_base, margins="both")

## Figure 2B: facet by ecological annotation.
## By eye, looks like there may be some enrichment of transposons on small plasmids
## in humans but unclear-- need statistics to be rigorous.
Fig2B <- Fig2_base + guides(color = "none") + facet_wrap(.~Annotation)


## Figure 2CD: show differences in mobility type.
## IMPORTANT: report in figure legend that NA datapoints were removed.
Fig2C_base <- transposon.annotated.PCN.data |>
  filter(!is.na(PredictedMobility)) |>
  filter(has_transposon==FALSE) |>
  make_PCN_base_plot() + ggtitle("transposon− plasmids")

## Get the legend.
Fig2CD_legend <- get_legend(Fig2C_base)
## now remove the legend from base figure.
Fig2C_base <- Fig2C_base + guides(color="none")
Fig2C <- ggExtra::ggMarginal(Fig2C_base, groupColour = TRUE, groupFill = TRUE, margins="both")

Fig2D_base <- transposon.annotated.PCN.data |>
  filter(!is.na(PredictedMobility)) |>
  filter(has_transposon==TRUE) |>
  make_PCN_base_plot() + ggtitle("transposon+ plasmids") + guides(color="none")
Fig2D <- ggExtra::ggMarginal(Fig2D_base, groupColour = TRUE, groupFill = TRUE, margins="both")

Fig2CD <- plot_grid(plot_grid(Fig2C, Fig2D, nrow=1),Fig2CD_legend, nrow=2,rel_heights=c(1,0.05))

Fig2CD

Fig2 <- plot_grid(plot_grid(Fig2A, Fig2B, nrow = 1), Fig2CD, nrow = 2)

## Draft figure 2, showing that transposons are largely on large conjugative plasmids,
## and rarely on small plasmids (but this is observed).
Fig2

###########################################################################
## Given these findings, what functions ARE found on smaller plasmids <10kB in size?
SIZE_THRESHOLD <- 10000

small.plasmid.proteins <- plasmid.proteins |> filter(replicon_length < SIZE_THRESHOLD)

small.plasmid.function.summary <- small.plasmid.proteins |>
  summarize(function_count = n(), .by = c(product)) |>
  arrange(desc(function_count))

## there are 1336 ARGs found on small plasmids.
small.plasmid.ARGs <- small.plasmid.proteins |> 
  filter(str_detect(product, antibiotic.keywords))

## what are these small plasmids with ARGs?
small.plasmids.with.ARGs <- PCN.data |>
  filter(SeqID %in% small.plasmid.ARGs$SeqID)
## There are 223 of these small plasmids with ARGs.

## there are 3,364 small plasmids.
small.plasmids <- PCN.data |>
  filter(replicon_length < SIZE_THRESHOLD) |>
  as_tibble()


## what about transposons on these guys?
transposons.on.small.plasmids.with.ARGs <- plasmid.MGEs |>
  filter(str_detect(product, transposon.keywords)) |> 
  filter(SeqID %in% small.plasmid.ARGs$SeqID)
## THIS IS A SUPER PROMISING RESULT!
## there are 197 cases of transposons on these 223 small plasmids with ARGs.

## now compare to baseline.
## There are 974 cases of transposons on small plasmids.
transposons.on.small.plasmids <- plasmid.MGEs |>
  filter(str_detect(product, transposon.keywords)) |> 
  filter(replicon_length < SIZE_THRESHOLD)

## 768 small plasmids have transposons, out of 3364 plasmids.
length(unique(transposons.on.small.plasmids$SeqID))
length(unique(small.plasmids$SeqID))

## 161 out of 223 small plasmids with ARGs have transposons.
length(unique(transposons.on.small.plasmids.with.ARGs$SeqID))
length(unique(small.plasmids.with.ARGs$SeqID))

## This looks like a super significant result, showing an association
## between ARGs and transposons on small plasmids!!!
## TODO: double-check this, and make this argument more rigorous.


################################################################################
## Analyze duplicate pairs that are found on small plasmids and on big plasmids.
## import duplicate proteins, filter for transposons or transposons, and plasmid_count >= 2.
## then search in the filtered_plasmid_proteins to see if these occur on separate plasmids.

duplicated.proteins <- read.csv("../results/duplicate-proteins.csv")

duplicated.ARGs <- duplicated.proteins |>
  filter(str_detect(product, antibiotic.keywords))

## how many of the small plasmid ARGs are duplicated?
## to estimate, filter duplicated.ARGs based on the (AnnotationAccession, product)
## rows in small.plasmid.ARGs.
small.plasmids.with.ARGs.filter.df <- small.plasmid.ARGs |>
  select(AnnotationAccession, product) |> distinct()

## there are 109 cases of like this.
duplicated.ARGs.with.copy.on.small.plasmid <- duplicated.ARGs |>
  semi_join(small.plasmids.with.ARGs.filter.df)

duplicated.ARG.with.copy.on.small.plasmid.filter <- duplicated.ARGs.with.copy.on.small.plasmid |>
  select(AnnotationAccession, product, sequence) |> distinct()
  
## now, check to see if these are duplicated on the same plasmid,
## or duplicated across plasmids (say between large and small).

candidate.mothership.plasmid.ARGs <- plasmid.ARGs |>
  semi_join(duplicated.ARG.with.copy.on.small.plasmid.filter) |>
  ## get PCN data too.
  left_join(PCN.data) |> 
  arrange(AnnotationAccession,product,replicon_length)

## write to file to examine by eye.
write.csv(candidate.mothership.plasmid.ARGs,
          "../results/candidate-mothership-plasmid-ARGs.csv",
          quote=F, row.names=F)
## These look quite promising!

## 86 genomes.
unique(candidate.mothership.plasmid.ARGs$AnnotationAccession)

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

