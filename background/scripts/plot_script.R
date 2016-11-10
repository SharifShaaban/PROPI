#!/bin/Rscript
#n = 50, method.args = list(span = 2/3, degree = 1)

#############################################################################
# - Uses ggplot2 to graph trends.                                           #
# - Calculates Spearmann Rho rank test.                                     #
#############################################################################

library(readr)
library(magrittr)
library(ggplot2)
library(cowplot)

tSUM <- read_tsv("./temp/length_logs/summary_length_log.summary.tab", col_names = FALSE)
tLYS <- read_tsv("./temp/length_logs/LYS_length_log.summary.tab", col_names = FALSE)
tREC <- read_tsv("./temp/length_logs/REC_length_log.summary.tab", col_names = FALSE)
tREG <- read_tsv("./temp/length_logs/REG_length_log.summary.tab", col_names = FALSE)
tSTR <- read_tsv("./temp/length_logs/STR_length_log.summary.tab", col_names = FALSE)
tEFF <- read_tsv("./temp/length_logs/EFF_length_log.summary.tab", col_names = FALSE)
tMET <- read_tsv("./temp/length_logs/MET_length_log.summary.tab", col_names = FALSE)
tTRN <- read_tsv("./temp/length_logs/TRN_length_log.summary.tab", col_names = FALSE)
tOTH <- read_tsv("./temp/length_logs/OTH_length_log.summary.tab", col_names = FALSE)
tIS1 <- read_tsv("./temp/length_logs/IS_length_log.summary.tab", col_names = FALSE)
tIS2 <- read_tsv("./temp/length_logs/IS_length_log.summary2.tab", col_names = FALSE)

colnames(tSUM) <- c("name", "mbr_num", "avg_length")
colnames(tLYS) <- c("name", "avg_length", "avg_occ", "prop")
colnames(tREC) <- c("name", "avg_length", "avg_occ", "prop")
colnames(tREG) <- c("name", "avg_length", "avg_occ", "prop")
colnames(tSTR) <- c("name", "avg_length", "avg_occ", "prop")
colnames(tEFF) <- c("name", "avg_length", "avg_occ", "prop")
colnames(tMET) <- c("name", "avg_length", "avg_occ", "prop")
colnames(tTRN) <- c("name", "avg_length", "avg_occ", "prop")
colnames(tOTH) <- c("name", "avg_length", "avg_occ", "prop")
colnames(tIS1) <- c("name", "avg_length", "avg_occ", "prop")
colnames(tIS2) <- c("name", "avg_length", "avg_occ", "prop")

# pch = point symbols, col = point colors, lty = type de lignes, lwd = weight of lines


myTheme <- theme(plot.title = element_text(size=10), axis.text = element_text(size=10),
    axis.title.y = element_text(size=9, margin = margin(0,4,0,1, unit = "pt")),
    axis.title.x = element_text(size=9, margin = margin(1,0,4,0, unit = "pt")))

xmax <- 16
xbre <- 15
xbha <- 15
xlab <- "Number of members"
ylab <- "Average prophage length"
ymax <- 160000
titl <- "Average prophage length\nVS Number of cluster members"
cort <- cor.test(tSUM$mbr_num, tSUM$avg_length, method = "spearman")
colr <- "yellow"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pSUM <- ggplot(tSUM, aes(x = mbr_num, y = avg_length)) + stat_smooth(colour = colr, size = 0.5) + geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 12, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)


xmax <- 160000
xbre <- 150000
xbha <- 75000
xlab <- "Average prophage length"
ylab <- "Pseudo proportion"


ymax <- 55
titl <- "Lysis"
cort <- cor.test(tLYS$avg_length, tLYS$prop, method = "spearman")
colr <- "blue"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pLYS <- ggplot(tLYS, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)


ymax <- 55
titl <- "Recombination and replication"
cort <- cor.test(tREC$avg_length, tREC$prop, method = "spearman")
colr <- "cyan2"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pREC <- ggplot(tREC, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)


ymax <- 80
titl <- "Regulation"
cort <- cor.test(tREG$avg_length, tREG$prop, method = "spearman")
colr <- "chocolate"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pREG <- ggplot(tREG, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)


ymax <- 150
titl <- "Structure"
cort <- cor.test(tSTR$avg_length, tSTR$prop, method = "spearman")
colr <- "coral"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pSTR <- ggplot(tSTR, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)


ymax <- 80
titl <- "Effector and virulence factor"
cort <- cor.test(tEFF$avg_length, tEFF$prop, method = "spearman")
colr <- "chartreuse"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pEFF <- ggplot(tEFF, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)


ymax <- 25
titl <- "Metabolism and transport"
cort <- cor.test(tMET$avg_length, tMET$prop, method = "spearman")
colr <- "burlywood"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pMET <- ggplot(tMET, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)


ymax <- 25
titl <- "tRNA"
cort <- cor.test(tTRN$avg_length, tTRN$prop, method = "spearman")
colr <- "cadetblue"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pTRN <- ggplot(tTRN, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)


ymax <- 150
titl <- "Other or ambiguous"
cort <- cor.test(tOTH$avg_length, tOTH$prop, method = "spearman")
colr <- "brown"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pOTH <- ggplot(tOTH, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)

ymax <- 25
titl <- "Insertion sequence (All)"
cort <- cor.test(tIS1$avg_length, tIS1$prop, method = "spearman")
colr <- "red"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pIS1 <- ggplot(tIS1, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)


ymax <- 25
titl <- "Insertion sequence (Average > 0)"
cort <- cor.test(tIS2$avg_length, tIS2$prop, method = "spearman")
colr <- "black"
rhov <- format(round(cort$estimate, digits = 2), nsmall = 2)
pval <- format(ceiling(cort$p.value * 100) / 100, nsmall = 2)
if ( grepl("-",rhov) == FALSE) {
    rhov = paste0("+", rhov)
} 

pIS2 <- ggplot(tIS2, aes(x = avg_length, y = prop)) + stat_smooth(colour = colr, size = 0.5) + 
    geom_point(colour = colr, size = 0.2) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    scale_y_continuous(breaks=seq(0, ymax, ymax)) + scale_x_continuous(breaks=seq(0, xbre, xbha)) +
    labs(list(title = titl, x = xlab, y = ylab)) +
    myTheme + annotate("text", x = 125000, y = ymax,
    label = paste0("rho = ", rhov,
        "\npval \u2264 ", pval), family = "mono", size = 3, hjust = 0.4, vjust = 1)

my_plots <- plot_grid(
    pSUM,
    pREG,
    pEFF,
    pMET,
    pTRN,
    pIS1,
    pIS2,
    pREC,
    pLYS,
    pSTR,
    pOTH,
    labels = LETTERS[1:11], ncol = 3
)

#save_plot("gene_freq_plot_ggplot2.pdf", my_plots, base_height = 8, base_aspect_ratio = 1.0, device = cairo_pdf)

ggsave(file="./output/gene_freq_plot_ggplot2.svg", plot=my_plots, height=8, width = 8)
