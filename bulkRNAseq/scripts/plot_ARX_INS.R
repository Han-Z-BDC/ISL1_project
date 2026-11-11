#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

d <- read.delim("arx_ins_tpm.tsv", stringsAsFactors = FALSE)
d$day <- factor(d$day, levels = c("D12","D15"))
d$genotype <- factor(d$genotype, levels = c("Ctr","KO"))
d$gene <- factor(d$gene, levels = c("ARX","INS"))

summary_df <- aggregate(TPM ~ day + genotype + gene, d, mean)
names(summary_df)[4] <- "mean_TPM"
sem <- aggregate(TPM ~ day + genotype + gene, d, function(x) sd(x)/sqrt(length(x)))
summary_df$sem <- sem$TPM

pal <- c(Ctr = "#1D5FB4", KO = "#C1502E")  # same activated/repressed-style pair used earlier in this project

p <- ggplot(summary_df, aes(x = day, y = mean_TPM, fill = genotype)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_TPM - sem, ymax = mean_TPM + sem),
                position = position_dodge(width = 0.7), width = 0.15, linewidth = 0.6) +
  geom_point(data = d, aes(x = day, y = TPM, group = genotype),
             position = position_dodge(width = 0.7), size = 2.2, shape = 21,
             fill = "white", color = "black", stroke = 0.6) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Expression (TPM)", fill = "Genotype",
       title = "ARX and INS expression, Ctr vs KO (GFP only, no rescue)") +
  theme_classic(base_size = 15, base_family = "Helvetica") +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(size = 13, hjust = 0.5),
        legend.position = "top")

ggsave("ARX_INS_expression.png", p, width = 7, height = 4.2, dpi = 220)
ggsave("ARX_INS_expression.pdf", p, width = 7, height = 4.2, device = cairo_pdf)
message("done")
