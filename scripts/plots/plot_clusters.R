library(cowplot)
library(ggDoubleHeat)
library(ggplot2)
library(grid)

setwd("/Users/rashatwi/Desktop/fls/files/analysis/clusters/ps")

data<-read.table('0.25_li2s2.csv', sep=',', header=TRUE)
x_axis_labels <- 1:6
data$Solvent <- factor(data$Solvent, levels = rev(unique(data$Solvent)), ordered=TRUE)
p <- ggplot(data = data, aes(x = Size, y = Solvent)) +
    geom_heat_grid(outside = Second, inside = First, outside_colors=c('#FFFFFF', '#E7E8F8', '#C6C8F0', '#7D81DF', '#3035CB'), inside_colors=c("#FFFFFF", "#FBE7E8", "#F9DDDD", "#F0A1A2", "#D90F12"), outside_name="Second Shell", inside_name="First Shell", r=4) +
    ggtitle(expression(Li[2]~S[2])) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 0.5,
                                  frame.colour = "black",
                                  frame.linewidth = 0.5)) +
    theme(text = element_text(size=20), axis.text.y = element_text(color="black"), panel.background = element_rect(fill = NA), axis.text.x = element_text(angle=0, hjust=0.5, color="black"), panel.border = element_rect(colour = "black", fill=NA, size=2), legend.position = "none", axis.title = element_text(face="bold")) +
    scale_x_continuous(labels=x_axis_labels, breaks=x_axis_labels, expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))
ggsave("s2.pdf", p, width=2.8)

data<-read.table('0.25_li2s4.csv', sep=',', header=TRUE)
x_axis_labels <- 1:6
data$Solvent <- factor(data$Solvent, levels = rev(unique(data$Solvent)), ordered=TRUE)
p <- ggplot(data = data, aes(x = Size, y = Solvent)) +
    geom_heat_grid(outside = Second, inside = First, outside_colors=c('#FFFFFF', '#E7E8F8', '#C6C8F0', '#7D81DF', '#3035CB'), inside_colors=c("#FFFFFF", "#FBE7E8", "#F9DDDD", "#F0A1A2", "#D90F12"), outside_name="Second Shell", inside_name="First Shell", r=4) +
    ggtitle(expression(Li[2]~S[4])) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 0.5,
                                  frame.colour = "black",
                                  frame.linewidth = 0.5)) +
    theme(text = element_text(size=20), axis.text.y = element_text(color="black"), panel.background = element_rect(fill = NA), axis.text.x = element_text(angle=0, hjust=0.5, color="black"), panel.border = element_rect(colour = "black", fill=NA, size=2), legend.position = "none", axis.title = element_text(face="bold")) +
    scale_x_continuous(labels=x_axis_labels, breaks=x_axis_labels, expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))
ggsave("s4.pdf", p, width=2.8)

data<-read.table('0.25_li2s8.csv', sep=',', header=TRUE)
x_axis_labels <- 1:6
data$Solvent <- factor(data$Solvent, levels = rev(unique(data$Solvent)), ordered=TRUE)
p <- ggplot(data = data, aes(x = Size, y = Solvent)) +
    geom_heat_grid(outside = Second, inside = First, outside_colors=c('#FFFFFF', '#E7E8F8', '#C6C8F0', '#7D81DF', '#3035CB'), inside_colors=c("#FFFFFF", "#FBE7E8", "#F9DDDD", "#F0A1A2", "#D90F12"), outside_name="Second Shell", inside_name="First Shell", r=4) +
    ggtitle(expression(Li[2]~S[8])) +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 0.5,
                                  frame.colour = "black",
                                  frame.linewidth = 0.5)) +
    theme(text = element_text(size=20), axis.text.y = element_text(color="black"), panel.background = element_rect(fill = NA), axis.text.x = element_text(angle=0, hjust=0.5, color="black"), panel.border = element_rect(colour = "black", fill=NA, size=2), legend.position = "none", axis.title = element_text(face="bold")) +
    scale_x_continuous(labels=x_axis_labels, breaks=x_axis_labels, expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))
ggsave("s8.pdf", p, width=2.8)


library("tidyverse"); theme_set(theme_minimal())
theme_update(panel.grid.minor = element_blank())
library("ggdensity")
library("patchwork")
library(ggplot2)
  library(ggpubr)
  library(dplyr)
library(ggpubr)
library(cowplot)

data1<-read.table('/Users/rashatwi/Desktop/data.csv', sep=',', header=TRUE)
p1 <- ggplot(data1, aes(x=FDeg)) +
#     geom_hdr(xlim = c(0, 100), fill="#777DF2") +
    geom_point(aes(y=BE_Salt), shape = 21, size=2, fill="#777DF2") +
    facet_wrap(vars(BE_Salt)) +
#     geom_hdr(xlim = c(0, 100), fill="#777DF2") +
    geom_point(aes(y=IP), shape = 21, size=2, fill="#777DF2") +
    facet_wrap(vars(IP)) +
    theme(text = element_text(size=26), axis.text.y = element_text(color="black"), panel.background = element_rect(fill = NA), axis.text.x = element_text(angle=0, hjust=0.5, color="black"), panel.border = element_rect(colour = "black", fill=NA, size=2)) + ylim(-40, -10)


library("tidyverse"); theme_set(theme_minimal())
theme_update(panel.grid.minor = element_blank())
library("ggdensity")
library("patchwork")
library(ggplot2)
  library(ggpubr)
  library(dplyr)
library(ggpubr)
library(cowplot)
data<-read.table('/Users/rashatwi/Desktop/data.csv', sep=',', header=TRUE)
ggplot(data, aes(FDeg, value, fill = variable)) +
  geom_hdr(xlim = c(0, 100)) +
  geom_point(shape = 21) +
  facet_wrap(vars(variable))



library("tidyverse"); theme_set(theme_minimal())
theme_update(panel.grid.minor = element_blank())
library("ggdensity")
library("patchwork")
library(ggplot2)
  library(ggpubr)
  library(dplyr)
library(ggpubr)
library(cowplot)
library(ggpmisc)
library(ggrepel)
library(broom)
library(grid)

data1 <- read.table('/Users/rashatwi/Desktop/relative_be.csv', sep=',', header=TRUE)
p1 <- ggplot(data1, aes(BE_Salt_R, BE_PS_R)) +
  geom_point(shape = 21, size=6, fill="#777DF2") +
  geom_quadrant_lines() +
  stat_dens2d_filter(color = "red", keep.fraction = 0.02) +
  stat_dens2d_labels(aes(label = NA), keep.fraction = 0.02,
                     geom = "text_repel", size = 2, colour = "red") +
  theme(strip.background = element_blank(), text = element_text(size=50), axis.text.y = element_text(color="black"), panel.background = element_rect(fill = NA), axis.text.x = element_text(angle=0, hjust=0.5, color="black"), panel.border = element_rect(colour = "black", fill=NA, size=4), legend.position = "none", axis.title = element_text(face="bold", size=34), plot.tag = element_text(face="bold"), axis.ticks.length=unit(-0.5, "cm")) +
  scale_x_continuous(limits = c(-14,14), breaks = seq(-14,14,7)) +
  scale_y_continuous(limits = c(-14,14), breaks = seq(-14,14,7)) +
  labs(x = "Relative BE - LiTFSI (kcal/mol)", y = "Relative BE - LiPS (kcal/mol)", tag = "a)")

data2 <- read.table('/Users/rashatwi/Desktop/be.csv', sep=',', header=TRUE)
p2 <- ggplot(data2, aes(FDeg, value)) +
    geom_hdr(xlim = c(0, 100), fill="#777DF2") +
    geom_point(shape = 21, size=5, fill="#C7CAFA") +
    facet_wrap(vars(variable)) +
    theme(strip.background = element_blank(), text = element_text(size=50), axis.text.y = element_text(color="black"), panel.background = element_rect(fill = NA), axis.text.x = element_text(angle=0, hjust=0.5, color="black"), panel.border = element_rect(colour = "black", fill=NA, size=4), legend.position = "none", axis.title = element_text(face="bold", size=34), plot.tag = element_text(face="bold"), axis.ticks.length=unit(-0.5, "cm")) +
    ylim(-40, -10) +
    labs(x = "Flourination Degree (%)", y = "BE - LiTFSI (kcal/mol)", tag="b)")

data3 <- read.table('/Users/rashatwi/Desktop/ip.csv', sep=',', header=TRUE)
p3 <- ggplot(data3, aes(FDeg, value)) +
    geom_hdr(xlim = c(0, 100), fill="#929292") +
    geom_point(shape = 21, size=5, fill="#FFFFFF") +
    facet_wrap(vars(variable)) +
    theme(strip.background = element_blank(), text = element_text(size=50), axis.text.y = element_text(color="black"), panel.background = element_rect(fill = NA), axis.text.x = element_text(angle=0, hjust=0.5, color="black"), panel.border = element_rect(colour = "black", fill=NA, size=4), legend.position = "none", axis.title = element_text(face="bold", size=34), plot.tag = element_text(face="bold"), axis.ticks.length=unit(-0.5, "cm")) +
    ylim(4, 8) +
    labs(x = "Flourination Degree (%)", y = "IP (ev)", tag="c)")

p <- p1 + p2 + p3
ggsave("be.pdf", p + theme(panel.margin=unit(2,"cm")), width=25, height=10)