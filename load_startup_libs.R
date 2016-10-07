setwd('~')
library(png)
library(reshape2)
library(akima)
library(ggplot2)
library(grid)
library(gridExtra)

gg.colors <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

spectral.colors <- colorRampPalette(c('#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#ffffbf',
                                      '#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'))

# env
jd<-'J:/hte_jcap_app_proto'
kd<-'K:'
unztemp<-'E:/unztemp'
