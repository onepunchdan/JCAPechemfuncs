setwd("~")
library(reshape2)
library(ggplot2)
library(stringr)
library(grid)
library(data.table)
library(viridis)
# library(Cairo) library(png) library(akima) library(gridExtra) library(gtable)
# library(ggthemes) library(repr) library(devtools)
# install_github('nathanvan/parallelsugar') library(parallelsugar)


gg.colors <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
    "yellow", "#FF7F00", "red", "#7F0000"))

spectral.colors <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61",
    "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

# env
jd <- "J:/hte_jcap_app_proto"
kd <- "K:"
# unztemp<-ifelse(dir.exists('D:/TEMP'), 'D:/TEMP', 'E:/TEMP')


# cpng<-function(ggplt, fn=NULL, w=800, h=800) { if(is.null(fn)){
# fn=paste0(deparse(substitute(ggplt)),'.png') } CairoPNG(filename=fn, width=w,
# height=h) print(ggplt) dev.off() }

# gpng<-function(ggplt, fn=NULL, w=800, h=800) { if(is.null(fn)){
# fn=paste0(deparse(substitute(ggplt)),'.png') } CairoPNG(filename=fn, width=w,
# height=h) print(grid.newpage()) print(grid.draw(ggplt)) dev.off() }

# theme_fivethirtyeight<-function (base_size = 20, base_family = 'sans') {
# (theme_foundation(base_size = base_size, base_family = base_family) +
# theme(line = element_line(colour = 'black'), rect = element_rect(fill =
# ggthemes_data$fivethirtyeight['ltgray'], linetype = 0, colour = NA), text =
# element_text(colour = ggthemes_data$fivethirtyeight['dkgray']), # axis.title =
# element_blank(), axis.text = element_text(), axis.ticks = element_blank(),
# axis.line = element_blank(), legend.background = element_rect(),
# legend.position = 'bottom', legend.direction = 'horizontal', legend.box =
# 'vertical', panel.grid = element_line(colour = NULL), panel.grid.major =
# element_line(colour = ggthemes_data$fivethirtyeight['medgray']),
# panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0, size =
# rel(1.5), face = 'bold'), plot.margin = unit(c(1, 1, 1, 1), 'lines'),
# strip.background = element_rect())) }
