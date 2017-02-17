library(gridExtra)
library(ggplot2)
library(grid)

gg.colors <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                 "yellow", "#FF7F00", "red", "#7F0000"))

spectral.colors <- colorRampPalette(c('#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#ffffbf',
                                      '#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'))

# ternshells <- function(dt) {
#     xcart <- function(a, b, c) {
#         return((1 / 2) * ((2 * b + c) / (a + b + c)))
#     }
#     ycart <- function(a, b, c) {
#         return((sqrt(3) / 2) * (c / (a + b + c)))
#     }
# 
#     dt[,`:=`(An=round(A/(A+B+C+D), 2), Bn=round(B/(A+B+C+D), 2), Cn=round(C/(A+B+C+D), 2), Dn=round(D/(A+B+C+D), 2))]
#     t1<-dt[Dn==0]
#     t2<-dt[An==0 & Dn!=0]
#     t3<-dt[Cn==0 & An!=0]
#     t4<-dt[Bn==0 & Cn!=0]
# 
#     ternsh<-rbindlist(list(t1[,`:=`(xc=xcart(An, Bn, Cn), yc=ycart(An, Bn, Cn))],
#                            t2[,`:=`(xc=0.5+xcart(Cn, Dn, Bn), yc=abs(ycart(Cn, Dn, Bn)-sqrt(3)/2))],
#                            t3[,`:=`(xc=1+xcart(Bn, An, Dn), yc=ycart(Bn, An, Dn))],
#                            t4[,`:=`(xc=1.5+xcart(Dn, Cn, An), yc=abs(ycart(Dn, Cn, An)-sqrt(3)/2))]))
#     setkey(ternsh, Sample)
# 
#     return(ternsh)
# }

ternshells <- function(dt, ch=c('A', 'B', 'C', 'D')) {
    xcart <- function(a, b, c) {
        return((1 / 2) * ((2 * b + c) / (a + b + c)))
    }
    ycart <- function(a, b, c) {
        return((sqrt(3) / 2) * (c / (a + b + c)))
    }
    
    
    dt[,`:=`(An=round(get(ch[1])/(get(ch[1])+get(ch[2])+get(ch[3])+get(ch[4])), 2), 
             Bn=round(get(ch[2])/(get(ch[1])+get(ch[2])+get(ch[3])+get(ch[4])), 2), 
             Cn=round(get(ch[3])/(get(ch[1])+get(ch[2])+get(ch[3])+get(ch[4])), 2), 
             Dn=round(get(ch[4])/(get(ch[1])+get(ch[2])+get(ch[3])+get(ch[4])), 2))]
    t1<-dt[Dn==0]
    t2<-dt[An==0 & Dn!=0]
    t3<-dt[Cn==0 & An!=0]
    t4<-dt[Bn==0 & Cn!=0]
    
    ternsh<-rbindlist(list(t1[,`:=`(xc=xcart(An, Bn, Cn), yc=ycart(An, Bn, Cn))],
                           t2[,`:=`(xc=0.5+xcart(Cn, Dn, Bn), yc=abs(ycart(Cn, Dn, Bn)-sqrt(3)/2))],
                           t3[,`:=`(xc=1+xcart(Bn, An, Dn), yc=ycart(Bn, An, Dn))],
                           t4[,`:=`(xc=1.5+xcart(Dn, Cn, An), yc=abs(ycart(Dn, Cn, An)-sqrt(3)/2))]))
    setkey(ternsh, Sample)
    
    return(ternsh)
}

addpolycoords <- function(ternshdt) {
    splevels <- sort(unique(ternshdt$An))
    sp <- (splevels[2] - splevels[1]) / sqrt(3)
    xd = round(c(sp*cos(pi/6), sp*cos(pi/2), sp*cos(5*pi/6), sp*cos(7*pi/6), sp*cos(3*pi/2), sp*cos(11*pi/6)), 3)
    yd = round(c(sp*sin(pi/6), sp*sin(pi/2), sp*sin(5*pi/6), sp*sin(7*pi/6), sp*sin(3*pi/2), sp*sin(11*pi/6)), 3)

    polydf<-data.frame(xp=unlist(lapply(ternshdt$xc, function(x) x+xd)),
                       yp=unlist(lapply(ternshdt$yc, function(y) y+yd)),
                       Sample=rep(ternshdt$Sample, each=6),
                       id=rep(seq(1, length(ternshdt$Sample), 1), each=6))

    return(as.data.table(merge(polydf, ternshdt, by='Sample')))
}

tsplot <- function(ternshdt, fom, clabs=c('A', 'B', 'C', 'D'), fw=NULL, fg=NULL, symmetric=F, bounds=NULL) {

    polyshells<-addpolycoords(ternshdt)

    yoff <- sqrt(3)/16
    ymax <- sqrt(3)/2 + yoff
    ymin <- 0 - yoff

    labdf <- data.frame(x=seq(0, 2.5, 0.5),
                        y=rep(c(ymin, ymax), 3),
                        lab=clabs[c(1, 3, 2, 4, 1, 3)])

    ternlines <- data.frame(x=c(0, 1, 0.5, 0.5, 1.5, 1, 1, 2, 1.5, 1.5, 2.5, 2),
                            y=c(0, 0, sqrt(3)/2, sqrt(3)/2, sqrt(3)/2, 0, 0, 0, sqrt(3)/2, sqrt(3)/2, sqrt(3)/2, 0),
                            idl=rep(1:4, each=3))
    maxfom<-max(ternshdt[,get(fom)])
    minfom<-min(ternshdt[,get(fom)])
    if(is.null(bounds)){
        ubound<-abs(ifelse(abs(maxfom)>abs(minfom), maxfom, minfom))
        lbound<-ifelse(symmetric, ubound*-1, minfom)
    }
    else{
        ubound<-bounds[2]
        lbound<-bounds[1]
    }
    
    tsp <- ggplot(data=polyshells) +
        geom_polygon(aes(x=xp, y=yp, fill=get(fom), group=id)) +
        scale_fill_gradientn(fom, colours=rev(spectral.colors(20)), limit=c(lbound, ubound)) +
        coord_fixed() +
        geom_text(data=labdf, aes(x=x, y=y, label=lab), size=5) +
        geom_polygon(data=ternlines, aes(x=x, y=y, group=idl), color='black', fill=NA, size=0.4, alpha=0.5) +
        theme_bw() +
        theme(
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            legend.key.height=unit(0.02, 'snpc'),
            legend.key.width=unit(0.15, 'snpc'),
            legend.position='bottom',
            panel.border=element_blank(),
            panel.grid=element_blank(),
            text=element_text(size=14)
        )

    if(!is.null(fw)) {
        tsp <- tsp + facet_wrap(fw)
    }
    else {
        if(!is.null(fg)) {
            tsp <- tsp + facet_grid(fg)
        }
    }

    return(tsp)
}
