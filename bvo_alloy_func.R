alloyplots<-function(fomdt, alloys=c('C', 'D', 'E', 'F', 'G', 'H'), titlestring='labelme') {
    pm<-getpm(72)
    pm[,`:=`(Bi=A*0.156, V=(A*0.144+B*0.11), C_=C*0.022, D_=D*0.022, E_=E*0.022, F_=F*0.022, G_=G*0.022, H_=H*0.022)]
    pm[,alloy_sum:=C_+D_+E_+F_+G_+H_]
    pm[,`:=`(BiV_ratio=round(Bi/V, 1))]
    setkey(fomdt, Sample)
    setkey(pm, Sample)

    # single alloys
    sal<-pm[fomdt][idx==2][A+B>0][D+E+F+G+H==0 | C+E+F+G+H==0 | C+D+F+G+H==0 | C+D+E+G+H==0 | C+D+E+F+H==0 | C+D+E+F+G==0]
    salm<-melt(sal[!is.na(pmax)][,.(pmed=median(pmax, na.rm=T), jmed=median(jsc, na.rm=T)), by=list(Bi, V, C_, D_, E_, F_, G_, H_, alloy_sum, BiV_ratio)], id.vars=c('Bi', 'V', 'pmed', 'jmed', 'alloy_sum', 'BiV_ratio'))
    salm[,variable:=factor(variable, levels=c('C_', 'D_', 'E_', 'F_', 'G_', 'H_'), labels=alloys)]

    # pures
    pur<-pm[fomdt][idx==2][A+B>0][C+D+E+F+G+H==0]
    purm<-melt(pur[!is.na(pmax)][,.(pmed=median(pmax, na.rm=T), jmed=median(jsc, na.rm=T)), by=list(Bi, V, C_, D_, E_, F_, G_, H_, alloy_sum, BiV_ratio)], id.vars=c('Bi', 'V', 'pmed', 'jmed', 'alloy_sum', 'BiV_ratio'))
    purm[,variable:=factor(variable, levels=c('C_', 'D_', 'E_', 'F_', 'G_', 'H_'), labels=alloys)]

    # merge pure and single alloys
    psal<-rbindlist(list(salm[value!=0], purm))
    psal[,alloy:='']

    # group into binaries

    cgrp<-pm[fomdt][idx==2][A+B>0][(C+D>0 & E+F+G+H==0) |
                                       (C+E>0 & D+F+G+H==0) |
                                       (C+F>0 & D+E+G+H==0) |
                                       (C+G>0 & D+E+F+H==0) |
                                       (C+H>0 & D+E+F+G==0)]

    dgrp<-pm[fomdt][idx==2][A+B>0][(D+C>0 & E+F+G+H==0) |
                                       (D+E>0 & C+F+G+H==0) |
                                       (D+F>0 & C+E+G+H==0) |
                                       (D+G>0 & C+E+F+H==0) |
                                       (D+H>0 & C+E+F+G==0)]

    egrp<-pm[fomdt][idx==2][A+B>0][(E+C>0 & D+F+G+H==0) |
                                       (E+D>0 & C+F+G+H==0) |
                                       (E+F>0 & C+D+G+H==0) |
                                       (E+G>0 & C+D+F+H==0) |
                                       (E+H>0 & C+D+F+G==0)]

    fgrp<-pm[fomdt][idx==2][A+B>0][(F+C>0 & D+E+G+H==0) |
                                       (F+D>0 & C+E+G+H==0) |
                                       (F+E>0 & C+D+G+H==0) |
                                       (F+G>0 & C+D+E+H==0) |
                                       (F+H>0 & C+D+E+G==0)]

    ggrp<-pm[fomdt][idx==2][A+B>0][(G+C>0 & D+E+F+H==0) |
                                       (G+D>0 & C+E+F+H==0) |
                                       (G+E>0 & C+D+F+H==0) |
                                       (G+F>0 & C+D+E+H==0) |
                                       (G+H>0 & C+D+E+F==0)]

    hgrp<-pm[fomdt][idx==2][A+B>0][(H+C>0 & D+E+F+G==0) |
                                       (H+D>0 & C+E+F+G==0) |
                                       (H+E>0 & C+D+F+G==0) |
                                       (H+F>0 & C+D+E+G==0) |
                                       (H+G>0 & C+D+E+F==0)]

    # calculate binary fractions per group
    cgrp[,alloy_frac:=C_/alloy_sum]
    dgrp[,alloy_frac:=D_/alloy_sum]
    egrp[,alloy_frac:=E_/alloy_sum]
    fgrp[,alloy_frac:=F_/alloy_sum]
    ggrp[,alloy_frac:=G_/alloy_sum]
    hgrp[,alloy_frac:=H_/alloy_sum]

    cmelt<-melt(cgrp[!is.na(pmax)][,.(pmed=median(pmax, na.rm=T), jmed=median(jsc, na.rm=T)), by=list(Bi, V, C_, D_, E_, F_, G_, H_, alloy_sum, alloy_frac, BiV_ratio)], id.vars=c('Bi', 'V', 'pmed', 'jmed', 'alloy_sum', 'alloy_frac', 'BiV_ratio'))
    dmelt<-melt(dgrp[!is.na(pmax)][,.(pmed=median(pmax, na.rm=T), jmed=median(jsc, na.rm=T)), by=list(Bi, V, C_, D_, E_, F_, G_, H_, alloy_sum, alloy_frac, BiV_ratio)], id.vars=c('Bi', 'V', 'pmed', 'jmed', 'alloy_sum', 'alloy_frac', 'BiV_ratio'))
    emelt<-melt(egrp[!is.na(pmax)][,.(pmed=median(pmax, na.rm=T), jmed=median(jsc, na.rm=T)), by=list(Bi, V, C_, D_, E_, F_, G_, H_, alloy_sum, alloy_frac, BiV_ratio)], id.vars=c('Bi', 'V', 'pmed', 'jmed', 'alloy_sum', 'alloy_frac', 'BiV_ratio'))
    fmelt<-melt(fgrp[!is.na(pmax)][,.(pmed=median(pmax, na.rm=T), jmed=median(jsc, na.rm=T)), by=list(Bi, V, C_, D_, E_, F_, G_, H_, alloy_sum, alloy_frac, BiV_ratio)], id.vars=c('Bi', 'V', 'pmed', 'jmed', 'alloy_sum', 'alloy_frac', 'BiV_ratio'))
    gmelt<-melt(ggrp[!is.na(pmax)][,.(pmed=median(pmax, na.rm=T), jmed=median(jsc, na.rm=T)), by=list(Bi, V, C_, D_, E_, F_, G_, H_, alloy_sum, alloy_frac, BiV_ratio)], id.vars=c('Bi', 'V', 'pmed', 'jmed', 'alloy_sum', 'alloy_frac', 'BiV_ratio'))
    hmelt<-melt(hgrp[!is.na(pmax)][,.(pmed=median(pmax, na.rm=T), jmed=median(jsc, na.rm=T)), by=list(Bi, V, C_, D_, E_, F_, G_, H_, alloy_sum, alloy_frac, BiV_ratio)], id.vars=c('Bi', 'V', 'pmed', 'jmed', 'alloy_sum', 'alloy_frac', 'BiV_ratio'))

    # apply alloy labels
    cmelt[,alloy:=alloys[1]]
    dmelt[,alloy:=alloys[2]]
    emelt[,alloy:=alloys[3]]
    fmelt[,alloy:=alloys[4]]
    gmelt[,alloy:=alloys[5]]
    hmelt[,alloy:=alloys[6]]

    # merge stuff
    all<-rbindlist(list(cmelt[variable!='C_'], dmelt[variable!='D_'], emelt[variable!='E_'], fmelt[variable!='F_'], gmelt[variable!='G_'], hmelt[variable!='H_']))
    all[,variable:=factor(variable, levels=c('C_', 'D_', 'E_', 'F_', 'G_', 'H_'), labels=alloys)]


    # next 3 blocks are for duplicating and labeling to satisfy ggplot facets
    # define alloy fraction for single alloys
    salm_c1<-copy(salm[value>0])
    salm_c1[, `:=`(alloy_frac=1, alloy=alloys[1])]
    salm_d1<-copy(salm[value>0])
    salm_d1[, `:=`(alloy_frac=1, alloy=alloys[2])]
    salm_e1<-copy(salm[value>0])
    salm_e1[, `:=`(alloy_frac=1, alloy=alloys[3])]
    salm_f1<-copy(salm[value>0])
    salm_f1[, `:=`(alloy_frac=1, alloy=alloys[4])]
    salm_g1<-copy(salm[value>0])
    salm_g1[, `:=`(alloy_frac=1, alloy=alloys[5])]
    salm_h1<-copy(salm[value>0])
    salm_h1[, `:=`(alloy_frac=1, alloy=alloys[6])]

    salm_1<-rbindlist(list(salm_c1[variable!=alloy], salm_d1[variable!=alloy], salm_e1[variable!=alloy], salm_f1[variable!=alloy], salm_g1[variable!=alloy], salm_h1[variable!=alloy]))

    # define alloy fraction for singe alloys
    salm_c0<-copy(salm[value>0])
    salm_c0[, `:=`(alloy_frac=0, alloy=alloys[1])]
    salm_d0<-copy(salm[value>0])
    salm_d0[, `:=`(alloy_frac=0, alloy=alloys[2])]
    salm_e0<-copy(salm[value>0])
    salm_e0[, `:=`(alloy_frac=0, alloy=alloys[3])]
    salm_f0<-copy(salm[value>0])
    salm_f0[, `:=`(alloy_frac=0, alloy=alloys[4])]
    salm_g0<-copy(salm[value>0])
    salm_g0[, `:=`(alloy_frac=0, alloy=alloys[5])]
    salm_h0<-copy(salm[value>0])
    salm_h0[, `:=`(alloy_frac=0, alloy=alloys[6])]

    salm_0<-rbindlist(list(salm_c0[variable==alloy], salm_d0[variable==alloy], salm_e0[variable==alloy], salm_f0[variable==alloy], salm_g0[variable==alloy], salm_h0[variable==alloy]))

    # overwrite variable names
    salm_c0b<-copy(salm_0)
    salm_c0b[,variable:=alloys[1]]
    salm_d0b<-copy(salm_0)
    salm_d0b[,variable:=alloys[2]]
    salm_e0b<-copy(salm_0)
    salm_e0b[,variable:=alloys[3]]
    salm_f0b<-copy(salm_0)
    salm_f0b[,variable:=alloys[4]]
    salm_g0b<-copy(salm_0)
    salm_g0b[,variable:=alloys[5]]
    salm_h0b<-copy(salm_0)
    salm_h0b[,variable:=alloys[6]]

    salm_0b<-rbindlist(list(salm_c0b[variable!=alloy], salm_d0b[variable!=alloy], salm_e0b[variable!=alloy], salm_f0b[variable!=alloy], salm_g0b[variable!=alloy], salm_h0b[variable!=alloy]))

    # merge binaries and single alloys
    bsal<-rbindlist(list(salm_1, salm_0b, all[,names(salm_1),with=F][value>0 & alloy_frac!=0]))

    # fix rounding stuffs
    bsal[,alloy_sum:=round(alloy_sum, 4)]
    
    # setup aesthetics and facet labels
    
    bsub<-bsal[alloy_sum %in% c(0.0022, 0.0044, 0.0066, 0.0110)]
    bsub[,sum_alloys:=factor(alloy_sum)]
    bsub[,Alab:=factor(alloy, levels=levels(factor(alloy)), labels=paste('A =', levels(factor(alloy))))]
    bsub[,vfrac:=round(V/(Bi+V),2)]
    bsub[,Yval:=floor(100*alloy_sum/(Bi+V+alloy_sum))/100]
    bsub[,vx_lab:=factor(vfrac, levels=levels(factor(vfrac)), labels=paste('x =', levels(factor(vfrac))))]
    psal[,vfrac:=round(V/(Bi+V),2)]
    psal[,Yval:=alloy_sum/(Bi+V+alloy_sum)]
    psal[,vx_lab:=factor(vfrac, levels=levels(factor(vfrac)), labels=paste('x =', levels(factor(vfrac))))]

    # make some plots
    splot<-ggplot(psal, aes(x=Yval, y=pmed*1E5, color=variable)) +
        theme_bw() +
        geom_line(aes(group=variable)) +
        geom_point() +
        facet_grid(alloy~vx_lab, switch='y') +
        xlab(expression(y~'in'~Bi[1-x]*V[x]*O[4+delta]:A[y])) +
        ylab(expression(P[max]~(mW/cm^2))) +
        scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0., 1)) +
        scale_color_manual(name=expression(A), values=gg.colors(6)) +
        theme(legend.position='right',
              strip.background=element_blank(),
              text=element_text(size=16),
              panel.spacing=unit(0.8, 'lines')) +
        ggtitle(titlestring)

    bplot<-ggplot(bsub, aes(x=alloy_frac, y=pmed*1E5, color=variable, shape=factor(Yval))) +
        theme_bw() +    
        geom_point() +
        geom_line(aes(group=interaction(variable, factor(floor(Yval*100))))) +
        facet_grid(Alab~vx_lab, switch='y') +
        xlab(expression(z~'in'~Bi[1-x]*V[x]*O[4+delta]:(A[1-z]*B[z])[y])) +
        ylab(expression(P[max]~(W))) +
        scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0., 1)) +
        scale_shape_manual(name=expression(y), values=1:4) +
        scale_color_manual(name=expression(B), values=gg.colors(6)) +
        theme(legend.position='right',
              strip.background=element_blank(),
              strip.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              text=element_text(size=16),
              panel.spacing=unit(0.8, 'lines'))

    sgrob<-ggplotGrob(splot)
    bgrob<-ggplotGrob(bplot)

    g<-rbind(sgrob, bgrob, size='first')
    g$widths<-unit.pmax(sgrob$widths, bgrob$widths)
    return(g)
}