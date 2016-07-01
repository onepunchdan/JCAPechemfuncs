index.toggle <- function(DT, toggle.key, thresh, flip=F) {
    DT.new <- copy(DT)
    thresh.bools <- if(flip) DT.new[,get(toggle.key)]<=thresh else DT.new[,get(toggle.key)]>thresh
    thresh.inds <- numeric()
    thresh.inds <- 1
    foreach(i=2:nrow(DT.new)) %do% {
        thresh.inds <- c(thresh.inds, 
                         if(thresh.bools[i]==thresh.bools[i-1]) thresh.inds[i-1] else thresh.inds[i-1]+1)
    }
    DT.new[, `:=`(thresh.ind=thresh.inds, thresh.bool=thresh.bools)]
    DT.new[, group.ind:=1:.N, by=thresh.ind]
    DT.new[, group.frac:=group.ind/.N, by=thresh.ind]
    return(DT.new)
}

calc.avg <- function(DT.inds, omit.start=0.45, omit.end=0.9) {
    DT.new <- copy(DT.inds[group.frac>=omit.start & group.frac<=omit.end,
                          .(time.s=mean(t.s), illum.bool=mean(thresh.bool),
                            Eavg.V=mean(Ewe.V), Iavg.A=mean(I.A)), by=thresh.ind])
    return(DT.new)
}

calc.iphoto <- function(DT.avg, omit.start=0.45, omit.end=0.9) {
    DT.new <- copy(if('Iavg.A' %in% names(DT.avg)) DT.avg else calc.avg(DT.avg, omit.start, omit.end))
    DT.new[, Iphoto.A:=Iavg.A-((shift(Iavg.A, 1, type='lag') + shift(Iavg.A, 1, type='lead'))/2)]
    return(na.omit(DT.new)[illum.bool==1, .(time.s, Iphoto.A)])
}

calc.polyfit <- function(DT.avg, o=3, sweep.dir=0, interp.step=1E-3, drop.start=0, drop.end=0) {
    DT.new <- copy(DT.avg)
    DT.new[, sweep.bool:=Eavg.V > shift(Eavg.V, 1, type='lag')]
    DT.new[1, sweep.bool:=DT.new[2, sweep.bool]]
    if(sweep.dir!=0) {
        DT.new <- if(sweep.dir==1) DT.new[sweep.bool==T] else DT.new[sweep.bool==F]
    }
    DT.new <- DT.new[thresh.ind>=unique(thresh.ind)[drop.start+1] & 
                         thresh.ind<=rev(unique(thresh.ind))[drop.end+1]]
    # separate fits for dark and illum signal
    Edark.V <- DT.new[illum.bool==0, Eavg.V]
    Eill.V <- DT.new[illum.bool==1, Eavg.V]
    Dark.lm <- lm(DT.new[illum.bool==0, Iavg.A]~poly(Edark.V, degree=o, raw=T))
    Ill.lm <- lm(DT.new[illum.bool==1, Iavg.A]~poly(Eill.V, degree=o, raw=T))
    interp.seq <- seq(from=min(DT.new$Eavg.V), to=max(DT.new$Eavg.V), by=interp.step)
    DT.fit <- data.table(Efit.V=interp.seq, Iphoto.A=predict(Ill.lm, data.frame(Eill.V=interp.seq)) - 
                              predict(Dark.lm, data.frame(Edark.V=interp.seq)))
    return(list(DT.avg=DT.new, Dark.lm=Dark.lm, Ill.lm=Ill.lm, DT.fit=DT.fit))
}

calc.loess <- function(DT.avg, o=2, sweep.dir=0, interp.step=1E-3, drop.start=0, drop.end=0) {
    DT.new <- copy(DT.avg)
    DT.new[, sweep.bool:=Eavg.V > shift(Eavg.V, 1, type='lag')]
    DT.new[1, sweep.bool:=DT.new[2, sweep.bool]]
    if(sweep.dir!=0) {
        DT.new <- if(sweep.dir==1) DT.new[sweep.bool==T] else DT.new[sweep.bool==F]
    }
    DT.new <- DT.new[thresh.ind>=unique(thresh.ind)[drop.start+1] & 
                         thresh.ind<=rev(unique(thresh.ind))[drop.end+1]]
    # separate fits for dark and illum signal
    Edark.V <- DT.new[illum.bool==0, Eavg.V]
    Eill.V <- DT.new[illum.bool==1, Eavg.V]
    Dark.loess <- loess(DT.new[illum.bool==0, Iavg.A]~Edark.V, degree=o)
    Ill.loess <- loess(DT.new[illum.bool==1, Iavg.A]~Eill.V, degree=o)
    interp.seq <- seq(from=min(DT.new$Eavg.V), to=max(DT.new$Eavg.V), by=interp.step)
    Iill.A <- predict(Ill.loess, data.frame(Eill.V=interp.seq))
    Idark.A <- predict(Dark.loess, data.frame(Edark.V=interp.seq))
    DT.fit <- data.table(Efit.V=interp.seq, Iill.A, Idark.A)
    return(list(DT.avg=DT.new, Dark.loess=Dark.loess, Ill.loess=Ill.loess, DT.fit=na.omit(DT.fit)[,Iphoto.A:=Iill.A-Idark.A]))
}

# calc.fitfoms <- function(DT.fit, iphoto.min=1E-8) {
#     root.inds <- which(DT.fit[,(Iphoto.A<=iphoto.min & shift(Iphoto.A, 1, type='lead', fill=F)>iphoto.min)])
# }
