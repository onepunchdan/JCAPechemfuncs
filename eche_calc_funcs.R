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

calc.polyfit <- function(DT.avg, o=3, sweep.dir=0, interp.step=1E-3, drop.start=0, drop.end=0, models=F) {
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
    Dark.fit <- lm(DT.new[illum.bool==0, Iavg.A]~poly(Edark.V, degree=o, raw=T))
    Ill.fit <- lm(DT.new[illum.bool==1, Iavg.A]~poly(Eill.V, degree=o, raw=T))
    interp.seq <- seq(from=min(DT.new$Eavg.V), to=max(DT.new$Eavg.V), by=interp.step)
    DT.fit <- data.table(Efit.V=interp.seq, Iphoto.A=predict(Ill.fit, data.frame(Eill.V=interp.seq)) - 
                              predict(Dark.fit, data.frame(Edark.V=interp.seq)))
    dlist <- list(DT.avg=DT.new, DT.fit=na.omit(DT.fit))
    if (models) {
        dlist <- c(dlist, Dark.fit=Dark.fit, Ill.fit=Ill.fit)
    }
    return(dlist)
}

calc.loess <- function(DT.avg, o=2, sweep.dir=0, interp.step=1E-3, drop.start=0, drop.end=0, models=F) {
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
    Dark.fit <- loess(DT.new[illum.bool==0, Iavg.A]~Edark.V, degree=o)
    Ill.fit <- loess(DT.new[illum.bool==1, Iavg.A]~Eill.V, degree=o)
    interp.seq <- seq(from=min(DT.new$Eavg.V), to=max(DT.new$Eavg.V), by=interp.step)
    DT.fit <- data.table(Efit.V=interp.seq, Iphoto.A=predict(Ill.fit, data.frame(Eill.V=interp.seq)) - 
                                     predict(Dark.fit, data.frame(Edark.V=interp.seq)))
    dlist <- list(DT.avg=DT.new, DT.fit=na.omit(DT.fit))
    if (models) {
        dlist <- c(dlist, Dark.fit=Dark.fit, Ill.fit=Ill.fit)
    }
    return(dlist)
}

calc.fitfoms <- function(dlist, iphoto.min=1E-8) {
    root.inds <- which(dlist$DT.fit[,(Iphoto.A<=iphoto.min & shift(Iphoto.A, 1, type='lead', fill=F)>iphoto.min)])
    if(length(root.inds)==0) {
        
    }
    return(root.inds)
}

# Voc calculation needs Eo and redox couple; O2/H2O sort E.avg[illum.bool==1] by Iavg.A ascending; H+/H2 sort E.avg[illum.bool==1] by Iavg.A descending
# use head(E.avg[sorted], 2) for slope and x-intercept.

# cv.fitloess
# which.min((sort(cv.fitloess$DT.avg[illum.bool==1][Eavg.V>min(cv.fitloess$DT.fit$Efit.V) & Eavg.V<max(cv.fitloess$DT.fit$Efit.V)]$Eavg.V)[1]-cv.fitloess$DT.fit$Efit.V)^2)
# which.min((sort(cv.fitloess$DT.avg[illum.bool==1][Eavg.V>min(cv.fitloess$DT.fit$Efit.V) & Eavg.V<max(cv.fitloess$DT.fit$Efit.V)]$Eavg.V)[2]-cv.fitloess$DT.fit$Efit.V)^2)
# (cv.fitloess$DT.fit$Iphoto.A[63]-cv.fitloess$DT.fit$Iphoto.A[23])/(cv.fitloess$DT.fit$Efit.V[63]-cv.fitloess$DT.fit$Efit.V[23])
# 
# cv.fitpoly
# which.min((sort(cv.fitpoly$DT.avg[illum.bool==1][Eavg.V>min(cv.fitpoly$DT.fit$Efit.V) & Eavg.V<max(cv.fitpoly$DT.fit$Efit.V)]$Eavg.V)[1]-cv.fitpoly$DT.fit$Efit.V)^2)
# which.min((sort(cv.fitpoly$DT.avg[illum.bool==1][Eavg.V>min(cv.fitpoly$DT.fit$Efit.V) & Eavg.V<max(cv.fitpoly$DT.fit$Efit.V)]$Eavg.V)[2]-cv.fitpoly$DT.fit$Efit.V)^2)
# (cv.fitpoly$DT.fit$Iphoto.A[82]-cv.fitpoly$DT.fit$Iphoto.A[42])/(cv.fitpoly$DT.fit$Efit.V[82]-cv.fitpoly$DT.fit$Efit.V[42])
# 
# pm28<-getpm(28)
# pm49<-getpm(49)
# lm(sort(unique(pm49$x)) ~ poly(1:length(sort(unique(pm49$x))), degree=1, raw=T))

# 
# calc.foms <- function(DT.data, )
 
conv.xrf <- function(DT.int, calib='48-800-2000-vacu-3.2__20160811-v1b-medTa.csv') {
    transitions <- names(DT.int)[!names(DT.int) %in% c('BatchLabel', 'StgLabel', 'StagX', 'StagY', 'StagZ', 'Sample')]
    caldt <- fread(file.path(kd, 'experiments', 'xrfs', 'calibration_libraries', calib))
    DT.new <- copy(DT.int)
    for(t in transitions) {
        calind <- match(t, sub('\\.', '', caldt$transition))
        ellab <- paste0(strtrim(t, nchar(t)-1), '.nmol')
        nmols <- caldt[calind]$nmol.CPS * DT.new[, t, with=F]
        DT.new[, eval(ellab):=eval(nmols)]
    }
    return(DT.new)
}

smoothxrf<-function(intdt, nx, ny) {
    intnames<-names(intdt)[!names(intdt) %in% c('BatchLabel', 'StgLabel', 'StagX', 'StagY', 'StagZ', 'Sample')]
    melter<-function(intcol) {
        interplist<-interp(x=100-intdt$StagX, y=intdt$StagY, z=unlist(intdt[, intcol, with=F]), nx=nx, ny=ny)
        meltdt<-as.data.table(melt(interplist$z, na.rm=T))
        names(meltdt)<-c('x', 'y', intcol)
        meltdt[, `:=`(x=interplist$x[x], y=interplist$y[y]), with=F]
        setkey(meltdt, x, y)
        return(meltdt)
    }
    intlist<-lapply(intnames, melter)
    outdt<-intlist[[1]]
    for(i in 2:length(intlist)){
        outdt<-outdt[intlist[[i]]]
    }
    return(outdt)
}