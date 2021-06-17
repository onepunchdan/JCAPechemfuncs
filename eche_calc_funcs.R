index.toggle <- function(DT, toggle.key, thresh, flip = F) {
    DT.new <- copy(DT)
    thresh.bools <- if (flip)
        DT.new[, get(toggle.key)] <= thresh else DT.new[, get(toggle.key)] > thresh
    thresh.inds <- numeric()
    thresh.inds <- 1
    for (i in 2:nrow(DT.new)) {
        thresh.inds <- c(thresh.inds, if (thresh.bools[i] == thresh.bools[i - 1]) thresh.inds[i -
            1] else thresh.inds[i - 1] + 1)
    }
    if (flip) {
        DT.new[, `:=`(thresh.bool, get(toggle.key) <= thresh)]
    } else {
        DT.new[, `:=`(thresh.bool, get(toggle.key) > thresh)]
    }
    # DT.new[, thresh.ind:=Eavg.V > shift(Eavg.V, 1, type='lag')]
    DT.new[, `:=`(thresh.ind = thresh.inds, thresh.bool = thresh.bools)]
    DT.new[, `:=`(thresh.ind = thresh.inds)]
    DT.new[, `:=`(group.ind, 1:.N), by = thresh.ind]
    DT.new[, `:=`(group.frac, group.ind/.N), by = thresh.ind]
    return(DT.new)
}

calc.avg <- function(DT.inds, omit.start = 0.45, omit.end = 0.9) {
    DT.new <- copy(DT.inds[group.frac >= omit.start & group.frac <= omit.end, .(time.s = mean(t.s),
        illum.bool = mean(thresh.bool), Eavg.V = mean(Ewe.V), Iavg.A = mean(I.A)),
        by = thresh.ind])
    return(DT.new)
}

calc.iphoto <- function(DT.avg, omit.start = 0.45, omit.end = 0.9) {
    DT.new <- copy(if ("Iavg.A" %in% names(DT.avg))
        DT.avg else calc.avg(DT.avg, omit.start, omit.end))
    DT.new[, `:=`(Iphoto.A, Iavg.A - ((shift(Iavg.A, 1, type = "lag") + shift(Iavg.A,
        1, type = "lead"))/2))]
    DT.new[, `:=`(Ephoto.V, Eavg.V)]
    return(na.omit(DT.new)[illum.bool == 1, .(time.s, Iphoto.A, Ephoto.V)])
}

calc.polyfit <- function(DT.avg, o = 3, sweep.dir = 0, interp.step = 0.001, drop.start = 0,
    drop.end = 0, models = F) {
    DT.new <- copy(DT.avg)
    DT.new[, `:=`(sweep.bool, Eavg.V > shift(Eavg.V, 1, type = "lag"))]
    DT.new[1, `:=`(sweep.bool, DT.new[2, sweep.bool])]
    if (sweep.dir != 0) {
        DT.new <- if (sweep.dir == 1)
            DT.new[sweep.bool == T] else DT.new[sweep.bool == F]
    }
    DT.new <- DT.new[thresh.ind >= unique(thresh.ind)[drop.start + 1] & thresh.ind <=
        rev(unique(thresh.ind))[drop.end + 1]]
    # separate fits for dark and illum signal
    Edark.V <- DT.new[illum.bool == 0, Eavg.V]
    Eill.V <- DT.new[illum.bool == 1, Eavg.V]
    Dark.fit <- lm(DT.new[illum.bool == 0, Iavg.A] ~ poly(Edark.V, degree = o, raw = T))
    Ill.fit <- lm(DT.new[illum.bool == 1, Iavg.A] ~ poly(Eill.V, degree = o, raw = T))
    interp.seq <- seq(from = min(DT.new$Eavg.V), to = max(DT.new$Eavg.V), by = interp.step)
    DT.fit <- data.table(Efit.V = interp.seq, Iphoto.A = predict(Ill.fit, data.frame(Eill.V = interp.seq)) -
        predict(Dark.fit, data.frame(Edark.V = interp.seq)))
    dlist <- list(DT.avg = DT.new, DT.fit = na.omit(DT.fit))
    if (models) {
        dlist <- c(dlist, Dark.fit = Dark.fit, Ill.fit = Ill.fit)
    }
    return(dlist)
}

calc.loess <- function(DT.avg, o = 2, sweep.dir = 0, interp.step = 0.001, drop.start = 0,
    drop.end = 0, models = F) {
    DT.new <- copy(DT.avg)
    DT.new[, `:=`(sweep.bool, Eavg.V > shift(Eavg.V, 1, type = "lag"))]
    DT.new[1, `:=`(sweep.bool, DT.new[2, sweep.bool])]
    if (sweep.dir != 0) {
        DT.new <- if (sweep.dir == 1)
            DT.new[sweep.bool == T] else DT.new[sweep.bool == F]
    }
    DT.new <- DT.new[thresh.ind >= unique(thresh.ind)[drop.start + 1] & thresh.ind <=
        rev(unique(thresh.ind))[drop.end + 1]]
    # separate fits for dark and illum signal
    Edark.V <- DT.new[illum.bool == 0, Eavg.V]
    Eill.V <- DT.new[illum.bool == 1, Eavg.V]
    Dark.fit <- loess(DT.new[illum.bool == 0, Iavg.A] ~ Edark.V, degree = o)
    Ill.fit <- loess(DT.new[illum.bool == 1, Iavg.A] ~ Eill.V, degree = o)
    interp.seq <- seq(from = min(DT.new$Eavg.V), to = max(DT.new$Eavg.V), by = interp.step)
    DT.fit <- data.table(Efit.V = interp.seq, Iphoto.A = predict(Ill.fit, data.frame(Eill.V = interp.seq)) -
        predict(Dark.fit, data.frame(Edark.V = interp.seq)))
    dlist <- list(DT.avg = DT.new, DT.fit = na.omit(DT.fit))
    if (models) {
        dlist <- c(dlist, Dark.fit = Dark.fit, Ill.fit = Ill.fit)
    }
    return(dlist)
}

calc.fitfoms <- function(dlist, iphoto.min = 1e-08) {
    root.inds <- which(dlist$DT.fit[, (Iphoto.A <= iphoto.min & shift(Iphoto.A, 1,
        type = "lead", fill = F) > iphoto.min)])
    if (length(root.inds) == 0) {
    }
    return(root.inds)
}



fpl.fit <- function(DT, weighy = F, plot = F) {
    # Self-starting...
    y <- DT$Iphoto.A
    x <- DT$Ephoto.V
    if (weighy) {
        log.ss <- nls(y ~ SSfpl(x, phi1, phi2, phi3, phi4), weights = 1/sqrt(y))
    } else {
        log.ss <- nls(y ~ SSfpl(x, phi1, phi2, phi3, phi4))
    }
    # C
    C_lower <- summary(log.ss)$coef[1]
    C_upper <- summary(log.ss)$coef[2]
    # a
    A <- summary(log.ss)$coef[3]
    # A <- exp((summary(log.ss)$coef[3]) * (1/summary(log.ss)$coef[4])) k
    K <- summary(log.ss)$coef[4]
    # K <- (1 / summary(log.ss)$coef[4])
    if (plot) {
        plot(y ~ x, main = "Four-Parameter Logistic Function", xlab = "Ephoto.V",
            ylab = "Iphoto.A")
        lines(seq(min(x), max(x), length.out = 1000), predict(log.ss, data.frame(x = seq(min(x),
            max(x), length.out = 1000))), col = "red")
    }
    r1 <- sum((x - mean(x))^2)
    r2 <- sum(residuals(log.ss)^2)
    r_sq <- (r1 - r2)/r1
    out <- data.frame(cbind(c(C_lower = C_lower, C_upper = C_upper, a = A, k = K,
        R.value = sqrt(r_sq))))
    names(out)[1] <- "Logistic Curve"
    # return(out)
    return(log.ss)
}

calc.logisfit <- function(DT.iph, ipb = 1e-08) {
    DT.new <- copy(DT.iph)
    pts <- nrow(DT.new)
    signs <- sapply(1:pts, function(x) ifelse(is.na(DT.new[x + 1]$Ephoto.V), sign(DT.new[x]$Ephoto.V -
        DT.new[x - 1]$Ephoto.V), sign(DT.new[x + 1]$Ephoto.V - DT.new[x]$Ephoto.V)))
    grp <- 1
    for (i in 2:pts) {
        ifelse(signs[i] == signs[i - 1], grp <- c(grp, grp[i - 1]), grp <- c(grp,
            grp[i - 1] + 1))
    }
    DT.new[, `:=`(sweep.dir, factor(signs, levels = c(-1, 1), labels = c("cathodic",
        "anodic")))]
    DT.new[, `:=`(sweep.idx, grp)]
    modlist <- lapply(unique(grp), function(x) tryCatch(fpl.fit(DT.new[sweep.idx ==
        x & Iphoto.A > 0]), error = function(e) NULL))
    names(modlist) <- unique(DT.new[, paste0(sweep.idx, "_", sweep.dir)])
    # return(modlist)
    calcpmax <- function(m) optimize(function(x) (eo - x) * predict(m, data.frame(x = x)),
        interval = c(0, 1), maximum = T)
    calcjsc <- function(m) predict(m, data.frame(x = eo))
    voclo
    vochi
    calcvoc <- function(m) optimize(function(x) (ipb - predict(m, data.frame(x = x)))^2,
        interval = c(0, 1))
    DT.fom <- rbindlist(lapply(1:length(modlist), function(i) {
        sweep.idx = as.numeric(strsplit(names(modlist)[i], "_")[[1]][1])
        sweep.dir = strsplit(names(modlist)[i], "_")[[1]][2]
        if (!is.null(modlist[[i]])) {
            opt.pmax = calcpmax(modlist[[i]])
            pred.jmp = predict(modlist[[i]], data.frame(x = opt.pmax$maximum))
            pred.jsc = calcjsc(modlist[[i]])
            opt.voc = calcvoc(modlist[[i]])
            return(list(idx = sweep.idx, dir = sweep.dir, pmax = opt.pmax$objective[1],
                vmp = eo - opt.pmax$maximum, jmp = pred.jmp[1], jsc = pred.jsc[1],
                voc = eo - opt.voc$minimum, ff = opt.pmax$objective[1]/(pred.jsc[1] *
                  (eo - opt.voc$minimum))))
        } else {
            return(list(idx = sweep.idx, dir = sweep.dir, pmax = NA, vmp = NA, jmp = NA,
                jsc = NA, voc = NA, ff = NA))
        }
    }))
    return(DT.fom)
}

calc.logisfit2 <- function(DT.iph, eo, ipb = 1e-08) {
    DT.new <- copy(DT.iph)
    pts <- nrow(DT.new)
    signs <- sapply(1:pts, function(x) ifelse(is.na(DT.new[x + 1]$Ephoto.V), sign(DT.new[x]$Ephoto.V -
        DT.new[x - 1]$Ephoto.V), sign(DT.new[x + 1]$Ephoto.V - DT.new[x]$Ephoto.V)))
    grp <- 1
    for (i in 2:pts) {
        ifelse(signs[i] == signs[i - 1], grp <- c(grp, grp[i - 1]), grp <- c(grp,
            grp[i - 1] + 1))
    }
    DT.new[, `:=`(sweep.dir, factor(signs, levels = c(-1, 1), labels = c("cathodic",
        "anodic")))]
    DT.new[, `:=`(sweep.idx, grp)]
    modlist <- lapply(unique(grp), function(x) tryCatch(fpl.fit(DT.new[sweep.idx ==
        x & Iphoto.A > 0]), error = function(e) NULL))
    names(modlist) <- unique(DT.new[, paste0(sweep.idx, "_", sweep.dir)])
    calcpmax <- function(m) optimize(function(x) (eo - x) * predict(m, data.frame(x = x)),
        interval = c(0, 1), maximum = T)
    calcjsc <- function(m) predict(m, data.frame(x = eo))
    calcvoc <- function(m) optimize(function(x) (ipb - predict(m, data.frame(x = x)))^2,
        interval = c(eo - 1, eo))
    DT.fom <- rbindlist(lapply(1:length(modlist), function(i) {
        sweep.idx = as.numeric(strsplit(names(modlist)[i], "_")[[1]][1])
        sweep.dir = strsplit(names(modlist)[i], "_")[[1]][2]
        if (!is.null(modlist[[i]])) {
            opt.pmax = calcpmax(modlist[[i]])
            pred.jmp = predict(modlist[[i]], data.frame(x = opt.pmax$maximum))
            pred.jsc = calcjsc(modlist[[i]])
            opt.voc = calcvoc(modlist[[i]])
            return(list(idx = sweep.idx, dir = sweep.dir, pmax = opt.pmax$objective[1],
                vmp = eo - opt.pmax$maximum, jmp = pred.jmp[1], jsc = pred.jsc[1],
                voc = eo - opt.voc$minimum, ff = opt.pmax$objective[1]/(pred.jsc[1] *
                  (eo - opt.voc$minimum))))
        } else {
            return(list(idx = sweep.idx, dir = sweep.dir, pmax = NA, vmp = NA, jmp = NA,
                jsc = NA, voc = NA, ff = NA))
        }
    }))
    return(list(fom = DT.fom, mod = modlist))
}

easylogis <- function(rawdt, eo, sweepidx = 2, interd = F, plot = F) {
    # rawdt<-readeche(file)[,.(t.s, Ewe.V, Ach.V, I.A, Toggle, sample_no)]
    snum <- if ("Sample" %in% names(rawdt))
        rawdt$Sample[1] else rawdt$sample_no[1]
    if ("t.s" %in% names(rawdt)) {
        idxdt <- index.toggle(rawdt, toggle.key = "Toggle", thresh = 0.5)
        iphdt <- calc.iphoto(idxdt)
        fit <- calc.logisfit2(iphdt, eo)
        fomdt <- fit$fom
        suppressWarnings({
            fomdt[, `:=`(Sample, snum)]
        })
        if (plot) {
            interval <- seq(round(min(rawdt$Ewe.V), 3), round(max(rawdt$Ewe.V), 3),
                0.001)
            fitplot <- ggplot() + geom_path(data = idxdt, aes(x = Ewe.V, y = I.A),
                color = "black") + geom_point(data = iphdt, aes(x = Ephoto.V, y = Iphoto.A),
                color = "blue") + geom_line(data = NULL, aes(x = interval, y = predict(fit$mod[[sweepidx]],
                data.frame(x = interval))), color = "red") + ggtitle(paste("Sample",
                snum)) + theme_bw()
            print(fitplot)
        }
    } else {
        idxdt <- NA
        iphdt <- NA
        fomdt <- data.table(idx = sweepidx, dir = NA, pmax = NA, vmp = NA, jmp = NA,
            jsc = NA, voc = NA, ff = NA, Sample = snum)
    }
    retlist <- list(raw = idxdt, iph = iphdt, fom = fomdt)
    retfom <- fomdt[idx == sweepidx]
    if (interd) {
        return(retlist)
    } else {
        return(retfom)
    }
}

conv.xrf <- function(DT.int, calib = "48-800-2000-vacu-3.2__20180726-v1-medTa.csv") {
    transitions <- names(DT.int)[!names(DT.int) %in% c("BatchLabel", "StgLabel",
        "StagX", "StagY", "StagZ", "Sample")]
    caldt <- fread(file.path(kd, "experiments", "xrfs", "user", "calibration_libraries",
        calib))
    DT.new <- copy(DT.int)
    for (t in transitions) {
        calind <- match(t, sub("\\.", "", caldt$transition))
        ellab <- paste0(strtrim(t, nchar(t) - 1), ".nmol")
        nmols <- caldt[calind]$nmol.CPS * DT.new[, t, with = F]
        DT.new[, `:=`(eval(ellab), eval(nmols))]
    }
    return(DT.new)
}
                                                        
calib.xrf <- function (DT.int, calib = "48-800-2000-vacu-3.2__20200123-v1-medTa.csv") 
{
    transitions <- grep('CPS', names(DT.int), value=T)
    caldt <- fread(file.path(kd, "experiments", "xrfs", "user", 
        "calibration_libraries", calib))
    DT.new <- copy(DT.int)
    for (t in transitions) {
        calind <- match(sub("\\.CPS", "", t), caldt$transition)
        ellab <- paste0(strtrim(t, nchar(t) - 4), ".nmol")
        nmols <- caldt[calind]$nmol.CPS * DT.new[, t, with = F]
        DT.new[, `:=`(eval(ellab), eval(nmols))]
    }
    return(DT.new)
}

smoothxrf <- function(intdt, nx, ny) {
    intnames <- names(intdt)[!names(intdt) %in% c("BatchLabel", "StgLabel", "StagX",
        "StagY", "StagZ", "Sample")]
    melter <- function(intcol) {
        interplist <- interp(x = 100 - intdt$StagX, y = intdt$StagY, z = unlist(intdt[,
            intcol, with = F]), nx = nx, ny = ny)
        meltdt <- as.data.table(melt(interplist$z, na.rm = T))
        names(meltdt) <- c("x", "y", intcol)
        meltdt[, `:=`(x = interplist$x[x], y = interplist$y[y]), with = F]
        setkey(meltdt, x, y)
        return(meltdt)
    }
    intlist <- lapply(intnames, melter)
    outdt <- intlist[[1]]
    for (i in 2:length(intlist)) {
        outdt <- outdt[intlist[[i]]]
    }
    return(outdt)
}

interpcomps <- function(xrfdt, fomdt, xrfpm, fompm) {
    setkey(xrfpm, Sample)
    setkey(fompm, Sample)
    setkey(xrfdt, Sample)
    setkey(fomdt, Sample)
    xrfsmps <- xrfpm[xrfdt]$Sample
    inx <- xrfpm[xrfdt]$x
    iny <- xrfpm[xrfdt]$y
    nmol_cols <- grep("\\.nmol", names(xrfdt), value = T)
    nmol_cols <- Filter(function(nc) all(!is.na(xrfdt[[nc]])), nmol_cols)
    fomsmps <- fompm[fomdt]$Sample
    outx <- sort(unique(fompm[fomdt]$x))
    outy <- sort(unique(fompm[fomdt]$y))
    outdt <- data.table(Sample = fompm[fomdt]$Sample, x = fompm[fomdt]$x, y = fompm[fomdt]$y)
    setkey(outdt, x, y)
    for (nc in nmol_cols) {
        ilist <- interp(x = inx, y = iny, z = xrfdt[[nc]], xo = outx, yo = outy,
            extrap = T, linear = F)
        idt <- as.data.table(interp2xyz(ilist, data.frame = T))
        setkey(idt, x, y)
        outdt[[nc]] <- idt[outdt]$z
    }
    return(outdt)
}
