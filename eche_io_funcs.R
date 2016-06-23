library(data.table)
library(stringr)
library(readr)
library(parallel)

readeche <- function(path) {
    dt <- fread(path)
    f <- basename(path)
    # use R-friendly names
    names(dt) <- sub('%column_headings=', '', names(dt))
    names(dt) <- sub('\\)', '', names(dt))
    names(dt) <- sub('\\(', '.', names(dt))
    # assign filename-parsed metadata to data table
    dt[, `:=`(sample_no=as.integer(str_match(f, '^Sample([0-9]+)_')[,2]),
              technique=str_match(f, '_([A-Z]+)[0-9]+\\.txt')[,2],
              technum=as.integer(str_match(f, '_[A-Z]+([0-9]+)\\.txt')[,2]),
              run=basename(sub(f, '', path))
              )]
    return(dt)
}

readopt <- function(path) {
    dt <- fread(path)
    f <- basename(path)
    # use R-friendly names
    names(dt) <- sub('\\)', '', names(dt))
    names(dt) <- sub('\\(', '_', names(dt))
    # assign filename-parsed metadata to data table
    dt[, `:=`(sample_no=as.integer(str_match(f, '^Sample([0-9]+)_')[,2]),
              technique=str_match(f, '_([A-Z]+)[0-9]+_(TRANS|DARK)')[,2],
              technum=as.integer(str_match(f, '_[A-Z]+([0-9]+)_(TRANS|DARK)')[,2]),
              trans=as.integer(str_match(f, '_(TRANS|DARK)([0-9]+)_')[,3]),
              time=as.character(strptime(str_match(f, '_([0-9]{8}\\.[0-9]{6})\\.[0-9]\\.opt')[,2], format='%Y%m%d.%H%M%S')),
              run=basename(sub(f, '', path))
              )]
    return(dt)
}

readechefiles <- function(l) {
    return(rbindlist(lapply(l, readeche)))
}

library(doParallel)
readechefiles2 <- function(l) {
  return(rbindlist(foreach(i=1:length(l), .export=c('readeche', 'fread', 'str_match')) %dopar% readeche(l[i])))
}

readechedir <- function(dir) {
  return(readechefiles(list.files(dir, pattern='.txt', full.names=T)))
}

readoptfiles <- function(l) {
    return(rbindlist(lapply(l, readopt)))
}

readoptdir <- function(dir) {
  return(readoptfiles(list.files(dir, pattern='.opt', full.names=T)))
}

readxrf <- function(f) {
    orbcsv<-fread(f, skip=6, drop=c('StagR', 'MatchName', 'MatchSig', 'CoatThk1', 'CoatThk2','CoatConc'))

    # check stage labels for sample number, must be unique, if not, generate 1:nrow(df)
    if(is.numeric(orbcsv$StgLabel)){
        orbcsv[,Sample:=StgLabel]
    }
    else{
        orbcsv[,Sample:=gsub("^\\s+|\\s+$", "", StgLabel)]
    }

    setkey(orbcsv, 'Sample')
    setnames(orbcsv, sub('\\s+', '', names(orbcsv)))

    # get column indices
    indInt<-which(names(orbcsv)=='Inte')
    indWt<-which(names(orbcsv)=='Wt%')
    indAt<-which(names(orbcsv)=='At%')
    indLbl<-which(names(orbcsv)=='StgLabel')
    extracol<-indLbl+4
    keepInt<-c(indWt:(indLbl-1), extracol)
    keepWt<-c((indInt+1):(indWt), indAt:(indLbl-1), extracol)
    keepAt<-c((indInt+1):(indAt), extracol)

    setnames(orbcsv, 'Inte', 'BatchLabel')

    return(list('Int'=orbcsv[, -keepInt, with=F],
                'Wt'=orbcsv[, -keepWt, with=F],
                'At'=orbcsv[, -keepAt, with=F]))
}

readpm <- function(f) {
    pmdt <- fread(f, skip=2)
    namesvec <- unlist(strsplit(readLines(f, n = 2)[2], ', '))
    namesvec <- namesvec[1:grep('code', namesvec)]
    setnames(pmdt, namesvec)
    setnames(pmdt, grep('Sample', names(pmdt), value = T), 'Sample')
    setnames(pmdt, grep('code', names(pmdt), value = T), 'code')
    setnames(pmdt, sub('\\(.*\\)', '', names(pmdt)))
    setkey(pmdt, 'Sample')
    return(pmdt)
}

getpm <- function(pmnum) {
    pmfiles <- list.files(file.path(jd, 'map'), full=T, pattern=paste0('^0*', pmnum, '-'))
    pmtxts <- grep('\\.txt', pmfiles, value=T)
    matches <- length(pmtxts)
    if(matches == 0) {
        print(paste('Platemap', pmnum, 'not found.'))
        return(NA)
    }
    if(matches > 1) {
        print(paste('Found', matches, 'ascii files for platemap', pmnum))
        print('Returning first match.')
    }
    else {
        return(readpm(pmtxts[1]))
    }
}

findnn <- function(pm, targetcode, nnd, channel='A', spancodes=F, sym=T) {
    pmdt <- readpm(pm)
    # determine nearest-neighbor distance by code:
    ddt <- pmdt[code %% 10 == 0,
                .(atpct = sort(unique(eval(parse(text = channel))))),
                by = code][,.(atdist = round(sqrt(2*(atpct[nnd + 1] - atpct[1])^2), 4)), by = code]
    # span like-codes for averaging or restrict to single code? create subtable with filtered codes
    if (spancodes) {
        subdt <- pmdt[code %in% ddt[atdist == ddt[code == targetcode, atdist], code]]
    }
    else {
        subdt <- pmdt[code == targetcode]
    }
    # per sample: return list of subtable indicies with composition distance < nnd & check symmetry
    if (sym) {
        subdt[, nninds := subdt[, .(list(intersect(
            which(round(sqrt((A - subdt[, A])^2 +
                                 (B - subdt[, B])^2 +
                                 (C - subdt[, C])^2 +
                                 (D - subdt[, D])^2), 4) <= ddt[code == targetcode, atdist]),
            which((if(A == min(subdt[, A])) subdt[, A] == min(subdt[, A]) else rep(T, nrow(subdt))) &
                      (if(B == min(subdt[, B])) subdt[, B] == min(subdt[, B]) else rep(T, nrow(subdt))) &
                      (if(C == min(subdt[, C])) subdt[, C] == min(subdt[, C]) else rep(T, nrow(subdt))) &
                      (if(D == min(subdt[, D])) subdt[, D] == min(subdt[, D]) else rep(T, nrow(subdt))))
            ))), by = Sample]$V1]
    }
    else {
        subdt[, nninds := subdt[, .(list(
            which(round(sqrt((A - subdt[, A])^2 +
                                 (B - subdt[, B])^2 +
                                 (C - subdt[, C])^2 +
                                 (D - subdt[, D])^2), 4) <= ddt[code == targetcode, atdist])
            )), by = Sample]$V1]
    }
    subdt[, nnsmps := lapply(nninds, function(x) Sample[unlist(x)])]
    return(subdt)
}

mean_omit_outlier <- function(x, sigma) {
    if (length(x)<=1) {
        return(x)
    }
    else{
        while (any(x > mean(x) + sigma * sd(x)) |
               any(x < mean(x) - sigma * sd(x))) {
            x <- x[(x <= mean(x) + sigma * sd(x)) & (x >= mean(x) - sigma * sd(x))]
        }
        return(mean(x))
    }
}

smoothfom <- function(fom, smoothdt, sigma) {
    if(is.character(fom)) {
        fdt <- fread(fom)
    }
    else {
        fdt <- fom
    }
    smpcol_name <- names(fdt)[names(fdt) %in% c('Sample', 'sample', 'sample_no')]
    setnames(fdt, smpcol_name, 'Sample')
    fomnames <- names(fdt)[!names(fdt)=='Sample']
    setkey(fdt, Sample)
    setkey(smoothdt, Sample)
    subdt <- fdt[smoothdt]
    subdt[, nnsmps_filt:=lapply(nnsmps, function(x) intersect(Sample, unlist(x)))]
    for(fom in fomnames){
        subdt[, paste0(fom, '_nn', sigma, 'sig') := sapply(1:nrow(subdt), function(n) subdt[Sample %in% unlist(nnsmps_filt[n]) & !is.na(get(fom)), mean_omit_outlier(get(fom), sigma)])]
    }
    return(subdt)
}


writefom <- function(df) {
    outdf<-as.data.frame(df)
    outdf[,!(names(outdf) %in% c('Sample', 'sample_no', 'Label', 'StgLabel')) & apply(outdf, 2, is.numeric)]<-format(outdf[,!(names(outdf) %in% c('Sample', 'sample_no', 'Label', 'StgLabel')) & apply(outdf, 2, is.numeric)], scientific=T, digits=3)
    write.table(outdf, file=paste0(deparse(substitute(df)),'.csv'),
                sep=',', row.names=FALSE, na="NaN", quote=FALSE) #,eol="\r\n",
}



## NON-DT functions
##

readana<-function(anapath){
    frl = function(fname) { # faster readlines
        s = file.info( fname )$size
        buf = readChar( fname, s, useBytes=T)
        strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
    }
    stxt <- frl(anapath)
    # analist<-readRCP(anapath)

    dflist<-list()
    anapathsplit<-unlist(strsplit(anapath, '/'))
    anafile<-anapathsplit[length(anapathsplit)]
    anafolder<-sub(anafile, '', anapath)
    # anas<-grep('ana__', names(analist), value=T)
    # for(i in anas){
    # fomfiles<-names(analist[[i]][['files_multi_run']][['fom_files']])
    csvlines <- unique(grep('.csv', stxt, fixed=T, value=T))
    fomfiles <- sapply(csvlines, function(x) trim(strsplit(x, ':', fixed=T)[[1]][1]))
    # print(as.character(fomfiles))
    for(f in fomfiles){
        ananame<-paste0('ana', '__', strsplit(f, '__')[[1]][2])
        fompath<-paste0(anafolder, f)
        fomcon<-file(fompath)
        fomhead<-as.numeric(unlist(strsplit(readLines(fomcon, n=1), '\t')))
        close(fomcon)
        fomtbl<-read.table(fompath, header=T, sep=',', na.strings='NaN', skip=fomhead[4]+1, nrows=fomhead[3]+1)
        # names(fomtbl)<-paste0('ana', strsplit(f, '__')[[1]][2], '_', names(fomtbl))
        # names(fomtbl)<-paste0(sub('__', '', i), '_', names(fomtbl))
        dflist[[ananame]]<-I(as.data.frame(fomtbl))
    }
    # }
    return(dflist)
}

# Takes a single .rcp file path and returns a list with nested lists for params,
# and filenames.
readnest <- function(fpath, ext, unzip=F) {
    if(unzip) {
        fn<-strsplit(fpath, '.', fixed=T)[[1]]
        rcpfile<-paste0(tail(strsplit(fn[1], '/')[[1]], n=1), '.', fn[2], ext)
        filecon<-unz(fpath, filename=rcpfile)
        filezcon<-gzcon(filecon)
        stxt <- readLines(filezcon)
        close(filezcon)
    }
    else {
        frl = function(fname) { # faster readlines
            s = file.info( fname )$size
            buf = readChar( fname, s, useBytes=T)
            strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
        }
        stxt <- frl(fpath)
    }
    # returns string w/o leading whitespace
    trim.leading <- function (x)  sub("^\\s+", "", x)

    # returns string w/o trailing whitespace
    trim.trailing <- function (x) sub("\\s+$", "", x)

    getKeyDepth = function(txtline) {
        key <- strsplit(txtline, ':')[[1]][1]
        print(key)
        counter <- 0
        while(grepl(paste0(rep('    ', counter + 1), collapse = ''), key)==TRUE) {
            counter <- counter + 1
        }
        return(counter)
    }

    getKV <- function(txtline) {
        v <- sub("^[^:]+:\\s*", "", txtline)
        v <- sub('\r', '', v)
        if(nchar(v)==0){
            k <- trim.leading(sub(':', '', txtline))
            k <- sub('\r', '', k)
            return(k)
        }
        else{
            k <- sub(v, "", txtline, fixed=T)
            k <- trim.leading(trim.trailing(sub(':', '', k)))
            k <- sub('\r', '', k)
            el <- list()
            el[[k]] <- v
            return(el)
        }
    }

    consecinds <- function (x) {
        if(length(x) == 0) {
            return(x)
        }
        else {
            cinds <- x[1]
            for(i in 2:length(x)) {
                if(x[i] == 1 + x[i-1]) {
                    cinds <- c(cinds, x[i])
                }
                else {
                    break
                }
            }
            return(cinds)
        }
    }

    buildnest <- function(inds, depth) {
        templist <- list()
        for(i in inds) {
            cval <- getKV(stxt[i])
            if(i == length(stxt)) {
                templist <- c(templist, cval)
            }
            else {
                if(keydepths[i] < keydepths[i + 1]) {
                    subinds <- which(keydepths == (depth + 1))
                    subinds <- subinds[subinds > i]
                    cdepths <- which(keydepths <= depth)
                    if(any(cdepths > i)) {
                        nextind <- min(cdepths[cdepths > i])
                        subinds <- subinds[subinds < nextind]
                    }
                    sublist<-buildnest(subinds, depth + 1)
                    templist[[cval]] <- sublist
                }
                else {
                    templist <- c(templist, cval)
                }
            }
        }
        return(templist)
    }

    keydepths <- as.numeric(sapply(stxt, getKeyDepth))
    print(keydepths)
    rootdepth <- 0
    rootinds <- which(keydepths==rootdepth)
    return(buildnest(rootinds, rootdepth))
}

readrcp <- function(fpath) {
    if(grepl('\\.zip', fpath)) {
        return(readnest(fpath, ext='.rcp', unzip=T))
    }
    else{
        return(readnest(fpath))
    }
}

readfom <- function(fompath, fn=F) {
    fomdt <- fread(fompath)
    if(fn){
        fomdt[,filename:=basename(fompath)]
    }
    if('Sample' %in% names(fomdt)) {
        fomdt[,sample_no:=Sample]
        fomdt[,Sample:=NULL]
    }
    return(fomdt)
}


getinfo <- function(plateno) {
    infopath <- file.path(jd, 'plate', plateno, paste0(plateno, '.info'))
    if(file.exists(infopath)) {
        return(readnest(infopath))
    }
    else {
        print(paste(infopath, 'not found.'))
        return(NA)
    }
}


# TESTING

# ladt<-fread('/home/dan/htehome/users/hte/catalysts on BVO/La/NiLaCoCe_29922_LampTrans_fixedcal_fixedspecinds_withTRANS64_AveAM15.csv')
# ladt<-ladt[,1:9, with=F]
#
#
# smoothladt<-smoothfom(ladt, testnnsym, 2)
# smoothladt
#
# system.time(testnnsym<-findnn('~/pms/0049-04-0830-mp.txt', targetcode=30, 1, spancodes=F))
# system.time(testnn<-findnn('~/pms/0049-04-0830-mp.txt', targetcode=30, 1, spancodes=F, sym=F))
# list.files('~/pms/')
#
# testnn[A==0 & B==0 & C==0]
# testnn[Sample %in% c(726, 730, 733,1000)]
# testnnsym[A==0 & B==0 & C==0]
# testnnsym[Sample %in% c(726, 730, 1000)]
#
# subpm[, .(A, B, C, D, list(which(if(A==min(subpm[,A])) subpm[,A]==min(subpm[,A]) else rep(T, nrow(subpm)) &
#                                      if(B==min(subpm[,B])) subpm[,B]==min(subpm[,B]) else rep(T, nrow(subpm)) &
#                                      if(C==min(subpm[,C])) subpm[,C]==min(subpm[,C]) else rep(T, nrow(subpm)) &
#                                      if(D==min(subpm[,D])) subpm[,D]==min(subpm[,D]) else rep(T, nrow(subpm))
#                                      ))), by = Sample]
