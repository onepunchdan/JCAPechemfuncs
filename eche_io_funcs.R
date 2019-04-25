library(data.table)
library(stringr)

readeche <- function(path) {

    dt <- fread(path, header=T, skip=14)
    f <- basename(path)
    # use R-friendly names
    names(dt) <- sub('%column_headings=', '', names(dt))
    names(dt) <- sub('\\)', '', names(dt))
    names(dt) <- sub('\\(', '.', names(dt))
    # assign filename-parsed metadata to data table
    suppressWarnings({
        dt[, `:=`(sample_no=as.integer(str_match(f, '^Sample([0-9]+)_')[,2]),
                  technique=str_match(f, '_([A-Z]+)[0-9]+\\.txt')[,2],
                  technum=as.integer(str_match(f, '_[A-Z]+([0-9]+)\\.txt')[,2]),
                  run=basename(sub(f, '', path))
        )]
    })
    return(dt)
}

scrubnames <- function(DT) {
    colnames <- names(DT)
    colnames <- sub('\\)', '', colnames)
    colnames <- sub('\\(', '.', colnames)
    colnames <- sub('%column_headings=', '', colnames)
    setnames(DT, colnames)
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
              #trans=as.integer(str_match(f, '_(TRANS|DARK)([0-9]+)_')[,3]),
              #time=strptime(str_match(f, '_([0-9]{8}\\.[0-9]{6})\\.[0-9]\\.opt')[,2], format='%Y%m%d.%H%M%S'),
              run=basename(sub(f, '', path))
              )]
    return(dt)
}

readechefiles <- function(l) {
    return(rbindlist(lapply(l, readeche)))
}

readechefiles2 <- function(l) {
  return(rbindlist(foreach(i=1:length(l), .export=c('readeche', 'fread', 'str_match')) %dopar% readeche(l[i])))
}

readechedir <- function(dir) {
  return(readechefiles(list.files(dir, pattern='.txt', full.names=T)))
}

readoptdir <- function(dir) {
  return(readoptfiles(list.files(dir, pattern='.opt', full.names=T)))
}

readxrf <- function(f) {
#     orbcsv<-fread(f, skip='Inte', drop=c('StagR', 'MatchName', 'MatchSig', 'CoatThk1', 'CoatThk2','CoatConc'))
    orbcsv<-as.data.table(read.csv(f, skip=6, check.names=F))
    orbcsv[,c('StagR', 'MatchName', 'MatchSig', 'CoatThk1', 'CoatThk2','CoatConc'):=NULL]
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
    
    if(length(indWt)==0){
        setnames(orbcsv, 'Inte', 'BatchLabel')
        return(list('Int'=orbcsv[, -extracol, with=F]))
    }
    else{
        keepInt<-c(indWt:(indLbl-1), extracol)
        keepWt<-c((indInt+1):(indWt), indAt:(indLbl-1), extracol)
        keepAt<-c((indInt+1):(indAt), extracol)
        setnames(orbcsv, 'Inte', 'BatchLabel')
        return(list('Int'=orbcsv[, -keepInt, with=F],
                    'Wt'=orbcsv[, -keepWt, with=F],
                    'At'=orbcsv[, -keepAt, with=F]))
    }
    
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
    outdf[,!(names(outdf) %in% c('Sample', 'sample_no', 'Label', 'StgLabel')) & apply(outdf, 2, is.numeric)]<-format(outdf[,!(names(outdf) %in% c('Sample', 'sample_no', 'Label', 'StgLabel')) & apply(outdf, 2, is.numeric)], scientific=T, digits=6)
    write.table(outdf, file=paste0(deparse(substitute(df)),'.csv'),
                sep=',', row.names=FALSE, na="NaN", quote=FALSE) #,eol="\r\n",
}


getexpfiles <- function(meta, type='pstat_files') {
    exp <- meta #readrcp(ep)
    runs <- grep('run__', names(exp), value=T)
    runs <- Filter(function(x) any(grepl('files_technique__', names(exp[[x]]))), runs)
    getfilelist <- function(run) {
        run_use <- exp[[run]]$run_use
        ftechs <- grep('files_technique__', names(exp[[run]]), value=T)
        makefiledt <- function(ft) {
            matlist <- lapply(exp[[run]][[ft]][[type]],
                              function(x) matrix(unlist(strsplit(x, ';')), ncol=5))
            newdt <- as.data.table(do.call(rbind, matlist))
            setnames(newdt, c('type', 'colnames', 'skip', 'nrows', 'Sample'))
            newdt[,filename:=names(matlist)]
            newdt[,skip:=as.integer(skip)]
            newdt[,nrows:=as.numeric(nrows)]
            newdt[,Sample:=as.numeric(Sample)]
            techstr <- unlist(strsplit(ft, '__'))[2]
            newdt[,tech:=sub('[0-9]+', '', techstr)]
            newdt[,technum:=as.numeric(sub('[A-Za-z]+', '', techstr))]
            newdt[,run_use:=run_use]
            return(newdt)
        }
        rbindlist(lapply(ftechs, makefiledt))
    }
    finaldt<-rbindlist(lapply(runs, function(x) getfilelist(x)[,`:=`(run_idx=x, run_path=file.path(jd, 'run', sub('^/', '', exp[[x]][['run_path']])))]))
    return(finaldt[order(technum, Sample)])
}

# readrawfromexp <- function(expdt) {
#     unzdir=file.path(unztemp, basename(dirname(expdt$run_path)), sub('\\.zip$','',basename(expdt$run_path)))
#     target=file.path(unzdir, expdt$filename)
#     if(!dir.exists(unzdir)) {
#         dir.create(unzdir, recursive=T)
#     }
#     if(!file.exists(target)) {
#         unzip(expdt$run_path, files=expdt$filename, exdir=unzdir, overwrite=T, junkpaths=T)
#     }
#     newdt<-fread(target, col.names=unlist(strsplit(expdt$colnames, ',')[[1]]))
#     newdt[, `:=`(Sample=as.numeric(trimws(expdt$Sample)),
#                  technum=as.numeric(trimws(expdt$technum)),
#                  tech=trimws(expdt$tech),
#                  filename=trimws(expdt$filename),
#                  type=trimws(expdt$type))]
#     names(newdt) <- sub('\\)', '', names(newdt))
#     names(newdt) <- sub('\\(', '.', names(newdt))
#     return(newdt[order(technum, Sample, t.s)])
# }

readrawfromexp <- function(expdt) {
    runs <- unique(expdt$run_idx)
    readzipped <- function(run) {
        rundt <- expdt[run_idx==run]
        run_path <- rundt[1]$run_path
        unztemp<-tempdir()
        unztemp<-gsub("\\\\", "/", unztemp)
        unzip(run_path, files=rundt$filename, exdir=unztemp, overwrite=T, junkpaths=T)
        newdt<-rbindlist(apply(rundt, 1, function(x) fread(file.path(unztemp, x['filename']),
                                                           # skip=x['skip'],
                                                           # nrows=x['nrows'],
                                                           col.names=unlist(strsplit(x['colnames'], ',')[[1]]))[, `:=`(Sample=as.numeric(trimws(x['Sample'])),
                                                                                                                       technum=as.numeric(trimws(x['technum'])),
                                                                                                                       tech=trimws(x['tech']),
                                                                                                                       filename=trimws(x['filename']),
                                                                                                                       type=trimws(x['type']))]))
        names(newdt) <- sub('\\)', '', names(newdt))
        names(newdt) <- sub('\\(', '.', names(newdt))
        do.call(file.remove, list(list.files(unztemp, full.names = TRUE)))
        return(newdt)
    }
    finaldt<-rbindlist(lapply(runs, function(x) readzipped(x)[,run_index:=x]))
    return(finaldt[order(technum, Sample)])
}                           
                              
readexpfiles <- function(expdt) {
    runs <- unique(expdt$run_idx)
    readzipped <- function(run) {
        rundt <- expdt[run_idx==run]
        run_path <- rundt[1]$run_path
        unztemp<-tempdir()
        unztemp<-gsub("\\\\", "/", unztemp)
        unzip(run_path, files=rundt$filename, exdir=unztemp, overwrite=T, junkpaths=T)
        runlist<-list()
        for(i in 1:dim(rundt)[1]){
            rdt<-fread(file.path(unztemp, rundt[i]$filename),
                   skip=rundt[i]$skip,
                   nrows=as.numeric(rundt[i]$nrows),
                   header=F,
                   colClasses=rep('numeric', length(unlist(strsplit(rundt[i]$colnames, ',')[[1]]))),
                   col.names=unlist(strsplit(rundt[i]$colnames, ',')[[1]]))[, `:=`(Sample=as.numeric(trimws(rundt[i]$Sample)),
                                                                                                                       technum=as.numeric(trimws(rundt[i]$technum)),
                                                                                                                       tech=trimws(rundt[i]$tech),
                                                                                                                       filename=trimws(rundt[i]$filename),
                                                                                                                       type=trimws(rundt[i]$type))]
            runlist<-c(runlist, list(rdt))
        }
        newdt<-rbindlist(runlist)
        names(newdt) <- sub('\\)', '', names(newdt))
        names(newdt) <- sub('\\(', '.', names(newdt))
        do.call(file.remove, list(list.files(unztemp, full.names = TRUE)))
        return(newdt)
    }
    finaldt<-rbindlist(lapply(runs, function(x) readzipped(x)[,run_index:=x]))
    return(finaldt[order(technum, Sample, t.s)])
}

getanafiles <- function(meta, type='inter_files') {
    ana <- meta #readrcp(ep)
    runs <- grep('ana__', names(ana), value=T)
    getfilelist <- function(run) {
        ana_name <- ana[[run]]$name
        ftechs <- grep('files_run__', names(ana[[run]]), value=T)
        makefiledt <- function(ft) {
            matlist <- lapply(ana[[run]][[ft]][[type]],
                              function(x) matrix(unlist(strsplit(x, ';')), ncol=5))
            newdt <- as.data.table(do.call(rbind, matlist))
            setnames(newdt, c('type', 'colnames', 'skip', 'nrows', 'Sample'))
            newdt[,filename:=names(matlist)]
            newdt[,skip:=as.integer(skip)]
            newdt[,nrows:=as.numeric(nrows)]
            newdt[,Sample:=as.numeric(Sample)]
            techstr <- ana[[run]]$technique
            newdt[,tech:=sub('[0-9]+', '', techstr)]
            newdt[,technum:=as.numeric(sub('[A-Za-z]+', '', techstr))]
            newdt[,ana_name:=ana_name]
            return(newdt)
        }
        rbindlist(lapply(ftechs, makefiledt))
    }
    finaldt<-rbindlist(lapply(runs, function(x) getfilelist(x)[,`:=`(ana_idx=x, ana_path=list.files(file.path(jd, 'analysis', 'eche'), pattern=ana$name, full.names=T)[1])]))
    return(finaldt[order(ana_idx, technum, Sample)])
}

                              
readanafiles <- function(anadt) {
    ana_path <- anadt[1]$ana_path
    unztemp<-tempdir()
    unztemp<-gsub("\\\\", "/", unztemp)
    unzip(ana_path, files=anadt$filename, exdir=unztemp, overwrite=T, junkpaths=T)
    newdt<-rbindlist(apply(anadt, 1, function(x) fread(file.path(unztemp, x['filename']),
                                                           skip=x['skip'],
#                                                            nrows=x['nrows'],
                                                           col.names=unlist(strsplit(x['colnames'], ',')[[1]]))[, `:=`(Sample=as.numeric(trimws(x['Sample'])),
                                                                                                                       technum=as.numeric(trimws(x['technum'])),
                                                                                                                       tech=trimws(x['tech']),
                                                                                                                       filename=trimws(x['filename'])
                                                                                                                      )]))
    names(newdt) <- sub('\\)', '', names(newdt))
    names(newdt) <- sub('\\(', '.', names(newdt))
    do.call(file.remove, list(list.files(unztemp, full.names = TRUE)))

    return(newdt[order(technum, Sample)])
}
                              
                               
readoptfiles <- function(expdt) {
    runs <- unique(expdt$run_idx)
    readzipped <- function(run) {
        rundt <- expdt[run_idx==run]
        run_path <- rundt[1]$run_path
        unztemp<-tempdir()
        unztemp<-gsub("\\\\", "/", unztemp)
        #command <- function(command){
        #  system(paste("cmd.exe /c", command))
        #}
        #command(paste('7z e -aoa', paste0('-o', unztemp), run_path))
        unzip(run_path, exdir=unztemp, overwrite=T, junkpaths=T)
        if('measure_time' %in% names(rundt)){
            newdt<-rbindlist(lapply(1:nrow(rundt), function(n) readopt(file.path(unztemp, rundt[n]$filename))[,`:=`(run_use=rundt[n]$run_use, filename=rundt[n]$filename, trans=rundt[n]$trans, time=rundt[n]$time, measure_time=rundt[n]$measure_time)]))
        } else {
            newdt<-rbindlist(apply(rundt, 1, function(x) readopt(file.path(unztemp, x['filename']))[,`:=`(run_use=x['run_use'], filename=x['filename'])]))
        }
        
        names(newdt) <- sub('\\)', '', names(newdt))
        names(newdt) <- sub('\\(', '.', names(newdt))
        do.call(file.remove, list(list.files(unztemp, full.names = TRUE)))
        return(newdt)
    }
    finaldt<-as.data.table(rbindlist(lapply(runs, function(x) readzipped(x)[,run_index:=x])))
    
    return(finaldt[order(run_use, technum, sample_no, trans, Wavelength_nm)])
}
                              
getpath <- function(timestamp, exper='eche', getexp=T) {
    processes='L:/processes'
    if(getexp) {
        subdir='experiment'
        ext='.exp'
    }
    else {
        subdir='analysis'
        ext='.ana'
    }
    path=list.files(file.path(processes, subdir, exper), pattern=timestamp, full.names=T)
    fpath=file.path(path, paste0(timestamp, ext))
    return(fpath)
}

findexpz <- function(timestamp, type='eche'){
    path=list.files(file.path('J:/hte_jcap_app_proto/experiment', type), pattern=paste0(timestamp, '.*\\.zip$'), full.names=T)
    fpath=path[which.max(as.numeric(gsub('^.*\\.copied-|P[DS]T\\.zip$', '', basename(path))))]
    return(fpath)
}

findanaz <- function(timestamp, type='eche'){
    path=list.files(file.path('J:/hte_jcap_app_proto/analysis', type), pattern=paste0(timestamp, '.*\\.zip$'), full.names=T)
    fpath=path[which.max(as.numeric(gsub('^.*\\.copied-|P[DS]T\\.zip$', '', basename(path))))]
    return(fpath)
}

findexp <- function(timestamp, type='eche'){
    path=list.files(file.path('L:/processes/experiment', type), pattern=paste0(timestamp), full.names=T)
    fpath=path[which.max(as.numeric(gsub('^.*\\.copied-|P[DS]T$|^.*\\.done$|^.*\\.run$', '', basename(path))))]
    return(list.files(fpath, pattern='\\.exp$', full.names=T))
}

findana <- function(timestamp, type='eche'){
    path=list.files(file.path('L:/processes/analysis', type), pattern=paste0(timestamp), full.names=T)
    fpath=path[which.max(as.numeric(gsub('^.*\\.copied-|P[DS]T$|^.*\\.done$|^.*\\.run$', '', basename(path))))]
    return(list.files(fpath, pattern='\\.ana$', full.names=T))
}

getrunfromexp <- function(exppath) {
    rcp<-readrcp(exppath)
    runs<-grep('run__', names(rcp), value=T)
    paths=list()
    for(run in runs){
        relrunpath<-dirname(rcp[[run]]$run_path)
        relrunpath<-substr(relrunpath, 2, nchar(relrunpath))
        paths<-c(paths, list.files(file.path('J:/hte_jcap_app_proto/run', relrunpath), pattern='\\.zip$', full.names=T))
    }
    return(paths)
}
                                            
readana<-function(ts, fomname='I.A_photo', dt=T){
    anapath<-trimws(getpath(ts, getexp = F))
    anafolder<-trimws(dirname(anapath))
    if(is.na(anapath) | is.na(anafolder) | anapath=='NA' | anapath==''){
        print(ts)
        return(NA)
    }
    meta<-readrcp(anapath)
    idxs<-grep('ana__', names(meta), value=T)
    csvlist<-list()
    for(i in idxs){
        blk<-meta[[i]]
        if('files_multi_run' %in% names(blk)) {
            fmr<-blk[['files_multi_run']]
            if('fom_files' %in% names(fmr)) {
                ff<-fmr[['fom_files']]
                fns<-names(ff)
                for(f in fns){
                    if(grep('csv_fom_file', ff[[f]])) {
                        mstr<-ff[[f]]
                        mparts<-strsplit(mstr, ';')[[1]]
                        if(grepl(fomname, mparts[2])){
                            csvlist<-c(csvlist, list(data.table(plate_id=meta$plate_ids, idx=i, dir=anafolder, fn=f, cols=mparts[2], skip=mparts[3], len=mparts[4])))
                        }
                    }
                }
            }
        }
        
    }
    if(dt){
#         csvdt<-rbindlist(csvlist)
#         return(rbindlist(apply(csvdt, 1, function(x) fread(file.path(x['dir'], x['fn']), skip=as.numeric(x['skip']), nrows=as.numeric(x['len']), col.names=strsplit(x['cols'], ',')[[1]])[,`:=`(plate_id=as.integer(x['plate_id']), ana=ts, idx=x['idx'])])))
        fomlist<-list()
        for(x in csvlist){
            if(file.exists(file.path(x$dir, x$fn))){
                fomdt <- fread(file.path(x$dir, x$fn), skip=as.numeric(x$skip), nrows=as.numeric(x$len), col.names=strsplit(x$cols, ',')[[1]])[,`:=`(plate_id=as.integer(x$plate_id), ana=ts, idx=x$idx)]
                setkey(fomdt, sample_no)
                fomlist<-c(fomlist, list(fomdt))
            }
        }
        if(length(fomlist)==0){
            return(NA)
        }
        else{
            return(rbindlist(fomlist, fill = T))
        }
    }
    else{
        return(csvdt)
    }
}

findsmps <- function(z, pm, fomname='I.A_photo', thresh=1E-8) {
    plate = z$plate_id
    anas = as.character(sapply(Filter(function(x) x$type=='eche', getinfo(plate)$analyses), function(x) sub('\\.copied.*$', '', basename(x$path))))
    if(length(anas)==0){
        return(NA)
        print(z)
    }
    else{
        analist <- list()
        for(a in anas) {
            smpdt<-readana(a, fomname=fomname)
            suppressWarnings({
                if(!is.na(smpdt)){
                    analist <- c(analist, list(smpdt))
                }
            })
        }
        anadt <- rbindlist(analist, fill = T)
    }
    setkey(anadt, sample_no)
    fomdt <- pm[anadt, nomatch=0]
    fomdt[,c('dx', 'dy', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'):=NULL]
    if(z$ymargin==0 & z$xmargin>0) {
        # vertical band
        fomdt[, in_range:=(x>(z$xmin-z$xmargin) & x<(z$xmax+z$xmargin))]
    }
    else if(z$xmargin==0 & z$ymargin>0) {
        # horizontal band
        fomdt[, in_range:=y>(z$ymin-z$ymargin) & y<(z$ymax+z$ymargin)]
    }
    else{
        m = tan(z$angle)
        uint = z$ymax - m*z$xmax
        lint = z$ymin - m*z$xmin
        # sloped band
        fomdt[, in_range:=y>(m*x+lint-z$ymargin) & y<(m*x+uint+z$ymargin)]
    }
    fomdt[,`:=`(has_photo=abs(I.A_photo)>=thresh,
                mpid=z$mpid,
                formula=z$formula,
                phasefrac=z$phasefrac,
                direction=z$direct,
                location=z$pstring,
                comments=z$comments
               )]

    return(fomdt)
}
                               
parsepos <- function(s) {
    s<-gsub('=','-',s)
    xy<-strsplit(s, ',')[[1]]
    xs<-strsplit(xy[1], '-')[[1]]
    ys<-strsplit(xy[2], '-')[[1]]
    xmin = as.numeric(xs[1])
    xmax = as.numeric(ifelse(length(xs)>1, xs[2], xs[1]))
    ymin = as.numeric(ys[1])
    ymax = as.numeric(ifelse(length(ys)>1, ys[2], ys[1]))
    return(c(xmin, xmax, ymin, ymax))
}


## NON-DT functions
##

# readana<-function(anapath){
#    frl = function(fname) { # faster readlines
#        s = file.info( fname )$size
#        buf = readChar( fname, s, useBytes=T)
#        strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
#    }
#    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
#    stxt <- frl(anapath)#

#    dflist<-list()
#    anapathsplit<-unlist(strsplit(anapath, '/'))
#    anafile<-anapathsplit[length(anapathsplit)]
#    anafolder<-sub(anafile, '', anapath)
#    csvlines <- unique(grep('.csv', stxt, fixed=T, value=T))
#    fomfiles <- sapply(csvlines, function(x) trim(strsplit(x, ':', fixed=T)[[1]][1]))
#    for(f in fomfiles){
#        ananame<-paste0('ana', '__', strsplit(f, '__')[[1]][2])
#        fompath<-paste0(anafolder, f)
#        fomcon<-file(fompath)
#        fomhead<-as.numeric(unlist(strsplit(readLines(fomcon, n=1), '\t')))
#        close(fomcon)
#        fomtbl<-read.table(fompath, header=T, sep=',', na.strings='NaN', skip=fomhead[4]+1, nrows=fomhead[3]+1)
#        dflist[[ananame]]<-I(as.data.frame(fomtbl))
#    }
#     }
#    return(dflist)
# }

# Takes a single .rcp file path and returns a list with nested lists for params,
# and filenames.
                                      
readnest <- function(fpath, ext, unzip=F) {
    if(unzip) {
        fn<-strsplit(fpath, '.', fixed=T)[[1]]
        rcpfile<-basename(paste0(tail(strsplit(fn[1], '\\', fixed=T)[[1]], n=1), '.', fn[2], ext))
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
    
    stxt<-stxt[sapply(stxt, nchar) > 0]
    
    # returns string w/o leading whitespace
    trim.leading <- function (x)  sub("^\\s+", "", x)

    # returns string w/o trailing whitespace
    trim.trailing <- function (x) sub("\\s+$", "", x)

    getKeyDepth = function(txtline) {
        key <- strsplit(txtline, ':')[[1]][1]
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
    rootdepth <- 0
    rootinds <- which(keydepths==rootdepth)
    return(buildnest(rootinds, rootdepth))
}

readrcp <- function(fpath) {
    if(grepl('\\.zip', fpath)) {
        if(basename(dirname(dirname(fpath)))=='experiment'){
            fext='.exp'
        } else if(basename(dirname(dirname(fpath)))=='analysis'){
            fext='.ana'
        } else {
            fext='.rcp'
        }
        rcplist<-readnest(fpath, ext=fext, unzip=T)
    }
    else{
        rcplist<-readnest(fpath)
    }
    return(rcplist)
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

getruns <- function(plateno) {
    infolist <- getinfo(plateno)

}

# library(doSNOW)
# cl<-makeCluster(4)
# registerDoSNOW(cl)
# system.time(for(i in names(getinfo(2154)$runs)) getinfo(2154)$runs[[i]]$path)
# system.time(lapply(names(getinfo(2154)$runs), function(i) getinfo(2154)$runs[[i]]$path))
# system.time(foreach(i=names(getinfo(2154)$runs)) %do% getinfo(2154)$runs[[i]]$path)
# system.time(foreach(i=names(getinfo(2154)$runs)) %dopar% getinfo(2154)$runs[[i]]$path)

### RCP list functions to write
# (1) get techniques
gettechs <- function(rcp) {
    techs <- sub('echem_params__', '', grep('echem_params__(.*[0-9])', names(rcp), value=T))
    return(techs)
}

# (2) get technique params
getparams <- function(rcp, tech) {
    param.ind <- grep(paste0('echem_params__', tech), names(rcp))
    params <- rcp[[param.ind]]
    return(params)
}

# (3) get list of files/samples
getpstatfiles <- function(rcp, tech) {
    file.ind <- grep(paste0('files_technique__', tech), names(rcp))
    flist <- rcp[[file.ind]]$pstat_files
    return(flist)
}

getspecfiles <- function(rcp, tech) {
    file.ind <- grep(paste0('files_technique__', tech), names(rcp))
    flist <- rcp[[file.ind]]$spectrum_files
    return(flist)
}


#library(foreach)
#library(doSNOW)

readpstatfiles.rcp <- function(rcp, tech, samples=NULL) {
    flist <- getpstatfiles(rcp, tech)
    slist <- as.integer(str_match(flist, ';([0-9]*)$')[,2])
    importsamples <- if(is.null(samples)) slist else samples
    rcpfolder <- sub(basename(rcp$rcppath), '', rcp$rcppath)
    DT.files <- data.table(Sample=slist, fn=names(flist), val=as.character(flist))
    DT.new <- rbindlist(foreach(i=importsamples, .packages=c('data.table', 'stringr')) %dopar% {
        fread(file.path(rcpfolder, DT.files[Sample==i, fn]),
              col.names=sub('\\)', '',
                            sub('\\(', '.',
                                strsplit(strsplit(DT.files[Sample==i, val], ';')[[1]][2], ',')[[1]])))[,`:=`(Sample=i, Eta.V=Ewe.V-as.numeric(rcp$reference_e0))]
    })
    return(DT.new)
}


splitcol2rows_mget <- function(dtInput, col2split, sep){
  dtInput <- dtInput[, .(tmp.add.col = unlist(strsplit(get(col2split),sep,T))), by=names(dtInput)]

  dtInput[, c(col2split):=NULL];
  setnames(dtInput, 'tmp.add.col', col2split); 
  return(dtInput);
}
                                            
getcompdt<-function(infodt){
    prec<-rbindlist(lapply(infodt$prints, function(r) as.data.table(data.frame(r))), fill=T)[id==infodt$screening_print_id]
    channellist=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
    els<-strsplit(as.character(prec$elements), split=',', fixed=T)[[1]]
    pm<-getpm(prec$map_id)
    setkey(pm, Sample)
                           
    concdt<-data.table(pm[,.(Sample, A, B, C, D, E, F, G, H)])
    setkey(concdt, Sample)
    concdt<-concdt[,1:(length(els)+1),with=F]
    names(concdt)[2:(length(els)+1)]<-els
                           
    if('concentration_values' %in% names(prec)){
    cels<-strsplit(as.character(prec$concentration_elements), split=',', fixed=T)[[1]]
    concs<-as.numeric(strsplit(as.character(prec$concentration_values), split=',', fixed=T)[[1]])
    loaddt<-data.table(cels, concs)

        if(length(cels)>=length(els)){
            compdt<-data.table(compstr=trimws(strsplit(gsub("([A-Z])", " \\1", els), ",")), channels=channellist)
            scompdt<-splitcol2rows_mget(compdt, 'compstr', ' ')
            lscompdt<-cbind(scompdt, ldt)
            setkey(lscompdt, channels)

            calccomps<-function(smp){
                channelconcs<-melt(pm[Sample==smp], measure.vars=channellist)[,.(variable, value)]
                setkey(channelconcs, variable)
                total<-channelconcs[lscompdt][,.(cels, conc=value*concs)][,.(total=sum(conc)),by=cels]
                frow<-data.table(Sample=smp, transpose(total)[2])
                names(frow)[2:dim(frow)[2]]<-total[[1]]
                return(frow)
            }

            if(all(lscompdt$compstr==lscompdt$cels)){
                concdt<-rbindlist(lapply(pm$Sample, calccomps))
            }
        }
    }
    return(concdt)
}

getalloys<-function(infometa){
    printblk<-infometa$prints[[names(infometa$prints)[sapply(names(infometa$prints), function(n) infometa$prints[[n]]$id==infometa$screening_print_id)]]]
    return(strsplit(printblk$elements, ',')[[1]][3:8])
}

get_rawlendt<-function(meta, anablocks){
    rawlen_fns<-character()
    for(blk in anablocks){
        frunblocks<-grep('files_run__', names(meta[[blk]]), value=T)
        for(fblk in frunblocks){
            rawlen_fns<-c(rawlen_fns, names(meta[[blk]][[fblk]]$inter_rawlen_files))
        }
    }
    rawldt<-data.table(f=rawlen_fns, Sample=(sapply(rawlen_fns, function(x) str_match(x, 'Sample([0-9]+)')[2])))
    rawldt<-rawldt[order(Sample)]
    irawdt<-rbindlist(lapply(1:dim(rawldt)[1], function(i) fread(file.path(dirname(getpath(anameta$name, getexp=F)), rawldt[i]$f))[,Sample:=rawldt[i]$Sample]))
    #exp<-getexpfiles(readrcp(getpath(meta$experiment_name)))
    #rawdt<-readexpfiles(exp[order(Sample)])
    return(irawdt)
}
                             
                             
#platemap 72 functions
genrcoords<-function(pm, v1, v2, v3, v4, v5, v6, v7, v8) {
    AB<-copy(pm[get(eval(v3))+get(eval(v4))+get(eval(v5))+get(eval(v6))==0][get(eval(v1))>0 & get(eval(v2))>0])
    AB[,rad:=(pi/6)+(pi/3)*(get(eval(v1))/(get(eval(v1))+get(eval(v2))))]
    
    BC<-copy(pm[get(eval(v1))+get(eval(v4))+get(eval(v5))+get(eval(v6))==0][get(eval(v2))>0 & get(eval(v3))>0])
    BC[,rad:=(11*pi/6)+(pi/3)*(get(eval(v2))/(get(eval(v2))+get(eval(v3))))]
    
    CD<-copy(pm[get(eval(v1))+get(eval(v2))+get(eval(v5))+get(eval(v6))==0][get(eval(v3))>0 & get(eval(v4))>0])
    CD[,rad:=(3*pi/2)+(pi/3)*(get(eval(v3))/(get(eval(v3))+get(eval(v4))))]
    
    DE<-copy(pm[get(eval(v1))+get(eval(v2))+get(eval(v3))+get(eval(v6))==0][get(eval(v4))>0 & get(eval(v5))>0])
    DE[,rad:=(7*pi/6)+(pi/3)*(get(eval(v4))/(get(eval(v4))+get(eval(v5))))]
    
    EF<-copy(pm[get(eval(v1))+get(eval(v2))+get(eval(v3))+get(eval(v4))==0][get(eval(v5))>0 & get(eval(v6))>0])
    EF[,rad:=(5*pi/6)+(pi/3)*(get(eval(v5))/(get(eval(v5))+get(eval(v6))))]
    
    FA<-copy(pm[get(eval(v2))+get(eval(v3))+get(eval(v4))+get(eval(v5))==0][get(eval(v6))>0 & get(eval(v1))>0])
    FA[,rad:=(pi/2)+(pi/3)*(get(eval(v6))/(get(eval(v6))+get(eval(v1))))]
    
    A<-copy(pm[get(eval(v2))+get(eval(v3))+get(eval(v4))+get(eval(v5))+get(eval(v6))==0][get(eval(v1))>0])
    A[,rad:=(pi/2)]
    B<-copy(pm[get(eval(v1))+get(eval(v3))+get(eval(v4))+get(eval(v5))+get(eval(v6))==0][get(eval(v2))>0])
    B[,rad:=(pi/6)]
    C<-copy(pm[get(eval(v1))+get(eval(v2))+get(eval(v4))+get(eval(v5))+get(eval(v6))==0][get(eval(v3))>0])
    C[,rad:=(11*pi/6)]
    D<-copy(pm[get(eval(v1))+get(eval(v2))+get(eval(v3))+get(eval(v5))+get(eval(v6))==0][get(eval(v4))>0])
    D[,rad:=(3*pi/2)]
    E<-copy(pm[get(eval(v1))+get(eval(v2))+get(eval(v3))+get(eval(v4))+get(eval(v6))==0][get(eval(v5))>0])
    E[,rad:=(7*pi/6)]
    F<-copy(pm[get(eval(v1))+get(eval(v2))+get(eval(v3))+get(eval(v4))+get(eval(v5))==0][get(eval(v6))>0])
    F[,rad:=(5*pi/6)]
    
    outdt<-unique(rbindlist(list(AB, BC, CD, DE, EF, FA, A, B, C, D, E, F)))
    outdt[,`:=`(xr=cos(rad), yr=sin(rad))]
    base<-copy(pm[get(eval(v1))+get(eval(v2))+get(eval(v3))+get(eval(v4))+get(eval(v5))+get(eval(v6))==0][get(eval(v7)) + get(eval(v8))>0])
    base[,`:=`(rad=NA, xr=0, yr=0)]
    outdt2<-rbindlist(list(outdt, base))
    return(outdt2)
}
                             
makelab<-function(pmdt, elements) {
    channels<-c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
    for(i in 1:dim(pmdt)[1]){
        chvals=transpose(pmdt[i][,channels[1:length(elements)],with=F])$V1
        chvals=chvals/sum(chvals)
        chinds=chvals>0
        pmdt[i, lab:=paste0(paste0(elements[chinds], '[', chvals[chinds], ']'), collapse='*')]
    }
    return(pmdt)
}
                             
findexpz <- function(timestamp, type='eche'){
    path=list.files(file.path('J:/hte_jcap_app_proto/experiment', type), pattern=paste0(timestamp, '.*\\.zip$'), full.names=T)
    fpath=path[which.max(as.numeric(gsub('^.*\\.copied-|P[DS]T\\.zip$', '', basename(path))))]
    return(fpath)
}

findanaz <- function(timestamp, type='eche'){
    path=list.files(file.path('J:/hte_jcap_app_proto/analysis', type), pattern=paste0(timestamp, '.*\\.zip$'), full.names=T)
    fpath=path[which.max(as.numeric(gsub('^.*\\.copied-|P[DS]T\\.zip$', '', basename(path))))]
    return(fpath)
}

findexp <- function(timestamp, type='eche'){
    path=list.files(file.path('L:/processes/experiment', type), pattern=paste0(timestamp), full.names=T)
    fpath=path[which.max(as.numeric(gsub('^.*\\.copied-|P[DS]T$|^.*\\.done$|^.*\\.run$', '', basename(path))))]
    return(list.files(fpath, pattern='\\.exp$', full.names=T))
}

findana <- function(timestamp, type='eche'){
    path=list.files(file.path('L:/processes/analysis', type), pattern=paste0(timestamp), full.names=T)
    fpath=path[which.max(as.numeric(gsub('^.*\\.copied-|P[DS]T$|^.*\\.done$|^.*\\.run$', '', basename(path))))]
    return(list.files(fpath, pattern='\\.ana$', full.names=T))
}

readrcp <- function(fpath) {
    if(grepl('\\.zip', fpath)) {
        if(basename(dirname(dirname(fpath)))=='experiment'){
            fext='.exp'
        } else if(basename(dirname(dirname(fpath)))=='analysis'){
            fext='.ana'
        } else {
            fext='.rcp'
        }
        rcplist<-readnest(fpath, ext=fext, unzip=T)
    }
    else{
        rcplist<-readnest(fpath)
    }
    return(rcplist)
}