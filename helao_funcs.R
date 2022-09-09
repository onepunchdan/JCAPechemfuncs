library(rjson)
library(yaml)
library(bit64)


readhlo <- function(hlo_path) {
    lines <- readLines(hlo_path)
    starti = which(lines=="%%")
    dtbl <- rbindlist(lapply((starti+1):length(lines), function(i) as.data.table(fromJSON(lines[[i]]))), fill=T)
    my_int_handler <- function(x) {as.integer64(x)}
    meta <- read_yaml(text=paste0(lines[1:(starti-1)], collapse="\n"), handlers=list(int=my_int_handler))
    return(list(meta=meta, data=dtbl))
}