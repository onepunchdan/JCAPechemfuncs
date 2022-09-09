get_V_fracs <- function(V_frac) {
    h <- hist(V_frac, breaks = quantile(V_frac, seq(0, 1, length.out = 4)), plot = F)
    f <- factor(cut(V_frac, breaks = h$breaks, include.lowest = T))
    mediandt <- data.table(actual = V_frac, grouped = f)
    mediandt <- mediandt[, .(V_frac_median = median(actual)), by = grouped]
    medianfracs <- round(as.numeric(as.character(factor(f, levels = mediandt$grouped, labels = mediandt$V_frac_median))), 2)
    return(medianfracs)
}

get_unary_loadings <- function(alloysum) {
    h <- hist(alloysum, breaks = quantile(alloysum, seq(0, 1, 0.1)), plot = F)
    f <- factor(cut(alloysum, breaks = h$breaks, include.lowest = T))
    mediandt <- data.table(actual = alloysum, grouped = f)
    mediandt <- mediandt[, .(alloysum_median = median(actual)), by = grouped]
    medianloadings <- round(as.numeric(as.character(factor(f, levels = mediandt$grouped, labels = mediandt$alloysum_median))), 3)
    return(medianloadings)
}

bingroupsizes <- c(90, 180, 450, 810)
get_binary_loadings <- function(alloysum) {
    h <- hist(alloysum, breaks = quantile(alloysum, c(0, cumsum(bingroupsizes)/sum(bingroupsizes))), plot = F)
    f <- factor(cut(alloysum, breaks = h$breaks, include.lowest = T))
    mediandt <- data.table(actual = alloysum, grouped = f)
    mediandt <- mediandt[, .(alloysum_median = median(actual)), by = grouped]
    medianloadings <- round(as.numeric(as.character(factor(f, levels = mediandt$grouped, labels = mediandt$alloysum_median))), 3)
    return(medianloadings)
}

# anats_ <- '20200226.133715'
# anaint_ = 1
# fom_name_ = "I.A_photothresh"
# aggfunc_ = mean
# min_ufom_ = 1e-7
# max_ufom_ = 6e-6
# min_bfom_ = 1e-7
# max_bfom_ = 6e-6
# ylabel_ = expression(I[1.23*V] ~ (A))
# omit_codes_ = c()
# title_add_ = ' in OER9 at 1.23V'

make_panels <- function (anats, anaint, fom_name, aggfunc, min_ufom, max_ufom, min_bfom, max_bfom, ylabel, omit_codes, title_add) {
    analst <- readanalst(anats)
    fdt <- analst[[paste0("ana__", anaint)]]
    setkey(fdt, sample_no)
    fdt_sub <- fdt[, .(sample_no, fom = get(fom_name))]

    aggstring <- suppressWarnings(sub('\\..*$', '', methods(aggfunc)[1]))
    anameta <- readrcp(findana(anats))
    plate_id <- ifelse(is.null(anameta$plate_ids), sub("^.*plate_id ", "", anameta[[paste0("ana__", anaint)]]$description), anameta$plate_ids)
    infometa <- getinfo(plate_id)
    pm <- getpm(infometa$screening_map_id)
    compdt <- getcompdt(infometa)
    compdt[, `:=`(V_frac, round(V/(Cu + V), 3))]
    compdt <- compdt[V_frac != "NaN"]
    compdt[, `:=`(V_frac_grp, as.factor(get_V_fracs(V_frac)))]
    compdt <- compdt[Sample %in% pm[!code %in% omit_codes]$Sample]

    alloys <- getalloys(infometa)
    binlist <- split(permutations(6, 2, alloys), row(permutations(6, 2, alloys)))
    ulist <- list()
    for (x in alloys) {
        omit <- alloys[alloys != x]
        sdt <- copy(compdt)
        sdt <- sdt[get(x) > 0 & get(omit[1]) + get(omit[2]) + get(omit[3]) + get(omit[4]) + get(omit[5]) == 0]
        sdt[, `:=`(A_load, get(x))]
        sdt[, `:=`(A_lab, x)]
        ulist <- c(ulist, list(sdt))
    }
    blist <- list()
    for (b in binlist) {
        omit <- alloys[!alloys %in% b]
        sdt <- copy(compdt)
        sdt <- sdt[get(omit[1]) + get(omit[2]) + get(omit[3]) + get(omit[4]) == 0]
        sdt <- sdt[get(b[1]) > 0 & get(b[2]) > 0]
        sdt[, `:=`(A_frac, get(b[1])/(get(b[1]) + get(b[2])))]
        sdt[, `:=`(A_lab, b[1])]
        sdt[, `:=`(B_alloy, b[2])]
        blist <- c(blist, list(sdt))
    }
    pdt <- copy(compdt)
    pdt <- na.omit(pdt[get(alloys[1]) + get(alloys[2]) + get(alloys[3]) + get(alloys[4]) + get(alloys[5]) + get(alloys[6]) == 0])
    udt <- rbindlist(ulist)
    bdt <- rbindlist(blist)
    udt[, `:=`(alloy_load, get_unary_loadings(alloy_sum))]
    bdt[, `:=`(alloy_load, get_binary_loadings(alloy_sum))]
    pdt_sub <- pdt[, .(Sample, V_frac_grp)]
    udt_sub <- udt[, .(Sample, V_frac_grp, A_lab, alloy_frac, alloy_load)]
    bdt_sub <- bdt[, .(Sample, V_frac_grp, A_lab, B_alloy, A_frac, alloy_frac, alloy_load)]
    setkey(pdt_sub, Sample)
    setkey(udt_sub, Sample)
    setkey(bdt_sub, Sample)
    analst <- readanalst(anats)
    fdt <- analst[[paste0("ana__", anaint)]]
    setkey(fdt, sample_no)
    fdt_sub <- fdt[, .(sample_no, fom = get(fom_name))]
    pfoms <- fdt_sub[pdt_sub][, .(fom = aggfunc(fom)), by = V_frac_grp]
    ufoms <- fdt_sub[udt_sub][, .(fom = aggfunc(fom)), by = list(V_frac_grp, A_lab, alloy_frac, alloy_load)]
    bfoms <- fdt_sub[bdt_sub][, .(fom, V_frac_grp, A_lab, B_alloy, A_frac, alloy_frac, alloy_load)]
    bloads <- unique(bfoms$alloy_load)
    bcombs <- unique(bfoms[, .(V_frac_grp, A_lab, B_alloy, alloy_load)])
    setkey(bcombs, V_frac_grp, A_lab, alloy_load)
    uends <- ufoms[alloy_load %in% bloads]
    uends[, `:=`(A_frac, 1)]
    setkey(uends, V_frac_grp, A_lab, alloy_load)
    endfoms <- uends[bcombs][, .(fom, V_frac_grp, A_lab, B_alloy, A_frac, alloy_frac, alloy_load)]
    ustarts <- ufoms[alloy_load %in% bloads]
    names(ustarts)[2] <- "B_alloy"
    ustarts[, `:=`(A_frac, 0)]
    setkey(ustarts, V_frac_grp, B_alloy, alloy_load)
    setkey(bcombs, V_frac_grp, B_alloy, alloy_load)
    startfoms <- ustarts[bcombs][, .(fom, V_frac_grp, A_lab, B_alloy, A_frac, alloy_frac, alloy_load)]
    ubfoms <- rbindlist(list(bfoms, endfoms, startfoms))[order(V_frac_grp, A_lab, B_alloy, alloy_load, A_frac)]
    ugrps <- unique(ufoms[, .(V_frac_grp, A_lab)])
    setkey(ugrps, V_frac_grp)
    setkey(pfoms, V_frac_grp)
    purefoms <- pfoms[ugrps][, .(V_frac_grp, A_lab, alloy_frac = 0, alloy_load = 0, fom)]
    pufoms <- rbindlist(list(ufoms, purefoms))[order(V_frac_grp, A_lab, alloy_load)]
    pufoms[, `:=`(blank_factor, "")]

    uplot <- ggplot(pufoms, aes(x = alloy_load, y = fom, color = A_lab)) + geom_point() + geom_path() + facet_grid(blank_factor ~ paste("X =", 
        V_frac_grp)) + theme_bw(base_size = 16) + theme(strip.background = element_blank()) + scale_x_continuous(breaks = seq(0, 0.02, 0.004)) + xlab(expression(y ~ 
        "in" ~ Cu[1 - ~x] * V[x] * O[delta]:A[y])) + ylab(ylabel) + scale_color_manual(name = expression(A), values = gg.colors(6)) + ggtitle(paste(paste(c("Cu", 
        "V", alloys), collapse = "-"), paste0(plate_id, title_add, ": "), paste(aggstring, fom_name))) + ylim(ifelse(is.na(min_ufom), min(na.omit(pufoms)$fom), min_ufom), ifelse(is.na(max_ufom), max(na.omit(pufoms)$fom), max_ufom))
    bplot <- ggplot(ubfoms, aes(x = 1 - A_frac, y = fom, color = B_alloy, shape = factor(alloy_load))) + geom_point() + geom_path() + facet_grid(A_lab ~ 
        paste("X =", V_frac_grp)) + theme_bw(base_size = 16) + scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + xlab(expression(z ~ "in" ~ 
        Cu[1 - ~x] * V[x] * O[delta]:(A[1 - ~z] * B[z])[y])) + ylab(ylabel) + theme(strip.background = element_blank(), strip.text.x = element_blank()) + 
        scale_shape_manual(name = expression(y), values = c(6, 4, 1, 2)) + scale_color_manual(name = expression(B), values = gg.colors(6)) + ylim(ifelse(is.na(min_bfom), min(na.omit(ubfoms)$fom), min_bfom), ifelse(is.na(max_bfom), max(na.omit(ubfoms)$fom), max_bfom))
    ugrob <- ggplotGrob(uplot)
    bgrob <- ggplotGrob(bplot)
    g <- rbind(ugrob, bgrob, size = "first")
    g$widths <- unit.pmax(ugrob$widths, bgrob$widths)

    library(gridExtra)
    options(repr.plot.width=12, repr.plot.height=16)
    ugrob <- ggplotGrob(uplot)
    bgrob <- ggplotGrob(bplot)
    g <- rbind(ugrob, bgrob, size = "first")
    g$widths <- unit.pmax(ugrob$widths, bgrob$widths)
    aggstr <- gsub('"', "", sub('\\).*$', "", sub('^.*\\(', "", deparse(aggfunc)[2])))
    figfn <- paste(plate_id, anats, 'ana_', anaint, aggstr, fom_name, sep='_')
    figpng <- paste0(figfn, '.pdf')
    ggsave(figpng, g, device="pdf", dpi=200, width=12, height=16, units="in")
}

