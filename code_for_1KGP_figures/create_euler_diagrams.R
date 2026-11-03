library(eulerr)
library(UpSetR)
library(jsonlite)

sum_snps <- function(pop_combo, pops, counts, cutoff) {
  pop_tf <- rep(FALSE, length(pops))
  pop_tf[pop_combo] <- TRUE
  populations <- pops[pop_combo]
  if (length(pop_combo) > 1) {
    shared_common_snps <- which(rowSums(counts[,populations]>=cutoff)==length(populations))
  } else {
    shared_common_snps <- which(counts[,populations]>=cutoff)
  }
  snp_counts <- data.frame(t(c(pop_tf, sum(counts[shared_common_snps, ncol(counts)]))))
  colnames(snp_counts) <- c(pops, "shared_common_snps")
  return(snp_counts)
}

generate_euler_plot <- function(pairwise_file, poplist, selected_pops=c(), cutoff=3, common_pop="") {
  data <- read.table(
    file=pairwise_file,
    col.names=c("geovar_code", "counts"),
    colClasses=c("character", "numeric")
  )
  pops <- read.table(file=poplist)$V1
  if (length(selected_pops) < 1) {
    selected_pops <- pops
  }
  
  split.pops <- sapply(data[,1], function(a){strsplit(as.character(a),"")[[1]]})
  split.pops <- apply(split.pops,2,as.numeric)
  split.pops <- as.data.frame(t(split.pops))
  row.names(split.pops) <- 1:nrow(split.pops)
  colnames(split.pops) <- c(as.character(pops))
  data <- cbind(split.pops, data)
  data <- data[,c(selected_pops,"geovar_code","counts")]
  
  if (common_pop!=""){
    data <- data[which(data[common_pop]>=cutoff),]
  }
  
  combo_func <- Map(combn, list(1:length(selected_pops)), seq_along(1:length(selected_pops)), simplify=FALSE)
  combinations <- unlist(combo_func, recursive=FALSE)
  
  shared_common_snps <- do.call(rbind,lapply(combinations, FUN=sum_snps, pops=selected_pops, counts=data, cutoff=cutoff))
  shared_common_snps[nrow(shared_common_snps),"unique_snps"] <- shared_common_snps[nrow(shared_common_snps),"shared_common_snps"]
  
  for (i in (nrow(shared_common_snps)-1):1) {
    shared_common_snps[i,"unique_snps"] <- shared_common_snps[i,"shared_common_snps"] - sum(
      shared_common_snps[
        apply(
          shared_common_snps[,1:(ncol(shared_common_snps)-2)] - shared_common_snps[rep(i,nrow(shared_common_snps)),1:(ncol(shared_common_snps)-2)], 
          MARGIN=1, 
          FUN=function(x){!any(x<0)}
        ),"unique_snps"
      ], na.rm=TRUE
    )
  }
  
  euler_data <- shared_common_snps$unique_snps
  names(euler_data) <- apply(shared_common_snps[selected_pops], 1, function(x){paste(names(x)[which(x==1)],collapse="&")})
  return(list("sets"=euler_data,"euler"=eulerr::euler(euler_data, shape = "ellipse"),"size"=sum(shared_common_snps$unique_snps)))
}

generate_d3_euler <- function(pops, common_pop="", supplemental="") {
  if (common_pop != "") {
    pops_euler <- generate_euler_plot("../data/biddanda_et_al_2019/new_1kg_nyc_hg38_filt_total.biallelic_snps.pops.ncat_000_001_005_010_100.geodist.total.txt", "../data/biddanda_et_al_2019/poplists/pops_panel.txt", selected_pops=pops, common_pop=common_pop)
  } else {
    pops_euler <- generate_euler_plot("../data/biddanda_et_al_2019/new_1kg_nyc_hg38_filt_total.biallelic_snps.pops.ncat_000_001_005_010_100.geodist.total.txt", "../data/biddanda_et_al_2019/poplists/pops_panel.txt", selected_pops=pops) 
  }
  if (supplemental == "") {
    supplemental <- read.table(file="../data/biddanda_et_al_2019/supplemental.txt", sep="\t", header=TRUE)
    supplemental <- supplemental[,c("Population.description", "Population.code..1KGP.", "NumberOfSampledIndividualsPhase3")]
    colnames(supplemental) <- c("description", "abbreviation", "sampled_individuals")
  }
  ellipses_coordinates <- pops_euler$euler$ellipses
  ellipses_coordinates$abbreviation <- row.names(ellipses_coordinates)
  ellipses_coordinates$common_variants <- sapply(pops, FUN=function(pop, euler) {return(sum(euler[grepl(pop, names(euler))]))}, euler=pops_euler$sets)
  ellipses_coordinates$unshared_common_variants <- pops_euler$sets[pops]
  ellipses_coordinates <- merge(ellipses_coordinates, supplemental, by="abbreviation")
  # Color Palette from https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7
  colors <- c("#D55E00", "#0072B2", "#CC79A7", "#009E73", "#5E159A", "#56B4E9", "#E69F00", "#000000")
  ellipses_coordinates$color <- colors[1:length(pops)]
  ellipses_coordinates$fill <- "none"
  ellipses_coordinates$stroke_dasharray <- "none"
  return(list("ellipses_coord"=ellipses_coordinates,"stress"=pops_euler$euler$stress, "diagError"=pops_euler$euler$diagError))
}

AMR <- c("ACB", "ASW", "CEU", "CLM", "GIH", "MXL", "PEL", "PUR")
colors <- c("#D55E00", "#0072B2", "#CC79A7", "#009E73", "#5E159A", "#56B4E9", "#E69F00", "#000000")

supplemental <- read.table(file="../data/biddanda_et_al_2019/supplemental.txt", sep="\t", header=TRUE)
supplemental <- supplemental[,c("Population.description", "Population.code..1KGP.", "NumberOfSampledIndividualsPhase3")]
colnames(supplemental) <- c("description", "abbreviation", "sampled_individuals")







# Figure 5A

measurable <- 2.9e9

data <- read.table(
  file="../data/biddanda_et_al_2019/new_1kg_nyc_hg38_filt_total.biallelic_snps.pops.ncat_000_001_005_010_100.geodist.total.txt",
  col.names=c("geovar_code", "counts"),
  colClasses=c("character", "numeric")
)
pops <- read.table(file="../data/biddanda_et_al_2019/poplists/pops_panel.txt")$V1
if (length(AMR) < 1) {
  AMR <- pops
}
split.pops <- sapply(data[,1], function(a){strsplit(as.character(a),"")[[1]]})
split.pops <- apply(split.pops,2,as.numeric)
split.pops <- as.data.frame(t(split.pops))
row.names(split.pops) <- 1:nrow(split.pops)
colnames(split.pops) <- c(as.character(pops))
data <- cbind(split.pops, data)
data <- data[,c(AMR,"geovar_code","counts")]
amr_num_variants <- sum(data$counts) - sum(data[which(rowSums(data[,1:length(AMR)]) == 0),"counts"])

amr_pops_euler <- generate_euler_plot("../data/biddanda_et_al_2019/new_1kg_nyc_hg38_filt_total.biallelic_snps.pops.ncat_000_001_005_010_100.geodist.total.txt", "../data/biddanda_et_al_2019/poplists/pops_panel.txt", selected_pops=AMR)
perspective <- amr_pops_euler$sets
perspective_common <- sum(perspective)
names(perspective) <- paste("Measurable&Variants&", names(perspective), sep="")
perspective <- c(perspective, "Measurable"=(measurable - amr_num_variants))
perspective <- c(perspective, "Measurable&Variants"=(amr_num_variants - perspective_common))
perspective_euler <- eulerr::euler(
  perspective,
  shape = "ellipse"
)
plot(perspective_euler, lwd=1, fill="transparent", labels=c())
persepective_ellipses_coordinates <- perspective_euler$ellipses
persepective_ellipses_coordinates$color <- c("#A9A9A9", "#A9A9A9", colors)
persepective_ellipses_coordinates$fill <- "none"
persepective_ellipses_coordinates$stroke_dasharray <- c(5,5,rep("none", length(AMR)))
toJSON(persepective_ellipses_coordinates)



# Figure 2C
simplified_perspective <- c("Measurable"=(measurable - amr_num_variants), "Measurable&Variants"=(amr_num_variants - perspective_common), "Measurable&Variants&Common"=perspective_common)
simplified_perspective_euler <- eulerr::euler(
  simplified_perspective,
  shape = "ellipse"
)
plot(simplified_perspective_euler, lwd=1, fill="transparent", labels=c())
simplified_persepective_ellipses_coordinates <- simplified_perspective_euler$ellipses
simplified_persepective_ellipses_coordinates$color <- c("#A9A9A9", "#A9A9A9", "none")
simplified_persepective_ellipses_coordinates$fill <- c("none", "none", "#56B4E9")
simplified_persepective_ellipses_coordinates$stroke_dasharray <- c(5,5,"none")
toJSON(simplified_persepective_ellipses_coordinates)



# Figure 4A
separated_amr <- amr_ellipses_coordinates$ellipses_coord$common_variants
names(separated_amr) <- amr_ellipses_coordinates$ellipses_coord$abbreviation
separated_amr_euler <- eulerr::euler(separated_amr, shape="ellipse")
separated_ellipses_coordinates <- separated_amr_euler$ellipses
separated_ellipses_coordinates$abbreviation <- row.names(separated_ellipses_coordinates)
separated_ellipses_coordinates <- merge(separated_ellipses_coordinates, amr_ellipses_coordinates$ellipses_coord[,-which(names(amr_ellipses_coordinates$ellipses_coord) %in% c("h", "k", "a", "b", "phi"))], by="abbreviation")
separated_ellipses_coordinates <- separated_ellipses_coordinates[order(separated_ellipses_coordinates$common_variants),]
heptagram <- data.frame(
  x=c(4000*cos(2*pi/7), 4000*cos(2*pi/7), 4000*cos(6*pi/7), 4000*cos(6*pi/7), 4000*cos(4*pi/7), 4000*cos(4*pi/7), 4000),
  y=c(-4000*sin(2*pi/7), 4000*sin(2*pi/7), 4000*sin(6*pi/7), -4000*sin(6*pi/7), 4000*sin(4*pi/7), -4000*sin(4*pi/7), 0)
)

octagon <- data.frame(
  x=c(0, 4000, 0, -4000, 4000/sqrt(2), 4000/sqrt(2), -4000/sqrt(2), -4000/sqrt(2)),
  y=c(-4000, 0, 4000, 0, 4000/sqrt(2), -4000/sqrt(2), -4000/sqrt(2), 4000/sqrt(2))
)

separated_ellipses_coordinates$h <- octagon$x
separated_ellipses_coordinates$k <- octagon$y
toJSON(separated_ellipses_coordinates)


# Figure 4B
amr_ellipses_coordinates <- generate_d3_euler(pops=AMR)
toJSON(amr_ellipses_coordinates$ellipses_coord)



# Figure 6
common_eulers <- lapply(AMR, FUN=generate_d3_euler, pops=AMR)

toJSON(common_eulers[[1]]$ellipses_coord) #ACB
toJSON(common_eulers[[2]]$ellipses_coord) #ASW
toJSON(common_eulers[[3]]$ellipses_coord) #CEU
toJSON(common_eulers[[4]]$ellipses_coord) #CLM
toJSON(common_eulers[[5]]$ellipses_coord) #GIH
toJSON(common_eulers[[6]]$ellipses_coord) #MXL
toJSON(common_eulers[[7]]$ellipses_coord) #PEL
toJSON(common_eulers[[8]]$ellipses_coord) #PUR

uniquely_common <- c()
for (e in 1:length(common_eulers)) {
  uniquely_common <- c(uniquely_common, common_eulers[[e]]$ellipses_coord[e, "unshared_common_variants"] / common_eulers[[e]]$ellipses_coord[e, "common_variants"])
}

round(uniquely_common * 100, 1)
