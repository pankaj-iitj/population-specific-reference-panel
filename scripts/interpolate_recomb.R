# The R-script to interpolate recombination rates from known HapMap variants

# Install GenomicRanges if not installed
if(!require("GenomicRanges", character.only = TRUE, quietly = TRUE)){
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
    require("BiocManager", character.only = TRUE, quietly = TRUE)
    BiocManager::install("GenomicRanges")
}


## Assuming subsetted VCF files are in a subdir called 'vcf'
getPos <- function(chr){
    file <- paste0("vcf/", chr, ".vcf")
    pos <- system(paste("cat", file, "| awk '$1 !~ /#/ {print $2}'"), intern = TRUE)
    gr <- GRanges(rep(chr, length(pos)), IRanges(as.numeric(pos), width = 1))
    gr
}


## If  
decodeGr <- function(pop, chr, chrgr, try.a.sub = FALSE){
    if(!file.exists("pyrho_recomb/hg38/ACB"))
        stop(paste("You must first download the hg38 pyrho recombination maps here:\nhttps://drive.google.com/drive/folders/1Tgt_7GsDO0-o02vcYSfwqHFd3JNF6R06",
                   "\nAnd uncompress in your working directory"), call. = FALSE)
    tmp <- read.table(paste0("pyrho_recomb/hg38/", pop, "/", pop, "_recombination_map_hapmap_format_hg38_chr_", sub("chr", "", chr), ".txt"), header = TRUE)
    if(try.a.sub) tmp <- tmp[1:10,]
    ## adjust to have decode data start at zero and end at the max extent of the data
    tmp[1,2] <- 0
    if(start(chrgr)[length(chrgr)] > tmp[nrow(tmp),3])
       tmp[nrow(tmp),3] <- start(chrgr)[length(chrgr)] + 1
    gr <- GRanges(tmp[,1], IRanges(tmp[,2], width = 1), cMperMb = tmp[,3], cM = tmp[,4])
    gr
}



interpThat <- function(inds, snpgr, decodeblocks) {
    blocklst <- split(start(snpgr), inds)
    blocks <- unique(inds)
    decodeblocks <- decodeblocks[blocks]
    cms <- do.call(c, lapply(seq(along = blocklst), function(x) {
                          cMperBp <- decodeblocks$cMperMb[x]/1e6
                          if(x == 1L)
                              startcM <- 0
                          else
                              startcM <- decodeblocks$cM[x-1]
                          out <- startcM + (blocklst[[x]] - start(decodeblocks[x])) * cMperBp
                      }))
    cms
}

writeOut <- function(dat, fname){
    f <-  gzfile(fname, "w")
    write.table(dat, f, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    close(f)
}


## master function

runIt <- function(chr, pop){
    if(!file.exists(pop)) dir.create(pop)
    fname <- paste0(pop, "/", chr, ".tab.gz")
    gr <- getPos(chr)
    decgr <- decodeGr(pop, chr, gr)
    inds <- follow(gr, decgr)
    cms <- interpThat(inds, gr, decgr)
    gr$cM <- cms
    out <- as(gr, "data.frame")[,c(1,2,6)]
    writeOut(out, fname)
}


## run them all
for(i in paste0("chr", 1:22))
    for(j in c("BEB","PJL","GIH","STU","ITU"))
        runIt(i,j)
