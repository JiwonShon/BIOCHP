library(Gviz)
library(data.table)
library(rtracklayer) 
library(GenomicRanges)


### Cobalt

##=== BIOCHP-SKKUM-0001-TP0
bg <- fread("/Volumes/JiwonSSD/Work/test/BIOCHP-SKKUM-0001-TP0-A01-WGS-2OJ361.cobalt.ratio.chr8_123-129M.tsv")
bg_gr <- GRanges(
  seqnames = paste0("chr", bg$chromosome),  # add "chr" as prefix 
  ranges = IRanges(start = bg$position, end = bg$position),
  # score = bg$tumorReadDepth
  score = log2(bg$tumorReadDepth + 1)  # log2 
  
)

head(bg_gr)
##=== BIOCHP-SKKUM-0002-TP0
bg2 <- fread("/Volumes/JiwonSSD/Work/test/BIOCHP-SKKUM-0002-TP0-A01-WGS-2OJ361.cobalt.ratio.chr8_123-129M.tsv")
bg_gr2 <- GRanges(
  seqnames = paste0("chr", bg2$chromosome),  
  ranges = IRanges(start = bg2$position, end = bg2$position),
  # score = bg2$tumorReadDepth
  score = log2(bg2$tumorReadDepth + 1) 
)
head(bg_gr2)

chr <- as.character(unique(seqnames(bg_gr)))
gen <-"hg19"

axisTrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr)


dt <- DataTrack(
  range = bg_gr,
  chromosome = chr,
  genome = gen,
  name = "DM Coverage",
  type = "h",       # histogram
  # ylim = c(0, 5000), # 최대값을 직접 지정하거나, max값을 미리 산출해서 지정
  # ylim = c(0, max(mcols(bg_gr)$score, na.rm=TRUE)),
  ylim = c(0, 15),
  col = "steelblue"
)

dt2 <- DataTrack(
  range = bg_gr2,
  chromosome = chr,
  genome = gen,
  name = "HSR Coverage",
  type = "h",       # histogram
  # ylim = c(0, 3000), # 최대값을 직접 지정하거나, max값을 미리 산출해서 지정
  # ylim = c(0, max(mcols(bg_gr2)$score, na.rm=TRUE)),
  ylim = c(0, 15),
  col = "steelblue"
)


# GTF 불러오기
annots <- import("/Users/jiwon/Desktop/research/AMP_AD/code/gencode.v44.annotation.gtf")
annots <- annots[seqnames(annots) == chr & annots$type %in% c("exon", "transcript")]
annots_tr <- annots[seqnames(annots) == chr & annots$type %in% c("transcript")]
annots_ex <- annots[seqnames(annots) == chr & annots$type %in% c("exon")]
mcols(annots)$symbol <- mcols(annots)$gene_name

# 두 transcript만 선택
selected <- annots[mcols(annots)$transcript_name %in% c("MYC-201", "PVT1-221")]
selected_ex <- selected[seqnames(selected) == chr & selected$type %in% c("exon")]
# View(selected_ex)

geneTrack <- GeneRegionTrack(
  selected_ex,
  genome = gen,
  chromosome = chr,
  name = "MYC & PVT1",
  #transcriptAnnotation = "symbol",
  transcriptAnnotation = "gene_name",
  featureAnnotation = "exon_number", # <- 이 부분 추가!
  collapseTranscripts = FALSE,
  # fill = "black",
  # col = NA,
  transcriptDirection = "arrow"
)

plotTracks(
  list(itrack, axisTrack, dt, dt2, geneTrack),
  from = 128700000, to = 129000000,
  main = "chr8:126Mb-129Mb Coverage",
  transcriptAnnotation = "gene_name",
  showId = TRUE
)

png("../../results/paper/plot/colo320_2d_coblat_vlog.png", width = 600, height = 300, res = 150, pointsize = 12)
plotTracks(
  list(itrack, axisTrack, dt, dt2, geneTrack),
  from = 127600000, to = 128300000,
  sizes = c(0.1, 0.25, 0.25, 0.25, 0.2),
  # main = "chr8:128.7Mb-129Mb Coverage",
  # transcriptAnnotation = "gene_name",
  showId = TRUE,
  fontfamily.title = "Arial",
  col.title = "black",
  background.title = "white",
  col.axis = "black",
  cex.title = 0.5,      # 트랙 타이틀(축 제목 등)
  cex.main = 0.5,       # Figure 전체 제목
  cex.axis = 0.5,       # y, x축 등 눈금
  cex = 0.8,            # 전체 scale(글씨)
  cex.legend = 0.8,     # 범례 글씨 
  margin = 5,
  collapseTranscripts = TRUE
)
dev.off()



png("../../results/paper/plot/colo320_3d_pvt1.png", width = 600, height = 300, res = 150, pointsize = 12)
plotTracks(
  list(itrack, axisTrack, dt, dt2, geneTrack),
  from = 128100000, to = 128200000,
  sizes = c(0.1, 0.25, 0.25, 0.25, 0.2),
  # main = "chr8:128.7Mb-129Mb Coverage",
  # transcriptAnnotation = "gene_name",
  showId = TRUE,
  fontfamily.title = "Arial",
  col.title = "black",
  background.title = "white",
  col.axis = "black",
  cex.title = 0.5,      # 트랙 타이틀(축 제목 등)
  cex.main = 0.5,       # Figure 전체 제목
  cex.axis = 0.5,       # y, x축 등 눈금
  cex = 0.8,            # 전체 scale(글씨)
  cex.legend = 0.8,     # 범례 글씨 
  margin = 5,
  collapseTranscripts = TRUE
)
dev.off()
