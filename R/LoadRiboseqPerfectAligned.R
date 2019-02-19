# Load perfect alignment between Ribo-Seq reads and genes and count reads
LoadRiboseqPerfectAligned <- function(bam, tx2gn, subset=NA, output=NA) {
  # bam       BAM file location
  # tx2gn     Name character vector of transcript (character names) to gene (character values) mapping
  # subset    If not NA, limit reads to a subset. A character vector of read IDS or a file with the IDs (one ID per line)

  require(ShortRead);
  require(GenomicRanges);
  require(GenomicAlignments);

  all <- readGAlignments(bam, use.names = TRUE);

  nms <- sort(unique(names(all)));
  ids <- 1:length(nms);
  names(ids) <- nms;

  ## Limit read ID to a subset
  if (!identical(subset, NA)) {
    if (file.exists(subset)) sub <- readLines(subset) else sub <- subset;
  } else sub <- NA;
  sub <- sub[!is.na(sub)];
  sub <- sub[sub %in% nms];

  tid <- sort(unique(names(tx2gn)));
  gid <- sort(unique(tx2gn));

  all <- list(sense=all[strand(all)=='+'], antisense=all[strand(all)=='-']);
  res <- lapply(all, function(aln) {
    rds <- as.vector(ids[names(aln)]);
    txs <- as.vector(seqnames(aln));
    gns <- as.vector(tx2gn[txs]);

    if (length(sub) > 0) {
      whc <- which(names(aln) %in% sub);
      rds <- rds[whc];
      txs <- txs[whc];
      gns <- gns[whc];
    };

    # read to transcript
    r2t <- sapply(split(txs, rds), unique);
    nm0 <- sapply(r2t, length);
    uni <- as.integer(names(nm0[nm0==1]));
    ind <- rds %in% uni;
    mp1 <- lapply(split(rds, txs), unique);
    mp2 <- lapply(split(rds[ind], txs[ind]), unique);
    nm1 <- sapply(mp1, length);
    nm2 <- sapply(mp2, length);

    ctx <- matrix(0, nrow = length(tid), ncol = 2, dimnames = list(tid, c('Total', 'Unique')));
    ctx[names(nm1), 1] <- nm1;
    ctx[names(nm2), 2] <- nm2;

    all.tx <- list(count=ctx, mapping=list(read2tx=r2t, tx2read=mp1));

    # read to gene
    r2g <- lapply(split(gns, rds), unique);
    nm0 <- sapply(r2g, length);
    uni <- as.integer(names(nm0[nm0==1]));
    ind <- rds %in% uni;
    mp1 <- lapply(split(rds, gns), unique);
    mp2 <- lapply(split(rds[ind], gns[ind]), unique);
    nm1 <- sapply(mp1, length);
    nm2 <- sapply(mp2, length);

    cgn <- matrix(0, nrow = length(gid), ncol = 2, dimnames = list(gid, c('Total', 'Unique')));
    cgn[names(nm1), 1] <- nm1;
    cgn[names(nm2), 2] <- nm2;

    all.gn <- list(count=cgn, mapping=list(read2gene=r2g, gene2read=mp1));

    list(tx=all.tx, gene=all.gn);
  });

  # transcript and gene level read counts
  cnt1 <- cbind(res$sense$tx$count, res$antisense$tx$count);
  cnt2 <- cbind(res$sense$gene$count, res$antisense$gene$count);
  colnames(cnt1) <- colnames(cnt2) <- c('Total', 'Unique', 'Total_AS', 'Unique_AS');
  cnt <- list(transcript=cnt1, gene=cnt2);

  mpp <- list(transcript=list(sense=res$sense$tx$mapping, antisense=res$antisense$tx$mapping),
              gene=list(sense=res$sense$gene$mapping, antisense=res$antisense$gene$mapping));

  if (!identical(output, NA)) {
    saveRDS(cnt, paste0(output, '_count.rds'));
    saveRDS(mpp, paste0(output, '_mapping.rds'));
    saveRDS(ids, paste0(output, '_index.rds'));
    saveRDS(all, paste0(output, '_aligned.rds'));
  };

  list(count=cnt, mapping=mpp, index=ids, subset=sub, aligned=all);
}
