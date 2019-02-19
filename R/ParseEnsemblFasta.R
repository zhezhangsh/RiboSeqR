ParseEnsemblFasta <- function(fa) {
  require(ShortRead);
  
  fa <- readFasta("cds/Homo_sapiens.GRCh38.cds.all.fa.gz");
  
  seq <- DNAStringSet(fa@sread);
  id  <- as.character(fa@id);
  
  id0 <- strsplit(id, ' description:');
  des <- sapply(id0, function(x) x[2]);
  des[is.na(des)] <- '';
  
  id1 <- sapply(id0, function(x) x[1]);
  id1 <- strsplit(id1, ' ');
  id1 <- do.call('rbind', id1);
  
  nm  <- id1[, 1];
  cls <- id1[, 2];
  loc <- strsplit(id1[, 3], ':');
  loc <- do.call('rbind', loc);
  chr <- loc[, 3];
  stt <- as.integer(loc[, 4]);
  end <- as.integer(loc[, 5]);
  str <- as.vector(c('-1'='-', '1'='+', '0'='*')[loc[, 6]]);
  
  gn  <- sub('gene:', '', id1[, 4]);
  tp1 <- sub('gene_biotype:', '', id1[, 5]);
  tp2 <- sub('transcript_biotype:', '', id1[, 6]);
  sym <- sub('gene_symbol:', '', id1[, 7]);
  
  anno <- data.frame(symbol=sym, gene=gn, chr=chr, start=stt, end=end, strand=str, class=cls, 
                     gene_type=tp1, tx_type=tp2, description=des, stringsAsFactors = FALSE);
  rownames(anno) <- names(seq) <- id1[, 1];
  
  list(anno=anno, seq=seq);
}