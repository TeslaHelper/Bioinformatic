# ==============================================================================
# Analysis 4: Background Universe Correction for GO/KEGG Enrichment
# Re-implements the clusterProfiler hypergeometric test (phyper-based) with a
# restricted universe (205-gene bridge network) instead of the whole proteome,
# using real Bioconductor annotation data (org.Hs.eg.db, GO.db).
# ==============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
})

hub_genes <- read.csv("hub_genes.csv", stringsAsFactors = FALSE)$Gene          # n = 15, the query
universe_205 <- read.csv("bridge_205_genes.csv", stringsAsFactors = FALSE)$Gene # n = 205, restricted universe

cat("Query set (hub genes):", length(hub_genes), "genes\n")
cat("Restricted universe (bridge network):", length(universe_205), "genes\n\n")

# ------------------------------------------------------------------------------
# Map gene symbols -> Entrez IDs
# ------------------------------------------------------------------------------
map_symbols <- function(symbols) {
  m <- suppressWarnings(AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", columns = "ENTREZID"))
  m <- m[!is.na(m$ENTREZID), ]
  m <- m[!duplicated(m$SYMBOL), ]
  m
}

hub_map <- map_symbols(hub_genes)
universe_map <- map_symbols(universe_205)

cat("Hub genes mapped to Entrez:", nrow(hub_map), "/", length(hub_genes), "\n")
cat("Universe genes mapped to Entrez:", nrow(universe_map), "/", length(universe_205), "\n")
missing_hub <- setdiff(hub_genes, hub_map$SYMBOL)
if (length(missing_hub) > 0) cat("  [Unmapped hub genes -- excluded]:", paste(missing_hub, collapse=", "), "\n")
missing_universe <- setdiff(universe_205, universe_map$SYMBOL)
if (length(missing_universe) > 0) cat("  [Unmapped universe genes -- excluded]:", length(missing_universe), "genes\n")
cat("\n")

hub_entrez <- unique(hub_map$ENTREZID)
universe_entrez <- unique(universe_map$ENTREZID)
# Ensure query is a subset of universe (required for a valid hypergeometric test)
universe_entrez <- union(universe_entrez, hub_entrez)
N_universe <- length(universe_entrez)
n_query <- length(hub_entrez)
cat("Final N (restricted universe, Entrez-mapped, union with query):", N_universe, "\n")
cat("Final n (query, Entrez-mapped):", n_query, "\n\n")

# ------------------------------------------------------------------------------
# GO annotation: gene -> GO terms (Biological Process only, matching original focus)
# ------------------------------------------------------------------------------
cat("Retrieving GO (BP) annotations for restricted universe (propagated GOALL, matching clusterProfiler's true-path-rule behavior)...\n")
go_map <- suppressWarnings(AnnotationDbi::select(org.Hs.eg.db, keys = universe_entrez, keytype = "ENTREZID",
                                    columns = c("GOALL", "ONTOLOGYALL")))
go_map <- go_map[go_map$ONTOLOGYALL == "BP" & !is.na(go_map$GOALL), ]
go_map <- unique(go_map[, c("ENTREZID", "GOALL")])
colnames(go_map)[2] <- "GO"

# Build term -> gene list (restricted to universe)
term2gene_universe <- split(go_map$ENTREZID, go_map$GO)
# Restrict to terms with reasonable size (matches typical clusterProfiler defaults: 3-500... but
# our universe is tiny (205), so relax upper bound to universe size)
term_sizes <- lengths(term2gene_universe)
valid_terms <- names(term_sizes)[term_sizes >= 3]
cat("GO BP terms represented in the 205-gene universe (size >=3):", length(valid_terms), "\n\n")

# ------------------------------------------------------------------------------
# Hypergeometric test: restricted universe
# phyper(x-1, K, N-K, n, lower.tail=FALSE)
#   x = query genes in term ; K = universe genes in term ; N = universe size ; n = query size
# ------------------------------------------------------------------------------
run_hypergeom <- function(term2gene, N, query_entrez, n_query, min_size = 3) {
  out <- data.frame()
  for (term in names(term2gene)) {
    genes_in_term <- term2gene[[term]]
    K <- length(genes_in_term)
    if (K < min_size) next
    x <- length(intersect(query_entrez, genes_in_term))
    if (x == 0) next
    pval <- phyper(x - 1, K, N - K, n_query, lower.tail = FALSE)
    out <- rbind(out, data.frame(Term = term, Count = x, BgCount = K,
                                   N = N, n = n_query, pvalue = pval,
                                   GeneRatio = paste0(x, "/", n_query),
                                   BgRatio = paste0(K, "/", N)))
  }
  if (nrow(out) > 0) {
    out$p.adjust <- p.adjust(out$pvalue, method = "BH")
    out <- out[order(out$pvalue), ]
  }
  out
}

cat("Running restricted-universe hypergeometric test (GO BP)...\n")
go_restricted <- run_hypergeom(term2gene_universe, N_universe, hub_entrez, n_query)

# Add human-readable GO term descriptions
go_terms_desc <- suppressWarnings(AnnotationDbi::select(GO.db, keys = go_restricted$Term,
                                                          keytype = "GOID", columns = "TERM"))
go_restricted <- merge(go_restricted, go_terms_desc, by.x = "Term", by.y = "GOID", all.x = TRUE)
go_restricted <- go_restricted[order(go_restricted$pvalue), ]

write.csv(go_restricted, "results/GO_BP_universe205_restricted.csv", row.names = FALSE)

cat("\nTop 15 GO BP terms (205-gene restricted universe):\n")
print(head(go_restricted[, c("TERM", "GeneRatio", "BgRatio", "pvalue", "p.adjust")], 15), row.names = FALSE)

# ------------------------------------------------------------------------------
# KEGG (legacy PATH mapping bundled in org.Hs.eg.db)
# ------------------------------------------------------------------------------
cat("\n\nRetrieving KEGG PATH annotations for restricted universe...\n")
path_map <- suppressWarnings(AnnotationDbi::select(org.Hs.eg.db, keys = universe_entrez, keytype = "ENTREZID", columns = "PATH"))
path_map <- path_map[!is.na(path_map$PATH), ]
path_map <- unique(path_map)
term2gene_kegg <- split(path_map$ENTREZID, path_map$PATH)
cat("KEGG pathways represented in the 205-gene universe:", length(term2gene_kegg), "\n\n")

kegg_restricted <- run_hypergeom(term2gene_kegg, N_universe, hub_entrez, n_query)
kegg_restricted <- kegg_restricted[order(kegg_restricted$pvalue), ]
write.csv(kegg_restricted, "results/KEGG_universe205_restricted.csv", row.names = FALSE)

cat("Top 15 KEGG pathways (205-gene restricted universe):\n")
print(head(kegg_restricted[, c("Term", "GeneRatio", "BgRatio", "pvalue", "p.adjust")], 15), row.names = FALSE)

cat("\n\n=================================================================\n")
cat("SUMMARY: term counts before/after universe restriction\n")
cat("=================================================================\n")
cat("GO BP terms tested (restricted universe):", nrow(go_restricted), "\n")
cat("GO BP terms significant at p.adjust<0.05 (restricted):", sum(go_restricted$p.adjust < 0.05), "\n")
cat("KEGG pathways tested (restricted universe):", nrow(kegg_restricted), "\n")
cat("KEGG pathways significant at p.adjust<0.05 (restricted):", sum(kegg_restricted$p.adjust < 0.05), "\n")
