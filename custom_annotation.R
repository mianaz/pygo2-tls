suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(Seurat)
  library(limma)
  library(msigdbr)
  library(fgsea)
  library(pbapply)
  library(dplyr)
  library(clusterProfiler)
  library(parallel)
})

##### API functions #####
call_llm_api <- function(base_url, model_name, system_prompt, user_prompt,
                         api_key, track_tokens = TRUE) {
  # Validate required parameters
  if (missing(base_url) || missing(model_name) || missing(api_key)) {
    stop("base_url, model_name, and api_key are required parameters")
  }

  headers <- c(
    "Content-Type" = "application/json",
    "Authorization" = paste("Bearer", api_key)
  )

  body <- list(
    model = model_name,
    messages = list(
      list(role = "system", content = system_prompt),
      list(role = "user", content = user_prompt)
    )
  )

  # Simple retry
  max_retries <- 3
  response <- NULL
  for (attempt in 1:max_retries) {
    response <- tryCatch({
      httr::POST(
        url = base_url,
        add_headers(.headers = headers),
        body = toJSON(body, auto_unbox = TRUE),
        encode = "json",
        timeout(30)
      )
    }, error = function(e) {
      if (attempt < max_retries) {
        message(paste("API request failed, retrying... (attempt", attempt, "of", max_retries, ")"))
        Sys.sleep(2^attempt)
        return(NULL)
      } else {
        stop("API request failed after ", max_retries, " attempts: ", e$message)
      }
    })
    if (!is.null(response)) break
  }

  if (http_error(response)) {
    error_content <- content(response, "text", encoding = "UTF-8")
    status <- status_code(response)
    error_msg <- tryCatch({
      error_json <- fromJSON(error_content)
      if (!is.null(error_json$error$message)) error_json$error$message else error_content
    }, error = function(e) error_content)
    stop(sprintf("API request failed [%d]: %s", status, error_msg))
  }

  parsed <- fromJSON(rawToChar(response$content), simplifyVector = FALSE)
  # Extract content robustly across providers
  result_content <- NULL
  if (!is.null(parsed$choices) && length(parsed$choices) > 0) {
    choice <- parsed$choices[[1]]
    if (!is.null(choice$message) && is.list(choice$message) && !is.null(choice$message$content)) {
      result_content <- choice$message$content
    } else if (!is.null(choice$text)) {
      result_content <- choice$text
    } else if (!is.null(choice$message) && is.character(choice$message)) {
      result_content <- choice$message
    }
  }
  if (is.null(result_content) && !is.null(parsed$output_text)) {
    result_content <- parsed$output_text
  }
  if (is.null(result_content) && !is.null(parsed$content)) {
    if (is.character(parsed$content)) {
      result_content <- parsed$content
    } else if (is.list(parsed$content) && length(parsed$content) > 0) {
      # Concatenate text parts if present (OpenAI responses-style)
      parts <- vapply(parsed$content, function(p) {
        if (is.list(p) && !is.null(p$text)) return(as.character(p$text))
        ""
      }, character(1))
      if (any(nzchar(parts))) result_content <- paste(parts[nzchar(parts)], collapse = "\n")
    }
  }
  if (is.null(result_content)) {
    stop("Unexpected API response format; could not extract content")
  }
  result_usage <- if (!is.null(parsed$usage)) parsed$usage else NULL

  if (track_tokens && !is.null(result_usage)) {
    message(sprintf("Tokens used - Prompt: %d, Completion: %d, Total: %d",
                    result_usage$prompt_tokens,
                    result_usage$completion_tokens,
                    result_usage$total_tokens))
  }

  return(result_content)
}

##### Main function #####
run_aiCelltype <- function(seurat_obj,
                          de_markers = NULL,
                          species = c("human", "mouse"),
                          tissue,
                          base_url = "https://api.openai.com/v1/chat/completions",
                          model_name = "gpt-4o",
                          api_key = NULL,
                          top_n = 10,
                          gsea_data = NULL,
                          top_gs = 5,
                          verbose = TRUE) {
  species <- match.arg(species)
  if (missing(tissue) || is.null(tissue)) stop("Tissue must be specified.")
  if (is.null(api_key)) {
    api_key <- Sys.getenv("OPENAI_API_KEY")
    if (!nzchar(api_key)) stop("API key not provided. Set api_key or OPENAI_API_KEY env var.")
  }
  if (verbose) message(paste("Using model:", model_name))

  # Prepare DE markers if not provided
  if (is.null(de_markers)) {
    if (verbose) message("Cluster markers not provided. Calculating cluster markers...")
    de_markers <- FindAllMarkers(seurat_obj, only.pos = FALSE,
                                 min.pct = 0.5,
                                 logfc.threshold = 0.25,
                                 return.thresh = 0.01)
  }

  # Correct gene symbols
  if (verbose) message("Correcting gene symbols...")
  unique_genes <- unique(de_markers$gene)
  check <- checkGeneSymbols(unique_genes, species = species)
  valid_genes <- check$Suggested.Symbol[!is.na(check$Suggested.Symbol)]
  de_markers <- de_markers[de_markers$gene %in% valid_genes, ]
  if (nrow(de_markers) == 0) stop("No valid genes remain after symbol correction.")

  # Prepare per-cluster data
  clusters <- unique(de_markers$cluster)
  combined_data <- lapply(clusters, function(cluster) {
    idx <- de_markers$cluster == cluster & de_markers$p_val_adj < 0.05 & !grepl("Rik$|^Gm", de_markers$gene)
    cluster_de <- de_markers[idx, , drop = FALSE]
    top_pos <- utils::head(cluster_de[cluster_de$avg_log2FC > 0, , drop = FALSE], top_n)
    top_neg <- utils::head(cluster_de[cluster_de$avg_log2FC < 0, , drop = FALSE], top_n)
    # Optional enriched pathways from supplied gsea_data
    enriched <- character(0)
    if (!is.null(gsea_data)) {
      key <- as.character(cluster)
      pathways_obj <- NULL
      if (is.list(gsea_data)) {
        if (!is.null(names(gsea_data)) && key %in% names(gsea_data)) {
          pathways_obj <- gsea_data[[key]]
        } else if (length(gsea_data) >= as.integer(cluster)) {
          pathways_obj <- gsea_data[[as.integer(cluster)]]
        }
      }
      if (!is.null(pathways_obj)) {
        if (is.data.frame(pathways_obj)) {
          df <- pathways_obj
          if ("padj" %in% names(df)) df <- df[df$padj < 0.05, , drop = FALSE]
          if ("NES" %in% names(df)) df <- df[order(-abs(df$NES)), , drop = FALSE]
          if (nrow(df) > 0) {
            if ("pathway" %in% names(df)) {
              enriched <- utils::head(df$pathway, top_gs)
            } else {
              enriched <- utils::head(df[[1]], top_gs)
            }
          }
        } else if (is.character(pathways_obj)) {
          enriched <- utils::head(pathways_obj, top_gs)
        }
      }
    }
    list(
      positive_markers = top_pos$gene,
      negative_markers = top_neg$gene,
      enriched_pathways = enriched
    )
  })
  names(combined_data) <- as.character(clusters)

  # Prompts
  system_prompt <- "You're an expert cell biologist. Predict cell types from markers and gsea results.\nOutput format: 'Cluster X: Cell_Type (0.xx)'\nConfidence: 0-0.2=uncertain, 0.2-0.4=low, 0.4-0.6=moderate, 0.6-0.8=high, 0.8-1.0=very high\nBe specific with cell subtypes when possible. Never predict 'unknown'."

  build_prompt <- function(cluster_ids, combined_data) {
    batch_prompt <- paste(
      paste("Species:", species),
      paste("Tissue:", tissue),
      "Predict cell types for the following clusters:",
      sep = "\n"
    )
    for (cl in cluster_ids) {
      batch_prompt <- paste(
        batch_prompt,
        paste0("\nCluster ", cl, ":"),
        paste("Positive markers:", paste(combined_data[[cl]]$positive_markers, collapse = ", ")),
        paste("Negative markers:", paste(combined_data[[cl]]$negative_markers, collapse = ", ")),
        if (!is.null(combined_data[[cl]]$enriched_pathways) && length(combined_data[[cl]]$enriched_pathways) > 0) {
          paste("Enriched pathways:", paste(combined_data[[cl]]$enriched_pathways, collapse = ", "))
        },
        sep = "\n"
      )
    }
    paste(batch_prompt,
          "\nProvide one line per cluster. Example: Cluster 0: CD8+ T cells (0.85)",
          sep = "\n")
  }

  cluster_batches <- list(names(combined_data))

  # Single prediction run
  predictions <- data.frame(
    cluster = character(),
    predicted_label = character(),
    llm_confidence = numeric(),
    stringsAsFactors = FALSE
  )

  for (batch in cluster_batches) {
    batch_prompt <- build_prompt(batch, combined_data)
    api_result <- tryCatch({
      call_llm_api(base_url, model_name, system_prompt, batch_prompt,
                   api_key = api_key)
    }, error = function(e) stop("API call failed: ", e$message))

    lines <- strsplit(api_result, "\n")[[1]]
    for (ln in lines) {
      match <- regexec("^Cluster (\\d+): (.+?)\\s*\\((0?\\.[0-9]+|1\\.0?)\\)\\s*$", ln, perl = TRUE)
      if (match[[1]][1] != -1) {
        cluster <- substr(ln, match[[1]][2], match[[1]][2] + attr(match[[1]], "match.length")[2] - 1)
        cell_type <- substr(ln, match[[1]][3], match[[1]][3] + attr(match[[1]], "match.length")[3] - 1)
        confidence <- as.numeric(substr(ln, match[[1]][4], match[[1]][4] + attr(match[[1]], "match.length")[4] - 1))
        predictions <- rbind(predictions,
                             data.frame(cluster = cluster,
                                        predicted_label = cell_type,
                                        llm_confidence = confidence))
      }
    }
  }

  if (verbose) {
    message("Cell type predictions completed:")
    print(predictions[, c("cluster", "predicted_label", "llm_confidence")])
  }
  return(predictions)
}

##### GSEA helper #####
run_gsea <- function(de_markers,
                     species = c("human", "mouse"),
                     method = "fgsea",
                     gene_set = "C8",
                     min_size = 15,
                     max_size = 500,
                     nperm = 1000,
                     pvalue_cutoff = 0.05,
                     ncores = detectCores() - 1,
                     verbose = TRUE) {

  species <- match.arg(species)
  method <- match.arg(method, c("fgsea", "clusterProfiler"))
  if (verbose) message(paste("Performing GSEA using", method, "..."))

  # Process gene sets
  if (is.character(gene_set)) {
    if (verbose) message(paste("Using msigdbr gene set categories:", paste(gene_set, collapse = ", ")))
    # Combine multiple gene set collections
    m_df <- do.call(rbind, lapply(gene_set, function(cat) {
      if (species == "mouse") {
        msigdbr(db_species = "MM",  # Mouse database
                species = "Mus musculus",
                collection = cat)
      } else {
        msigdbr(db_species = "HS",  # Human database
                species = "Homo sapiens",
                collection = cat)
      }
    }))
    m_list <- split(m_df$gene_symbol, m_df$gs_name)
  } else if (is.list(gene_set)) {
    if (verbose) message("Using custom gene set list")
    m_list <- gene_set
    m_df <- data.frame(
      gs_name = rep(names(m_list), sapply(m_list, length)),
      gene_symbol = unlist(m_list)
    )
  } else {
    stop("gene_set must be either a character vector specifying msigdbr categories or a named list of gene sets")
  }

  # Perform GSEA
  enrich_results <- pblapply(unique(de_markers$cluster), function(cluster) {
    cluster_genes <- de_markers[de_markers$cluster == cluster, ]
    ranks <- setNames(cluster_genes$avg_log2FC, cluster_genes$gene)
    if (length(ranks) < min_size) {
      if (verbose) message(paste("Skipping GSEA for cluster", cluster, ": insufficient genes"))
      return(NULL)
    }
    tryCatch({
      if (method == "fgsea") {
        result <- fgsea(pathways = m_list,
                        stats = ranks,
                        minSize = min_size,
                        maxSize = max_size,
                        nperm = nperm)
        result <- result %>%
          as.data.frame() %>%
          dplyr::filter(padj < pvalue_cutoff) %>%
          dplyr::arrange(desc(abs(NES)))

        # Add collection information to results
        result$collection <- sapply(result$pathway, function(p) {
          collection_info <- m_df[m_df$gs_name == p, "gs_cat", drop = TRUE][1]
          return(collection_info)
        })

      } else if (method == "clusterProfiler") {
        gene_list <- sort(ranks, decreasing = TRUE)
        gse <- GSEA(gene_list, TERM2GENE = m_df,
                    minGSSize = min_size,
                    maxGSSize = max_size,
                    pvalueCutoff = pvalue_cutoff,
                    pAdjustMethod = "BH",
                    nPermSimple = nperm)
        result <- as.data.frame(gse) %>%
          dplyr::rename(pathway = ID, pval = pvalue, padj = p.adjust, ES = enrichmentScore) %>%
          dplyr::select(pathway, pval, padj, ES, NES, everything()) %>%
          dplyr::arrange(desc(abs(NES)))

        # Add collection information to results
        result$collection <- sapply(result$pathway, function(p) {
          collection_info <- m_df[m_df$gs_name == p, "gs_cat", drop = TRUE][1]
          return(collection_info)
        })
      }
      if (verbose) message(paste("Completed GSEA for cluster", cluster))
      return(result)
    }, error = function(e) {
      warning(paste("GSEA failed for cluster", cluster, ":", e$message))
      return(NULL)
    })
  }, cl = ncores)

  names(enrich_results) <- unique(de_markers$cluster)
  enrich_results <- enrich_results[!sapply(enrich_results, is.null)]
  return(enrich_results)
}
