rm(list = ls())
graphics.off()

library(data.table)
library(dplyr)
library(ggplot2)
library(grid)

# Input

input_dir <- "snptest"
output_dir <- "snptest_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

models <- list(
  Model1 = list(
    pattern = "^chr[0-9]+_01_sex_maternal_age\\.snptest(\\.gz)?$",
    title = "SNPTEST model 1: Sex + maternal age"
  ),
  Model2 = list(
    pattern = "^chr[0-9]+_02_weight_category_exact_ga\\.snptest(\\.gz)?$",
    title = "SNPTEST model 2: Weight category + gestational age"
  ),
  Model3 = list(
    pattern = "^chr[0-9]+_03_gd_pre_spb_latency\\.snptest(\\.gz)?$",
    title = "SNPTEST model 3: Gestational diabetes/preeclampsia + SPB/PPROM + latency"
  ),
  Model4 = list(
    pattern = "^chr[0-9]+_04_mothers_age_nr_preg\\.snptest(\\.gz)?$",
    title = "SNPTEST model 4: Mother's age + number of pregnancies"
  ),
  Model5 = list(
    pattern = "^chr[0-9]+_05_prev_preg_mother_ancestry\\.snptest(\\.gz)?$",
    title = "SNPTEST model 5: Previous pregnancy complications + maternal ancestry"
  ),
  Model6 = list(
    pattern = "^chr[0-9]+_06_w_ga_weight_sex\\.snptest(\\.gz)?$",
    title = "SNPTEST model 6: Gestational age + birth weight + sex"
  )
)


genomewide_p <- 5e-8
suggestive_p  <- 5e-5

genomewide_line <- -log10(genomewide_p)
suggestive_line <- -log10(suggestive_p)

find_col <- function(df, candidates) {
  nm <- names(df)
  hit <- nm[tolower(nm) %in% tolower(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

to_numeric_clean <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", ".", "NA", "NaN", "nan", "Inf", "-Inf")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

to_chr_clean <- function(chr_value, source_file = NULL) {
  chr_txt <- trimws(as.character(chr_value))
  chr_txt <- gsub("^chr", "", chr_txt, ignore.case = TRUE)
  chr_num <- suppressWarnings(as.numeric(chr_txt))

  if (!is.null(source_file)) {
    src_txt <- as.character(source_file)
    chr_from_file <- suppressWarnings(
      as.numeric(sub("^.*chr([0-9]+).*$", "\\1", src_txt, ignore.case = TRUE))
    )
    chr_num[is.na(chr_num) & !is.na(chr_from_file)] <- chr_from_file[is.na(chr_num) & !is.na(chr_from_file)]
  }

  chr_num
}

choose_best_p_col <- function(df, candidates) {
  nm <- names(df)
  hits <- nm[tolower(nm) %in% tolower(candidates)]
  if (length(hits) == 0) return(NA_character_)

  valid_counts <- sapply(hits, function(col) {
    p <- to_numeric_clean(df[[col]])
    sum(!is.na(p) & p > 0 & p <= 1)
  })

  if (all(valid_counts == 0)) return(hits[1])
  hits[which.max(valid_counts)]
}

read_one_snptest_file <- function(file) {
  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) base::gzfile(file, "rt") else base::file(file, "rt")
  lines <- readLines(con, n = 1000, warn = FALSE)
  close(con)

  header_line <- grep(
    "^(alternate_ids|rsid|chromosome|SNP|snp_id|ID)[[:space:]]",
    lines,
    ignore.case = TRUE
  )[1]

  if (is.na(header_line)) {
    stop("Could not find SNPTEST results header in file: ", file)
  }

  header_text <- lines[header_line]

  df <- fread(
    file,
    header = TRUE,
    skip = header_text,
    data.table = FALSE,
    fill = TRUE,
    na.strings = c("NA", "NaN", "nan", ".", "-9"),
    check.names = FALSE,
    showProgress = FALSE
  )

  first_col <- names(df)[1]
  df <- df[!grepl("^(#|Analysis:|SNPTEST)", as.character(df[[first_col]]), ignore.case = TRUE), , drop = FALSE]
  df <- df[as.character(df[[first_col]]) != first_col, , drop = FALSE]

  df
}

read_model_files <- function(pattern, input_dir = ".") {
  files <- list.files(path = input_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop("No files matched pattern in folder '", input_dir, "': ", pattern)
  }

  message("  Matched ", length(files), " files")

  df <- rbindlist(
    lapply(files, function(f) {
      x <- read_one_snptest_file(f)
      x$source_file <- basename(f)
      x
    }),
    fill = TRUE
  )

  as.data.frame(df)
}

standardise_snptest <- function(df) {
  snp_col <- find_col(df, c("rsid", "SNP", "snp", "rs_id", "alternate_ids", "alternate_id", "id", "variant", "ID"))
  chr_col <- find_col(df, c("chromosome", "CHR", "chrom", "chr"))
  bp_col  <- find_col(df, c("position", "BP", "pos", "base_pair_location"))
  p_col   <- choose_best_p_col(df, c(
    "frequentist_add_pvalue",
    "frequentist_score_pvalue",
    "frequentist_1_pvalue",
    "frequentist_pvalue",
    "pvalue",
    "p_value",
    "pval",
    "p.value",
    "P",
    "p"
  ))

  if (is.na(snp_col) || is.na(chr_col) || is.na(bp_col) || is.na(p_col)) {
    message("Available columns are:")
    message(paste(names(df), collapse = ", "))
    stop("Could not identify required SNPTEST columns. Need SNP ID, chromosome, position, and p-value.")
  }

  message("  Using columns: SNP=", snp_col, ", CHR=", chr_col, ", BP=", bp_col, ", P=", p_col)

  source_col <- if ("source_file" %in% names(df)) df[["source_file"]] else NULL

  out0 <- data.frame(
    SNP = as.character(df[[snp_col]]),
    CHR = to_chr_clean(df[[chr_col]], source_col),
    BP  = to_numeric_clean(df[[bp_col]]),
    P   = to_numeric_clean(df[[p_col]]),
    stringsAsFactors = FALSE
  )

  message("  Rows read: ", nrow(out0))
  message("  Rows with valid chromosome: ", sum(!is.na(out0$CHR) & out0$CHR %in% 1:22))
  message("  Rows with valid position: ", sum(!is.na(out0$BP) & out0$BP > 0))
  message("  Rows with valid p-value: ", sum(!is.na(out0$P) & out0$P > 0 & out0$P <= 1))

  if (sum(!is.na(out0$P) & out0$P > 0 & out0$P <= 1) == 0) {
    p_cols <- names(df)[grepl("pvalue|p_value|pval|p.value|^p$|^P$", names(df), ignore.case = TRUE)]
    if (length(p_cols) > 0) {
      message("  Candidate p-value columns and number of valid p-values:")
      for (pc in p_cols) {
        pv <- to_numeric_clean(df[[pc]])
        message("    ", pc, ": ", sum(!is.na(pv) & pv > 0 & pv <= 1))
      }
    }
    message("  First 5 values in chosen p-value column:")
    message(paste(head(as.character(df[[p_col]]), 5), collapse = ", "))
  }

  out <- out0 %>%
    filter(!is.na(CHR), CHR %in% 1:22) %>%
    filter(!is.na(BP), BP > 0) %>%
    filter(!is.na(P), P > 0, P <= 1) %>%
    arrange(CHR, BP)

  if (nrow(out) == 0) {
    stop("No valid SNPTEST rows after filtering. Check diagnostic counts above, especially valid p-values.")
  }

  out
}

prepare_manhattan_data <- function(df) {
  chr_offsets <- df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(offset = cumsum(chr_len) - chr_len)

  plot_df <- df %>%
    left_join(chr_offsets, by = "CHR") %>%
    mutate(
      BPcum = BP + offset,
      logP = -log10(P)
    )

  axis_df <- plot_df %>%
    group_by(CHR) %>%
    summarise(center = (min(BPcum) + max(BPcum)) / 2, .groups = "drop")

  list(data = plot_df, axis = axis_df)
}

prepare_qq_data <- function(df) {
  qq_df <- df %>%
    filter(!is.na(P), P > 0, P <= 1) %>%
    arrange(P)

  if (nrow(qq_df) == 0) {
    stop("No valid p-values for QQ plot")
  }

  n <- nrow(qq_df)

  qq_df <- qq_df %>%
    mutate(
      expected = -log10(ppoints(n)),
      observed = -log10(P)
    )

  lambda <- median(qchisq(1 - qq_df$P, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)

  list(data = qq_df, lambda = lambda)
}

make_manhattan_plot <- function(df, model_title, y_max_common = NULL) {
  prepared <- prepare_manhattan_data(df)
  plot_df <- prepared$data
  axis_df <- prepared$axis

  y_max <- if (is.null(y_max_common)) {
    max(9, ceiling(max(plot_df$logP, suggestive_line, genomewide_line, na.rm = TRUE) + 1))
  } else {
    y_max_common
  }

  ggplot(plot_df, aes(x = BPcum, y = logP)) +
    geom_point(
      aes(color = factor(CHR %% 2)),
      size = 0.45,
      alpha = 0.85,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("grey35", "grey70")) +
    geom_hline(
      aes(yintercept = genomewide_line, linetype = "Genome-wide"),
      color = "red",
      linewidth = 0.45
    ) +
    geom_hline(
      aes(yintercept = suggestive_line, linetype = "Suggestive"),
      color = "blue",
      linewidth = 0.45
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("Genome-wide" = "solid", "Suggestive" = "dashed"),
      labels = c(
        "Genome-wide" = expression("Genome-wide significance (" * italic(p) < 5 %*% 10^-8 * ")"),
        "Suggestive"  = expression("Suggestive significance (" * italic(p) < 5 %*% 10^-5 * ")")
      )
    ) +
    guides(linetype = guide_legend(override.aes = list(color = c("red", "blue")))) +
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df$CHR,
      expand = expansion(mult = c(0.005, 0.01))
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, y_max, by = 2),
      expand = c(0, 0)
    ) +
    labs(
      title = model_title,
      x = "Chromosome",
      y = expression(-log[10](italic(p))),
      linetype = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 7)),
      axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 6)),
      axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 6)),
      axis.text.x = element_text(size = 8, margin = margin(t = 2)),
      axis.text.y = element_text(size = 8),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.35),
      axis.ticks.length = unit(0.12, "cm"),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.text = element_text(size = 7.8),
      legend.key.width = unit(1.1, "cm"),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.35),
      legend.margin = margin(t = 2, r = 4, b = 2, l = 4),
      plot.margin = margin(10, 10, 8, 8)
    )
}

make_qq_plot <- function(df, model_title, max_axis_common = NULL) {
  prepared <- prepare_qq_data(df)
  qq_df <- prepared$data
  lambda <- prepared$lambda

  max_axis <- if (is.null(max_axis_common)) {
    ceiling(max(c(qq_df$expected, qq_df$observed), na.rm = TRUE) * 10) / 10
  } else {
    max_axis_common
  }

  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(color = "grey35", size = 0.45, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.45) +
    annotate(
      "text",
      x = max_axis * 0.10,
      y = max_axis * 0.92,
      label = paste0("lambda = ", sprintf("%.3f", lambda)),
      hjust = 0,
      size = 3.0
    ) +
    scale_x_continuous(limits = c(0, max_axis), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, max_axis), expand = c(0, 0)) +
    labs(
      title = model_title,
      x = expression(Expected~~-log[10](italic(p))),
      y = expression(Observed~~-log[10](italic(p)))
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 7)),
      axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 6)),
      axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 6)),
      axis.text = element_text(size = 8),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.35),
      axis.ticks.length = unit(0.12, "cm"),
      plot.margin = margin(10, 10, 8, 8)
    )
}

save_plot_grid <- function(plot_list, filename, ncol = 2, width = 16, height = 18, dpi = 300) {
  if (length(plot_list) == 0) {
    stop("No plots available to save: ", filename)
  }

  nrow <- ceiling(length(plot_list) / ncol)

  png(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow, ncol)))

  for (i in seq_along(plot_list)) {
    row_i <- ceiling(i / ncol)
    col_i <- ((i - 1) %% ncol) + 1
    pushViewport(viewport(layout.pos.row = row_i, layout.pos.col = col_i))
    grid.draw(ggplotGrob(plot_list[[i]]))
    popViewport()
  }

  popViewport()
  dev.off()
}

clean_data <- list()

for (model_name in names(models)) {
  message("Processing ", model_name)

  raw_df <- tryCatch(
    read_model_files(models[[model_name]]$pattern, input_dir),
    error = function(e) {
      message("  Skipping ", model_name, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(raw_df)) next

  df <- tryCatch(
    standardise_snptest(raw_df),
    error = function(e) {
      message("  Skipping ", model_name, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(df)) next

  message("  Valid SNPs: ", nrow(df))
  message("  Minimum p-value: ", signif(min(df$P, na.rm = TRUE), 4))

  fwrite(
    df,
    file.path(output_dir, paste0(model_name, "_snptest_clean_for_plotting.tsv")),
    sep = "\t"
  )

  clean_data[[model_name]] <- df
}

if (length(clean_data) == 0) {
  stop("No valid SNPTEST datasets were available for plotting.")
}

all_logp_max <- max(unlist(lapply(clean_data, function(x) -log10(x$P))), suggestive_line, genomewide_line, na.rm = TRUE)
y_max_common <- max(9, ceiling(all_logp_max + 1))

all_expected_observed_max <- max(unlist(lapply(clean_data, function(x) {
  qq <- prepare_qq_data(x)$data
  c(qq$expected, qq$observed)
})), na.rm = TRUE)
max_axis_common <- ceiling(all_expected_observed_max * 10) / 10

manhattan_plots <- list()
qq_plots <- list()

for (model_name in names(clean_data)) {
  manhattan_plots[[model_name]] <- make_manhattan_plot(
    clean_data[[model_name]],
    models[[model_name]]$title,
    y_max_common = y_max_common
  )

  qq_plots[[model_name]] <- make_qq_plot(
    clean_data[[model_name]],
    models[[model_name]]$title,
    max_axis_common = max_axis_common
  )
}

out_manhattan <- file.path(output_dir, "combined_snptest_manhattan.png")
out_qq <- file.path(output_dir, "combined_snptest_qq.png")

save_plot_grid(
  manhattan_plots,
  filename = out_manhattan,
  ncol = 2,
  width = 18,
  height = 18,
  dpi = 300
)
message("Saved: ", out_manhattan)

save_plot_grid(
  qq_plots,
  filename = out_qq,
  ncol = 2,
  width = 14,
  height = 18,
  dpi = 300
)
message("Saved: ", out_qq)

message("Done.")
