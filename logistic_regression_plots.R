rm(list = ls())
graphics.off()

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(patchwork)

# Input

files <- list(
  Model1 = "logistic_01_sex_maternal_age.assoc.logistic",
  Model2 = "logistic_02_weight_category_exact_ga.assoc.logistic",
  Model3 = "logistic_03_gd_pre_spb_latency.assoc.logistic",
  Model4 = "logistic_04_mothers_age_nr_preg.assoc.logistic",
  Model5 = "logistic_05_prev_preg_mother_ancestry.assoc.logistic",
  Model6 = "logistic_06_w_ga_weight_sex.assoc.logistic"
)

titles <- c(
  "Model 1: sex + maternal age",
  "Model 2: weight category + birth weight relative to GA",
  "Model 3: GD/preeclampsia + SPB & PPROM + latency period",
  "Model 4: maternal age + number of pregnancies",
  "Model 5: previous pregnancy + maternal ancestry",
  "Model 6: GA + weight + sex"
)
names(titles) <- names(files)

genomewide_p <- 5e-8
suggestive_p <- 5e-5

genomewide_line <- -log10(genomewide_p)
suggestive_line <- -log10(suggestive_p)

n_label <- 20

read_plink_logistic <- function(file) {
  df <- tryCatch(
    fread(
      file,
      header = TRUE,
      fill = TRUE,
      data.table = FALSE,
      na.strings = c("NA", "NaN", ".", "nan")
    ),
    error = function(e) {
      stop(paste("Could not read file:", file, "|", e$message))
    }
  )
  
  if (nrow(df) == 0) {
    stop(paste("File is empty:", file))
  }
  
  df
}

prepare_manhattan_data <- function(df) {
  required_cols <- c("SNP", "CHR", "BP", "P", "TEST")
  
  if (!all(required_cols %in% names(df))) {
    stop(
      paste(
        "Missing required columns:",
        paste(setdiff(required_cols, names(df)), collapse = ", ")
      )
    )
  }
  
  df <- df %>%
    filter(TEST == "ADD") %>%
    mutate(
      SNP = as.character(SNP),
      CHR = as.numeric(as.character(CHR)),
      BP  = as.numeric(as.character(BP)),
      P   = as.numeric(as.character(P))
    ) %>%
    filter(!is.na(SNP), SNP != "") %>%
    filter(!is.na(CHR), !is.na(BP), !is.na(P)) %>%
    filter(CHR %in% 1:22) %>%
    filter(P > 0, P <= 1) %>%
    arrange(CHR, BP)
  
  if (nrow(df) == 0) {
    stop("No valid rows after filtering TEST == 'ADD'")
  }
  
  chr_sizes <- df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(tot = cumsum(chr_len) - chr_len)
  
  df <- df %>%
    left_join(chr_sizes, by = "CHR") %>%
    mutate(
      BPcum = BP + tot,
      logP = -log10(P)
    )
  
  axis_df <- df %>%
    group_by(CHR) %>%
    summarise(center = (min(BPcum) + max(BPcum)) / 2, .groups = "drop")
  
  list(data = df, axis = axis_df)
}

prepare_qq_data <- function(df) {
  required_cols <- c("P", "TEST")
  
  if (!all(required_cols %in% names(df))) {
    stop(
      paste(
        "Missing required columns:",
        paste(setdiff(required_cols, names(df)), collapse = ", ")
      )
    )
  }
  
  qq_df <- df %>%
    filter(TEST == "ADD") %>%
    mutate(P = as.numeric(as.character(P))) %>%
    filter(!is.na(P), P > 0, P <= 1) %>%
    arrange(P)
  
  if (nrow(qq_df) == 0) {
    stop("No valid rows after filtering TEST == 'ADD' for QQ plot")
  }
  
  n <- nrow(qq_df)
  
  qq_df <- qq_df %>%
    mutate(
      observed = -log10(P),
      expected = -log10(ppoints(n))
    )
  
  lambda <- median(qchisq(1 - qq_df$P, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
  
  list(data = qq_df, lambda = lambda)
}

make_manhattan_plot <- function(file, model_name, n_label = 8) {
  if (!file.exists(file)) {
    message("File not found: ", file)
    return(NULL)
  }
  
  df <- tryCatch(
    read_plink_logistic(file),
    error = function(e) {
      message("Read error in ", model_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(df)) return(NULL)
  
  prepared <- tryCatch(
    prepare_manhattan_data(df),
    error = function(e) {
      message("Column/data error in ", model_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(prepared)) return(NULL)
  
  plot_df <- prepared$data
  axis_df <- prepared$axis
  
  if (nrow(plot_df) == 0) {
    message("No valid rows for ", model_name)
    return(NULL)
  }
  
  label_df <- plot_df %>%
    filter(logP > suggestive_line) %>%
    arrange(P) %>%
    slice_head(n = n_label) %>%
    mutate(label_y = pmax(logP + 0.30, suggestive_line + 0.45))
  
  plot_df <- plot_df %>%
    mutate(
      chr_group = factor(CHR %% 2),
      top_signal = ifelse(SNP %in% label_df$SNP, "yes", "no")
    )
  
  p <- ggplot(plot_df, aes(x = BPcum, y = logP)) +
    geom_point(
      color = ifelse(plot_df$CHR %% 2 == 0, "grey35", "grey70"),
      size = 0.80,
      alpha = 0.85
    ) +
    geom_point(
      data = subset(plot_df, top_signal == "yes"),
      color = "red3",
      size = 1.5,
      alpha = 0.95
    ) +
    geom_label_repel(
      data = label_df,
      aes(label = SNP),
      size = 3.0,
      min.segment.length = 0,
      box.padding = 0.30,
      point.padding = 0.20,
      label.padding = 0.15,
      max.overlaps = Inf,
      segment.size = 0.30,
      fill = "white",
      color = "black"
    ) +
    geom_hline(
      aes(yintercept = genomewide_line, linetype = "Genome-wide"),
      color = "red",
      linewidth = 0.6
    ) +
    geom_hline(
      aes(yintercept = suggestive_line, linetype = "Suggestive"),
      color = "blue",
      linewidth = 0.6
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "Genome-wide" = "solid",
        "Suggestive"  = "dashed"
      ),
      labels = c(
        "Genome-wide" = expression("Genome-wide significance (" * italic(p) < 5 %*% 10^-8 * ")"),
        "Suggestive"  = expression("Suggestive significance (" * italic(p) < 5 %*% 10^-5 * ")")
      )
    ) +
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df$CHR,
      expand = expansion(mult = c(0.005, 0.01))
    ) +
    scale_y_continuous(
      limits = c(0, 9),
      breaks = seq(0, 9, by = 2),
      expand = c(0, 0)
    ) +
    labs(
      title = model_name,
      x = "Chromosome",
      y = expression(-log[10](italic(P))),
      linetype = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      axis.title.x = element_text(
        size = 13,
        face = "bold",
        margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 13,
        face = "bold",
        margin = margin(r = 8)
      ),
      axis.text.x = element_text(
        size = 10,
        margin = margin(t = 3)
      ),
      axis.text.y = element_text(size = 10),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.4),
      axis.ticks.length = unit(0.16, "cm"),
      legend.position = "inside",
      legend.position.inside = c(0.97, 0.995),
      legend.justification.inside = c(1, 1),
      legend.direction = "vertical",
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.6, "cm"),
      legend.margin = margin(t = 4, r = 6, b = 4, l = 6),
      legend.background = element_rect(
        fill  = "white",
        color = "black",
        linewidth = 0.4
      ),
      plot.margin = margin(14, 14, 8, 8)
    ) +
    guides(
      linetype = guide_legend(override.aes = list(color = c("red", "blue")))
    )
  
  return(p)
}

make_qq_plot <- function(file, model_name) {
  if (!file.exists(file)) {
    message("File not found: ", file)
    return(NULL)
  }
  
  df <- tryCatch(
    read_plink_logistic(file),
    error = function(e) {
      message("Read error in ", model_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(df)) return(NULL)
  
  prepared <- tryCatch(
    prepare_qq_data(df),
    error = function(e) {
      message("QQ data error in ", model_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(prepared)) return(NULL)
  
  qq_df <- prepared$data
  lambda <- prepared$lambda
  
  max_axis <- ceiling(max(c(qq_df$expected, qq_df$observed), na.rm = TRUE) * 10) / 10
  
  p <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(
      color = "grey35",
      size = 0.8,
      alpha = 0.8
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "red",
      linewidth = 0.6
    ) +
    annotate(
      "text",
      x = max_axis * 0.12,
      y = max_axis * 0.92,
      label = paste0("lambda = ", sprintf("%.3f", lambda)),
      hjust = 0,
      size = 4
    ) +
    scale_x_continuous(
      limits = c(0, max_axis),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, max_axis),
      expand = c(0, 0)
    ) +
    labs(
      title = model_name,
      x = expression(Expected~~-log[10](italic(P))),
      y = expression(Observed~~-log[10](italic(P)))
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      axis.title.x = element_text(
        size = 13,
        face = "bold",
        margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 13,
        face = "bold",
        margin = margin(r = 8)
      ),
      axis.text = element_text(size = 10),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.4),
      axis.ticks.length = unit(0.16, "cm"),
      plot.margin = margin(14, 14, 8, 8)
    )
  
  return(p)
}


manhattan_plots <- list()
qq_plots <- list()

for (nm in names(files)) {
  message("Processing: ", nm)
  
  p_manhattan <- make_manhattan_plot(files[[nm]], titles[[nm]], n_label = n_label)
  if (!is.null(p_manhattan)) {
    manhattan_plots[[nm]] <- p_manhattan
  }
  
  p_qq <- make_qq_plot(files[[nm]], titles[[nm]])
  if (!is.null(p_qq)) {
    qq_plots[[nm]] <- p_qq
  }
}

if (length(manhattan_plots) > 0) {
  combined_manhattan <- wrap_plots(manhattan_plots, ncol = 2) +
    plot_annotation(title = "Manhattan plots across all logistic regression models") &
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      )
    )
  
  ggsave(
    filename = "ALL_models_manhattan.png",
    plot = combined_manhattan,
    width = 16,
    height = 18,
    dpi = 300,
    bg = "white"
  )
  
  message("Saved: ALL_models_manhattan.png")
} else {
  message("No Manhattan plots were generated.")
}

if (length(qq_plots) > 0) {
  combined_qq <- wrap_plots(qq_plots, ncol = 2) +
    plot_annotation(title = "QQ plots across all logistic regression models") &
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      )
    )
  
  ggsave(
    filename = "ALL_models_qq.png",
    plot = combined_qq,
    width = 14,
    height = 16,
    dpi = 300,
    bg = "white"
  )
  
  message("Saved: ALL_models_qq.png")
} else {
  message("No QQ plots were generated.")
}
