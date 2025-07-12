library(dplyr)
library(kableExtra)

supp_table <-  function(df, row_col,
                        metrics = NULL,
                        caption = "Supplementary Table S1: Performance metrics for all models.",
                        label = "supp:performance") {
  all_metrics <- unique(sub("_(.*)", "", setdiff(names(df), row_col)))
  if (is.null(metrics)) metrics <- all_metrics
  if (length(metrics) > 3) warning("More than 3 metrics passed — only the first 3 will be shown.")
  metrics <- metrics[1:3]
  
  panel_labels <- c("Panel A", "Panel B", "Panel C")
  
  # Start outer table float
  cat("\\begin{table}[ht]\n\\centering\n")
  
  for (i in seq_along(metrics)) {
    metric <- metrics[i]
    panel <- panel_labels[i]
    
    sub_cols <- grep(paste0("^", metric, "_"), names(df), value = TRUE)
    sub_df <- df[, c(row_col, sub_cols)]
    colnames(sub_df) <- c(row_col, sub(paste0("^", metric, "_"), "", sub_cols))
    
    # Format numeric columns to 3 digits
    for (col in colnames(sub_df)[-1]) {
      sub_df[[col]] <- formatC(as.numeric(sub_df[[col]]), format = "f", digits = 3)
    }
    
    header_row <- setNames(c(1, ncol(sub_df) - 1), c(" ", metric))
    
    # Begin subtable
    cat("\\begin{subtable}[t]{\\linewidth}\n\\centering\n\\vspace{0pt}\n")
    
    # Render table
    kbl <- kable(sub_df, format = "latex", booktabs = TRUE, align = "l", caption = NULL) %>%
      add_header_above(header_row) %>%
      kable_styling(latex_options = c("hold_position"))
    
    print(kbl)
    
    # Add subcaption
    cat(paste0("\\caption{", panel, ": ", metric, " metrics}\n"))
    cat("\\end{subtable}\n\n")
  }
  
  # Final caption and label
  cat(paste0("\\caption{", caption, "}\n"))
  cat(paste0("\\label{", label, "}\n"))
  cat("\\end{table}\n")
}
summary_df <- samplerun_eval %>% 
  group_by(mode) %>% 
  summarise(across(where(is.numeric),
                   list(mean = ~mean(.),
                        median = ~median(.),
                        q25 = ~quantile(., 0.25),
                        q75 = ~quantile(., 0.75)),
                   .names = "{.col}_{.fn}"))
summary_df$mode <- mode_label[summary_df$mode]
summary_df%<>%mutate_if(is.numeric,formatC,format = "f", digits = 3)
# generate_metric_tables(summary_df, row_col = "mode")

sink("Figures/supp_table.tex")
supp_table(summary_df, row_col = "mode")
sink()
