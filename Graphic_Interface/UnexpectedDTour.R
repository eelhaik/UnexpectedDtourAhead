# Install/load required packages
if (!require(shiny)) install.packages("shiny")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(plotly)) install.packages("plotly")
if (!require(readr)) install.packages("readr")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(patchwork)) install.packages("patchwork")
if (!require(shinyWidgets)) install.packages("shinyWidgets")
library(shinyWidgets)
library(shiny)
library(ggplot2)
library(plotly)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)

## ------------------
## Helper functions
## ------------------

# Common helper to extract a block of lines from the file given a figure marker and header pattern.
extract_figure_block <- function(lines, figure_marker, header_pattern, stop_pattern = "^(Figure|\\s*$)") {
  fig_idx <- grep(figure_marker, lines, ignore.case = TRUE)
  if (length(fig_idx) == 0) {
    stop(paste("No marker found for", figure_marker))
  }
  header_idx <- grep(paste0("^\\s*", header_pattern), lines)
  header_idx <- header_idx[header_idx > fig_idx[1]]
  if (length(header_idx) == 0) {
    stop("No header line found after the marker.")
  }
  header_idx <- header_idx[1]
  stop_idx <- grep(stop_pattern, lines)
  stop_idx <- stop_idx[stop_idx > header_idx]
  end_idx <- if (length(stop_idx) == 0) length(lines) else stop_idx[1] - 1
  lines[header_idx:end_idx]
}

## ----- Source1 (Figures 7, 8, 9, 10, S2) helper functions -----

parseDataWithFix <- function(block, expected_cols = NULL) {
  cleaned_block <- block[!grepl("^\\s*$", block)]
  df <- read.table(text = cleaned_block, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  if (!is.null(expected_cols)) {
    if(ncol(df) > expected_cols) {
      df <- df[, 1:expected_cols, drop = FALSE]
    }
    if(ncol(df) != expected_cols) {
      warning(sprintf("Parsed %d columns, but expected %d.", ncol(df), expected_cols))
    }
  }
  df
}

parseFigure10Block <- function(lines, skipN = 0) {
  df <- read.table(text = lines, header = FALSE, skip = skipN,
                   fill = TRUE, stringsAsFactors = FALSE)
  df
}

## ----- Source2 (Figures 4, 5, 6, 11, S3) helper functions -----

parseDataFig4 <- function(all_lines) {
  block_lines <- extract_figure_block(all_lines, "Figure 4", "C_dist")
  table_text <- paste(block_lines, collapse = "\n")
  df <- read.table(text = table_text, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  expected_names <- c("C_dist", "N_P1", "N_P2", "D_N", "ABBA", "BABA", "D")
  if(ncol(df) < length(expected_names)){
    stop("Extracted Figure 4 table has fewer columns than expected.")
  }
  colnames(df)[1:length(expected_names)] <- expected_names
  df$C_dist <- as.numeric(df$C_dist)
  df$D_N <- as.numeric(df$D_N)
  df$D <- as.numeric(df$D)
  df <- df[complete.cases(df), ]
  df
}

parseDataFig5_main <- function(all_lines) {
  block_lines <- extract_figure_block(all_lines, "Figure 5, main plot", "C_dist")
  table_text <- paste(block_lines, collapse = "\n")
  df_main <- read.table(text = table_text, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  colnames(df_main) <- c("C_dist", "NP1_1", "NP1_2", "NP1_3", "NP1_S", "NP1_4", "NP2", "D_value", "Extra")
  df_main$C_dist  <- as.numeric(df_main$C_dist)
  df_main$D_value <- as.numeric(df_main$D_value)
  df_main <- df_main[complete.cases(df_main$C_dist, df_main$D_value), ]
  df_main
}

parseDataFig5_inset <- function(all_lines) {
  block_lines <- extract_figure_block(all_lines, "Figure 5, inset", "C_dist")
  table_text <- paste(block_lines, collapse = "\n")
  df_inset <- read.table(text = table_text, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  colnames(df_inset) <- c("C_distance", "NP1_1", "NP1_2", "NP1_3", "NP1_S", "NP1_4", "NP2", "D_value", "Extra")
  df_inset$C_dist  <- as.numeric(df_inset$C_distance)
  df_inset$D_value <- as.numeric(df_inset$D_value)
  df_inset <- df_inset[complete.cases(df_inset$C_dist, df_inset$D_value), ]
  df_inset
}

parseDataFig6 <- function(all_lines) {
  block_lines <- extract_figure_block(all_lines, "Figure 6", "Triplet")
  table_text <- paste(block_lines, collapse = "\n")
  new_names <- c("Triplet",
                 "D_P1P2_poly", "SE_P1P2_poly",
                 "D_P1P1_poly", "SE_P1P1_poly",
                 "D_P1P1S_poly", "SE_P1P1S_poly",
                 "D_P1P2_mono", "SE_P1P2_mono",
                 "D_P1P1_mono", "SE_P1P1_mono",
                 "D_P1P1S_mono", "SE_P1P1S_mono")
  df6 <- read.table(text = table_text, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  if(ncol(df6) < length(new_names)){
    stop("Figure 6 block has fewer columns than expected.")
  }
  colnames(df6)[1:length(new_names)] <- new_names
  df6 <- df6 %>% filter(!is.na(Triplet) & Triplet != "")
  df6_long <- df6 %>%
    pivot_longer(
      cols = starts_with("D_"),
      names_to = "ComparisonType",
      values_to = "Dvalue"
    ) %>%
    pivot_longer(
      cols = starts_with("SE_"),
      names_to = "ComparisonTypeSE",
      values_to = "SEvalue"
    ) %>%
    filter(sub("SE_", "", ComparisonTypeSE) == sub("D_", "", ComparisonType)) %>%
    mutate(
      ComparisonType = sub("^D_", "", ComparisonType),
      ComparisonTypeSE = sub("^SE_", "", ComparisonTypeSE)
    ) %>%
    separate(ComparisonType, into = c("Scenario", "ArchaicsCategory"), sep = "_", extra = "merge", fill = "right") %>%
    mutate(
      Comparison = case_when(
        Scenario == "P1P2" ~ "YorFrc",
        Scenario == "P1P1" ~ "YorYor",
        Scenario == "P1P1S" ~ "YorSpike",
        TRUE ~ Scenario
      ),
      MutationType = if_else(grepl("_TV$", Triplet), "TV", "TS")
    )
  df6_long <- df6_long %>% mutate(
    Dvalue = as.numeric(Dvalue),
    SEvalue = as.numeric(SEvalue)
  )
  df6_long
}

parseDataFig11 <- function(all_lines) {
  extract_table_block_11 <- function(lines, start_pattern = "^slope\\s+") {
    start_idx <- grep(start_pattern, lines)
    if (!length(start_idx)) {
      stop(paste("No start pattern found:", start_pattern))
    }
    stop_candidates <- grep("^(Figure|\\s*$)", lines)
    stop_candidates <- stop_candidates[stop_candidates > start_idx[1]]
    if (!length(stop_candidates)) {
      stop("No stopping line found after start_pattern. Adjust if needed.")
    }
    end_idx <- stop_candidates[1] - 1
    table_txt <- lines[start_idx[1]:end_idx]
    df <- read.table(text = table_txt, header = TRUE, fill = TRUE, strip.white = TRUE)
    df
  }
  df_fig11_slope <- extract_table_block_11(all_lines, "^slope\\s+")
  df_fig11_het   <- extract_table_block_11(all_lines, "^Het_dif\\s+")
  list(slope = df_fig11_slope, het = df_fig11_het)
}

parseDataFigS3 <- function(all_lines) {
  figS3_idx <- grep("Figure S3", all_lines, ignore.case = TRUE)
  if(length(figS3_idx) == 0){
    stop("No 'Figure S3' marker found in file")
  }
  header_idx <- grep("^C_dist", all_lines)
  header_idx <- header_idx[header_idx > figS3_idx[1]][1]
  if(is.na(header_idx)){
    stop("No header line starting with 'C_dist' found after Figure S3")
  }
  block_lines <- all_lines[header_idx:length(all_lines)]
  table_lines <- block_lines[grepl("^(C_dist)|(^-?[0-9])", block_lines)]
  table_text <- paste(table_lines, collapse="\n")
  figS3_df <- read.table(text = table_text, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  expected_names <- c("C_dist",
                      "N_P1_C_A_NA", "N_P2_C_A_NA",
                      "N_P1_C_A_NB", "N_P2_C_A_NB",
                      "N_P1_C_B_NA", "N_P2_C_B_NA",
                      "Dprime_C_A_NA", "Dprime_C_A_NB", "Dprime_C_B_NA")
  if(ncol(figS3_df) >= length(expected_names)){
    colnames(figS3_df)[1:length(expected_names)] <- expected_names
  }
  figS3_df
}

figS3_long <- function(all_lines) {
  df <- parseDataFigS3(all_lines)
  df_long <- df %>%
    select(C_dist, Dprime_C_A_NA, Dprime_C_A_NB, Dprime_C_B_NA) %>%
    pivot_longer(cols = -C_dist, names_to = "Scenario", values_to = "BranchAsym")
  df_long$Scenario <- recode(df_long$Scenario,
                             Dprime_C_A_NA = "C=A, N=A",
                             Dprime_C_A_NB = "C=A, N=B",
                             Dprime_C_B_NA = "C=B, N=A")
  df_long
}

## ---------------------
## UI (Chronological order of figures)
## ---------------------

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* Card styling for sections */
      .section-card {
        background-color: white;
        border-radius: 5px;
        padding: 15px;
        margin-bottom: 30px;  /* spacing between sections */
        box-shadow: 0 2px 5px rgba(0,0,0,0.1);
      }
      
      /* Button styling */
      .btn {
        font-size: 14px;
        padding: 10px 15px;
        width: 100%;
        margin-bottom: 12px;
        border-radius: 4px;
        transition: all 0.2s ease;
      }
      
      .btn:hover {
        background-color: #3498db;
        color: white;
      }
      
      /* Section headings */
      h4.section-heading {
        border-bottom: 2px solid #3498db;
        padding-bottom: 8px;
        margin-bottom: 20px;
        color: #2c3e50;
      }
      
      /* File input styling */
      .file-input-container {
        margin-bottom: 25px;
      }
      
      .form-control {
        height: 38px;
        font-size: 14px;
      }
      
      /* Plot container */
      .plot-container {
        background-color: white;
        border-radius: 5px;
        padding: 15px;
        box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        height: 650px;
      }
    "))
  ),
  
  # App title
  div(class = "title-panel",
      h2("Unexpected D-tour: Combined Shiny App", style = "margin-top: 20px;")
  ),
  
  fluidRow(
    style = "margin-top: 30px;",
    # Left sidebar with controls
    column(3,
           # File upload section
           div(class = "section-card file-input-container",
               fileInput("datafile", "Upload Text File", 
                         accept = c(".txt"),
                         buttonLabel = "Browse...",
                         placeholder = "No file selected")
           ),
           
           # Main figures section
           div(class = "section-card",
               h4("Main Figures", class = "section-heading"),
               actionButton("fig4Button", "Plot Figure 4", class = "btn btn-default"),
               actionButton("fig5Button", "Plot Figure 5", class = "btn btn-default"),
               actionButton("fig6Button", "Plot Figure 6", class = "btn btn-default"),
               actionButton("fig7Button", "Plot Figure 7", class = "btn btn-default"),
               actionButton("fig8Button", "Plot Figure 8", class = "btn btn-default"),
               actionButton("fig9Button", "Plot Figure 9", class = "btn btn-default"),
               actionButton("fig10Button", "Plot Figure 10", class = "btn btn-default"),
               actionButton("fig11Button", "Plot Figure 11", class = "btn btn-default")
           ),
           
           # Supplementary figures section
           div(class = "section-card",
               h4("Supplementary Figures", class = "section-heading"),
               actionButton("figS2Button", "Plot Figure S2", class = "btn btn-default"),
               actionButton("figS3Button", "Plot Figure S3", class = "btn btn-default")
           )
    ),
    
    # Right side plot area
    column(9,
           div(class = "plot-container",
               plotlyOutput("plotOutput", height = "620px")
           )
    )
  ),
  
  # Dropdown instructions panel in the bottom right corner
  absolutePanel(
    bottom = 30, right = 30, fixed = TRUE,
    dropdownButton(
      label = "Instructions",
      icon = icon("info-circle"),
      status = "primary",
      circle = TRUE,
      width = "300px",
      up = TRUE,
      right = TRUE,
      # Content of the dropdown
      h5("Instructions", style = "margin-top: 0;"),
      p("1. Upload the combined text file from outputSource1 and outputSource2."),
      p("2. Click any plot button to display the corresponding figure."),
      p("3. Interact with the plot using:"),
      tags$ul(
        tags$li("Hover to see data values"),
        tags$li("Zoom in/out with mouse wheel"),
        tags$li("Pan by clicking and dragging"),
        tags$li("Double-click to reset the view"),
        tags$li("Click the camera icon to save as image")
      ),
      p(tags$b("Note:"), "Figure 9 requires running the code twice for each population pair.")
    )
  )
)

## ---------------------
## Server
## ---------------------

server <- function(input, output, session) {
  
  currentFigure <- reactiveVal(NULL)
  
  observeEvent(input$fig4Button, { currentFigure("fig4") })
  observeEvent(input$fig5Button, { currentFigure("fig5") })
  observeEvent(input$fig6Button, { currentFigure("fig6") })
  observeEvent(input$fig7Button, { currentFigure("fig7") })
  observeEvent(input$fig8Button, { currentFigure("fig8") })
  observeEvent(input$fig9Button, { currentFigure("fig9") })
  observeEvent(input$fig10Button, { currentFigure("fig10") })
  observeEvent(input$fig11Button, { currentFigure("fig11") })
  observeEvent(input$figS2Button, { currentFigure("figS2") })
  observeEvent(input$figS3Button, { currentFigure("figS3") })
  
  fileLines <- reactive({
    req(input$datafile)
    readLines(input$datafile$datapath)
  })
  
  plotFigures <- reactive({
    req(currentFigure(), fileLines())
    all_lines <- fileLines()
    figChoice <- currentFigure()
    
    if (figChoice == "fig4") {
      df <- parseDataFig4(all_lines)
      p4a <- ggplot(df, aes(x = C_dist, y = D_N)) +
        ggtitle("Figure 4") +
        geom_point(color = "black", size = 2) +
        labs(x = "D'_C (C_dist)", y = "D'_N") +
        theme_minimal()
      p4b <- ggplot(df, aes(x = C_dist, y = D)) +
        geom_point(color = "black", size = 2) +
        labs(x = "D'_C (C_dist)", y = "D") +
        theme_minimal()
      combined <- subplot(ggplotly(p4a), ggplotly(p4b), nrows = 1, shareX = TRUE, shareY = TRUE)
      return(combined)
      
    } else if (figChoice == "fig5") {
      df_main <- parseDataFig5_main(all_lines)
      df_inset <- parseDataFig5_inset(all_lines)
      pMain <- ggplot(df_main, aes(x = C_dist, y = D_value)) +
        ggtitle("Figure 5") +
        geom_point(color = "blue") +
        labs(x = "D'(C-Y1,C-X)", y = "D'(N-Y1,N-X)") +
        theme_minimal()
      pInset <- ggplot(df_inset, aes(x = C_dist, y = D_value)) +
        geom_point(color = "red") +
        labs(x = "D'(C-Y1,C-X)", y = "D'(N-Y1,N-X)") +
        theme_minimal()
      combined <- subplot(ggplotly(pMain), ggplotly(pInset), nrows = 1)
      return(combined)
      
    } else if (figChoice == "fig6") {
      df6_long <- parseDataFig6(all_lines)
      df_poly <- df6_long %>% filter(ArchaicsCategory == "poly")
      df_mono <- df6_long %>% filter(ArchaicsCategory == "mono")
      
      plot_figure6_panel <- function(data) {
        custom_colors <- c(
          "YorFrc_TV"   = "red",
          "YorFrc_TS"   = "gold",
          "YorYor_TV"   = "black",
          "YorYor_TS"   = "gray60",
          "YorSpike_TV" = "darkblue",
          "YorSpike_TS" = "lightblue"
        )
        ggplot(data, aes(x = Triplet, y = Dvalue,
                         color = interaction(Comparison, MutationType, sep = "_"),
                         group = interaction(Comparison, MutationType, sep = "_"))) +
          ggtitle("Figure 6") +
          geom_point() +
          geom_line() +
          geom_errorbar(aes(ymin = Dvalue - SEvalue, ymax = Dvalue + SEvalue),
                        width = 0.3, alpha = 0.5) +
          scale_color_manual(values = custom_colors) +
          labs(x = "Triplet", y = "D") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      }
      
      p1 <- plot_figure6_panel(df_poly)
      p2 <- plot_figure6_panel(df_mono)
      
      combined <- subplot(ggplotly(p1), ggplotly(p2), nrows = 1, shareX = FALSE, shareY = FALSE) %>%
        layout(annotations = list(
          list(x = 0.25, y = 1, text = "Panel (a)", 
               showarrow = FALSE, xref = "paper", yref = "paper"),
          list(x = 0.75, y = 1, text = "Panel (b)", 
               showarrow = FALSE, xref = "paper", yref = "paper")
        ),
        legend = list(orientation = "h", x = 0, y = -0.15))
      return(combined)
      
    } else if (figChoice == "fig11") {
      fig11_data <- parseDataFig11(all_lines)
      df_het <- fig11_data$het
      df_slope <- fig11_data$slope
      
      p_het <- ggplot(df_het, aes(x = Het_dif)) +
        geom_point(aes(y = D_ALL,   color = "D_ALL")) +
        geom_smooth(aes(y = D_ALL,  color = "D_ALL"), method = "lm", se = FALSE) +
        geom_point(aes(y = D_MajorA, color = "D_MajorA")) +
        geom_point(aes(y = D_MajorB, color = "D_MajorB")) +
        labs(
          x = "Heterozygosity Difference (Bantu - Japanese)",
          y = "D(...)"
        ) +
        scale_color_manual(
          values = c(
            "D_ALL"   = "black",
            "D_MajorA" = "blue",
            "D_MajorB" = "red"
          )
        ) +
        theme_minimal()
      
      p_slope <- ggplot(df_slope, aes(x = slope)) +
        geom_point(aes(y = D_ALL,   color = "D_ALL")) +
        geom_smooth(aes(y = D_ALL,  color = "D_ALL"), method = "lm", se = FALSE) +
        geom_point(aes(y = D_MajorA, color = "D_MajorA")) +
        geom_point(aes(y = D_MajorB, color = "D_MajorB")) +
        ggtitle("Figure 11") +
        labs(
          x = "Slope (Het vs. Dist. from Africa)",
          y = "D(...)"
        ) +
        scale_color_manual(
          values = c(
            "D_ALL"   = "black",
            "D_MajorA" = "blue",
            "D_MajorB" = "red"
          )
        ) +
        theme_minimal()
      
      combined <- subplot(ggplotly(p_het), ggplotly(p_slope), nrows = 1)
      return(combined)
      
    } else if (figChoice == "figS3") {
      df <- parseDataFigS3(all_lines)
      df_long <- figS3_long(all_lines) %>%
        mutate(Scenario = recode(Scenario,
                                 "C=A, N=A" = "Altai = chimpanzee = A",
                                 "C=A, N=B" = "Altai = A, chimpanzee = B",
                                 "C=B, N=A" = "Altai = B, chimpanzee = A"
        )) %>%
        rename(
          HetDiff = C_dist,
          BranchAsymmetry = BranchAsym,
          SiteClass = Scenario
        )
      
      
      p <- ggplot(df_long, aes(x = HetDiff, y = BranchAsymmetry, color = SiteClass)) +
        geom_point() +
        ggtitle("Figure S3") +
        labs(x = "Het Diff (Yoruba - French)", y = "Branch Asymmetry (D')") +
        theme_minimal() +
        scale_color_manual(values = c(
          "Altai = chimpanzee = A" = "blue",
          "Altai = A, chimpanzee = B" = "black",
          "Altai = B, chimpanzee = A" = "red"
        ))
      return(ggplotly(p))
      
    }
    
    else if (figChoice == "fig7") {
      block <- extract_figure_block(all_lines, "Figure 7", "Triplet")
      df <- parseDataWithFix(block)
      colnames(df)[1] <- "Triplet"
      df_long <- pivot_longer(df, cols = -Triplet, names_to = "Type", values_to = "Value")
      df_long$Type <- recode(df_long$Type,
                             "TS"   = "AABAA_TS",
                             "TS.1" = "PAAAA_TS",
                             "TS.2" = "PBBBA_TS",
                             "TV"   = "AABAA_TV",
                             "TV.1" = "PAAAA_TV",
                             "TV.2" = "PBBBA_TV"
      )
      color_map <- c(
        "AABAA_TS" = "red",
        "PAAAA_TS" = "purple",
        "PBBBA_TS" = "darkblue",
        "AABAA_TV" = "yellow",
        "PAAAA_TV" = "pink",
        "PBBBA_TV" = "lightblue"
      )
      
      p <- ggplot(df_long, aes(x = Triplet, y = Value, color = Type, group = Type)) +
        geom_line(size = 1) + 
        geom_point() +
        ggtitle("Figure 7") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_color_manual(values = color_map)
      
      return(ggplotly(p))
      
    } else if (figChoice == "fig8") {
      block <- extract_figure_block(all_lines, "Figure 8", "Triplet")
      df <- parseDataWithFix(block, expected_cols = 7)
      colnames(df)[1] <- "Triplet"
      colnames(df)[2:7] <- c("propn_ABBA", "propn_BABA", "diff_ABBA",
                             "propn_PBBBA", "propn_PAAAA", "diff_PBBBA")
      p <- ggplot(df, aes(x = diff_ABBA, y = diff_PBBBA)) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.8) +
        xlab("Relative Contribution to ABBAs (propn_ABBA - propn_BABA)") +
        ylab("Relative Mutation Rate (propn_PBBBA - propn_PAAAA)") +
        ggtitle("Figure 8") +
        theme_minimal()
      return(ggplotly(p))
      
    } else if (figChoice == "fig9") {
      block <- extract_figure_block(all_lines, "Figure 9", "Triplet")
      df <- parseDataWithFix(block, expected_cols = 5)
      colnames(df)[1] <- "Triplet"
      colnames(df)[2:5] <- c("TV_D_major_A", "TV_D_major_B", "TS_D_major_A", "TS_D_major_B")
      df$Triplet <- as.character(df$Triplet)
      df_tv <- df %>%
        select(Triplet, TV_D_major_A, TV_D_major_B) %>%
        pivot_longer(cols = c("TV_D_major_A", "TV_D_major_B"), names_to = "LineID", values_to = "D")
      df_ts <- df %>%
        select(Triplet, TS_D_major_A, TS_D_major_B) %>%
        pivot_longer(cols = c("TS_D_major_A", "TS_D_major_B"), names_to = "LineID", values_to = "D")
      df_tv <- df_tv %>% mutate(Major = ifelse(grepl("A", LineID), "A", "B"),
                                Mutation = "TV")
      df_ts <- df_ts %>% mutate(Major = ifelse(grepl("A", LineID), "A", "B"),
                                Mutation = "TS")
      df_combined <- bind_rows(df_tv, df_ts) %>%
        mutate(LineID2 = paste(Major, Mutation, sep = "_"))
      triplets_order <- c("AAA", "AAC", "AAG", "AAT",
                          "CAA", "CAC", "CAG", "CAT",
                          "GAA", "GAC", "GAG", "GAT",
                          "TAA", "TAC", "TAG", "TAT",
                          "ACA", "ACC", "ACG", "ACT",
                          "CCA", "CCC", "CCG", "CCT",
                          "GCA", "GCC", "GCG", "GCT",
                          "TCA", "TCC", "TCG", "TCT")
      df_combined$Triplet <- factor(df_combined$Triplet, levels = triplets_order, ordered = TRUE)
      color_map <- c("B_TV" = "purple",
                     "B_TS" = "darkblue",
                     "A_TV" = "red",
                     "A_TS" = "darkgreen")
      p <- ggplot(df_combined, aes(x = Triplet, y = D, group = LineID2, color = LineID2)) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        scale_color_manual(values = color_map) +
        theme_bw(base_size = 13) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              legend.title = element_blank()) +
        labs(title = "Figure 9", x = "Mutating Triplet", y = "D value")
      return(ggplotly(p))
      
    } else if (figChoice == "fig10") {
      block <- extract_figure_block(all_lines, "Figure 10", "Triplet")
      skipN <- 1
      df_10_noheader <- parseFigure10Block(block, skipN = skipN)
      
      df_10_noheader[, 2:14] <- lapply(df_10_noheader[, 2:14], as.numeric)
      
      stopifnot(ncol(df_10_noheader) >= 13)
      colnames(df_10_noheader)[1] <- "Triplet"
      colnames(df_10_noheader)[2:13] <- c("AP_TV", "AP_TS", "PA_TV", "PA_TS", "BP_TV", "BP_TS",
                                          "PB_TV", "PB_TS", "PP_TV", "PP_TS", "ALL_TV", "ALL_TS")
      if (ncol(df_10_noheader) >= 14) {
        colnames(df_10_noheader)[14] <- "TOTAL"
      }
      df_10_noheader$Triplet <- as.character(df_10_noheader$Triplet)
      df_10 <- df_10_noheader %>% filter(!is.na(Triplet) & Triplet != "")
      
      df_long <- df_10 %>%
        pivot_longer(
          cols = c(AP_TV, AP_TS, PA_TV, PA_TS, BP_TV, BP_TS,
                   PB_TV, PB_TS, PP_TV, PP_TS, ALL_TV, ALL_TS),
          names_to = c("Scenario", "MutationType"),
          names_pattern = "(.*)_(TV|TS)$",
          values_to = "Dstar"
        ) %>%
        mutate(Scenario = recode(Scenario,
                                 "AP" = "A/P", "PA" = "P/A", "BP" = "B/P",
                                 "PB" = "P/B", "PP" = "P/P", "ALL" = "ALL"))
      
      df_scenario <- df_long %>% filter(Scenario != "ALL")
      df_overall  <- df_long %>% filter(Scenario == "ALL") %>%
        select(Triplet, MutationType, OverallDstar = Dstar)
      df_final <- df_scenario %>% left_join(df_overall, by = c("Triplet", "MutationType"))
      
      df_final$Scenario <- factor(df_final$Scenario,
                                  levels = c("A/P", "P/A", "B/P", "P/B", "P/P"),
                                  ordered = TRUE)
      df_final$MutationType <- factor(df_final$MutationType,
                                      levels = c("TV", "TS"), ordered = TRUE)
      tripletLevels <- unique(df_10$Triplet)
      df_final$Triplet <- factor(df_final$Triplet, levels = tripletLevels, ordered = TRUE)
      
      p <- ggplot(df_final, aes(x = as.numeric(Triplet))) +
        geom_line(aes(y = OverallDstar, group = 1), color = "black", size = 0.5) +
        geom_line(aes(y = Dstar, group = 1), color = "red", size = 0.8) +
        geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
        facet_grid(Scenario ~ MutationType, scales = "free_y", switch = "y") +
        scale_x_continuous(breaks = seq_along(levels(df_final$Triplet)), 
                           labels = levels(df_final$Triplet)) +
        labs(x = "Triplet",
             y = "D*",
             title = "Figure 10") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      
      return(ggplotly(p))
    }
    else if (figChoice == "figS2") {
      block <- extract_figure_block(all_lines, "Figure S2", "Triplet")
      df <- parseDataWithFix(block)
      colnames(df)[1] <- "Triplet"
      df_long <- df %>% pivot_longer(cols = -Triplet, names_to = "Type", values_to = "Value")
      
      p <- ggplot(df_long, aes(x = Triplet, y = Value, color = Type, group = Type)) +
        geom_line(size = 1) + 
        geom_point() +
        ggtitle("Figure S2") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_color_manual(values = c(
          "AABAA" = "red",
          "AAABA" = "blue",
          "PAAAA" = "black"
        ))
      
      return(ggplotly(p))
    }
    
  })
  
  output$plotOutput <- renderPlotly({
    req(input$datafile, currentFigure())
    plotFigures()
  })
}

shinyApp(ui, server)
