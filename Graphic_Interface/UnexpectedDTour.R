#------------------------------------------------------------------------------
# Load Required Packages
#------------------------------------------------------------------------------
if (!require(shiny)) install.packages("shiny")
library(shiny) 

if (!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

if (!require(plotly)) install.packages("plotly")
library(plotly)

if (!require(readr)) install.packages("readr")
library(readr)

if (!require(dplyr)) install.packages("dplyr")
library(dplyr)

if (!require(tidyr)) install.packages("tidyr")
library(tidyr)

if (!require(patchwork)) install.packages("patchwork")
library(patchwork)

if (!require(shinyWidgets)) install.packages("shinyWidgets")
library(shinyWidgets)

## -----------------------------------------------------------------------------
## Helper Functions - General Purpose
## -----------------------------------------------------------------------------

#' Extract a Block of Lines from a Text File
#'
#' This function finds a specific figure marker in a vector of text lines,
#' then locates a header pattern after that marker, and extracts all lines
#' between the header and the next figure marker or empty line (stop pattern).
#'
#' @param lines A character vector where each element is a line from the input file.
#' @param figure_marker A string pattern to identify the start of the relevant figure section. Case-insensitive.
#' @param header_pattern A string pattern representing the header row of the data table within the section.
#' @param stop_pattern A regex pattern defining where the data block ends (e.g., start of the next figure or an empty line). Defaults to "^(Figure|\\s*$)".
#' @return A character vector containing the lines of the extracted data block, including the header.
#' @examples
#' \dontrun{
#'   file_lines <- c("Figure 1", "Some text", "Header A B", "1 2", "3 4", "", "Figure 2")
#'   extract_figure_block(file_lines, "Figure 1", "Header A B")
#' }
extract_figure_block <- function(lines, figure_marker, header_pattern, stop_pattern = "^(Figure|\\s*$)") {
  # Find the line number(s) containing the figure marker (case-insensitive).
  fig_idx <- grep(figure_marker, lines, ignore.case = TRUE)
  # Stop if the figure marker isn't found.
  if (length(fig_idx) == 0) {
    stop(paste("No marker found for", figure_marker))
  }
  
  # Find line number(s) matching the header pattern (allowing leading whitespace).
  header_idx <- grep(paste0("^\\s*", header_pattern), lines)
  # Filter header indices to only include those *after* the first figure marker found.
  header_idx <- header_idx[header_idx > fig_idx[1]]
  # Stop if no header line is found after the marker.
  if (length(header_idx) == 0) {
    stop("No header line found after the marker.")
  }
  # Take the first header index found after the figure marker.
  header_idx <- header_idx[1]
  
  # Find line number(s) matching the stop pattern (next figure or empty line).
  stop_idx <- grep(stop_pattern, lines)
  # Filter stop indices to only include those *after* the header line.
  stop_idx <- stop_idx[stop_idx > header_idx]
  
  # Determine the end line index for the block.
  # If no stop pattern is found after the header, take all lines until the end of the file.
  # Otherwise, take lines up to (but not including) the line with the stop pattern.
  end_idx <- if (length(stop_idx) == 0) length(lines) else stop_idx[1] - 1
  
  # Return the subset of lines representing the data block (header included).
  lines[header_idx:end_idx]
}


## -----------------------------------------------------------------------------
## Helper Functions - Specific to Source1 (Figures 7, 8, 9, 10, S2)
## -----------------------------------------------------------------------------

#' Parse Data Block with Potential Formatting Issues
#'
#' Reads a block of text lines into a data frame using read.table.
#' It removes empty lines first and attempts to handle potential extra columns
#' by truncating if `expected_cols` is provided. Issues a warning if the
#' number of columns doesn't match the expectation.
#'
#' @param block A character vector representing the lines of the data block.
#' @param expected_cols Optional: An integer specifying the number of columns expected in the data frame.
#' @return A data frame parsed from the input block.
parseDataWithFix <- function(block) {
  # Clean up extra white spaces
  block <- gsub("\\s+", " ", block)
  block <- trimws(block)
  
  # Remove any lines that are completely empty or contain only whitespace.
  cleaned_block <- block[!grepl("^\\s*$", block)]
  
  # Read the cleaned block into a data frame. `fill = TRUE` adds NA for missing values in rows.
  # `stringsAsFactors = FALSE` prevents character columns from being converted to factors.
  df <- read.table(text = cleaned_block, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, sep = " ")
  
  # Return the parsed data frame.
  df
}

#' Parse Figure 10 Data Block
#'
#' A specialized parser for Figure 10 data, using `read.table` without a header
#' and allowing skipping initial lines if needed.
#'
#' @param lines A character vector representing the lines of the data block for Figure 10.
#' @param skipN Number of lines to skip at the beginning of the block before reading data.
#' @return A data frame parsed from the input block.
parseFigure10Block <- function(lines, skipN = 0) {
  # Read the data block using read.table.
  # `header = FALSE` assumes no header row in the data itself.
  # `skip = skipN` skips the specified number of initial lines.
  # `fill = TRUE` handles rows with potentially fewer columns.
  # `stringsAsFactors = FALSE` keeps character data as characters.
  df <- read.table(text = lines, header = FALSE, skip = skipN,
                   fill = TRUE, stringsAsFactors = FALSE)
  # Return the parsed data frame.
  df
}


## -----------------------------------------------------------------------------
## Helper Functions - Specific to Source2 (Figures 4, 5, 6, 11, S3)
## -----------------------------------------------------------------------------

#' Parse Data for Figure 4
#'
#' Extracts the data block for Figure 4 using `extract_figure_block`,
#' parses it into a data frame, assigns expected column names, converts
#' specific columns to numeric, and removes rows with missing values in key columns.
#'
#' @param all_lines A character vector containing all lines from the uploaded file.
#' @return A data frame containing the processed data for Figure 4.
parseDataFig4 <- function(all_lines) {
  # Extract the relevant lines for Figure 4 using the general helper function.
  # Marker: "Figure 4", Header starts with: "C_dist"
  block_lines <- extract_figure_block(all_lines, "Figure 4", "C_dist")
  # Collapse the lines back into a single string with newline separators for read.table.
  table_text <- paste(block_lines, collapse = "\n")
  # Read the text into a data frame. Assumes no header within the extracted block itself (header=FALSE).
  df <- read.table(text = table_text, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, skip = 1) # Add skip = 1
  # Define the expected column names.
  expected_names <- c("C_dist", "N_P1", "N_P2", "D_N", "ABBA", "BABA", "D")
  # Check if the parsed data frame has at least the expected number of columns.
  if(ncol(df) < length(expected_names)){
    stop("Extracted Figure 4 table has fewer columns than expected.")
  }
  # Assign the expected names to the first N columns.
  colnames(df)[1:length(expected_names)] <- expected_names
  # Convert specific columns to numeric type.
  df$C_dist <- as.numeric(df$C_dist)
  df$D_N <- as.numeric(df$D_N)
  df$D <- as.numeric(df$D)
  # Remove rows where any of the key columns (C_dist, D_N, D) have NA values.
  df <- df[complete.cases(df[, c("C_dist", "D_N", "D")]), ] # Select relevant columns for complete.cases check
  # Return the cleaned data frame.
  df
}


#' Parse Data for Figure 5 (Main Plot)
#'
#' Extracts data for the main plot of Figure 5, parses it, assigns column names,
#' converts key columns to numeric, and removes rows with missing values in those columns.
#'
#' @param all_lines A character vector containing all lines from the uploaded file.
#' @return A data frame containing the processed data for the main plot of Figure 5.
parseDataFig5_main <- function(all_lines) {
  # Extract the relevant lines for Figure 5's main plot.
  # Marker: "Figure 5, main plot", Header starts with: "C_dist"
  block_lines <- extract_figure_block(all_lines, "Figure 5, main plot", "C_dist")
  # Collapse lines into a single string.
  table_text <- paste(block_lines, collapse = "\n")
  # Read the text into a data frame, assuming the first line is the header (header=TRUE).
  df_main <- read.table(text = table_text, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  # Assign specific column names. Note: Adjust if the actual header names differ.
  colnames(df_main) <- c("C_dist", "NP1_1", "NP1_2", "NP1_3", "NP1_S", "NP1_4", "NP2", "D_value", "Extra") # Assumes 9 columns
  # Convert key columns to numeric.
  df_main$C_dist  <- as.numeric(df_main$C_dist)
  df_main$D_value <- as.numeric(df_main$D_value)
  # Remove rows with NA in C_dist or D_value.
  df_main <- df_main[complete.cases(df_main$C_dist, df_main$D_value), ]
  # Return the cleaned data frame.
  df_main
}

#' Parse Data for Figure 5 (Inset Plot)
#'
#' Extracts data for the inset plot of Figure 5, parses it, assigns column names,
#' converts key columns to numeric, and removes rows with missing values.
#'
#' @param all_lines A character vector containing all lines from the uploaded file.
#' @return A data frame containing the processed data for the inset plot of Figure 5.
parseDataFig5_inset <- function(all_lines) {
  # Extract the relevant lines for Figure 5's inset plot.
  # Marker: "Figure 5, inset", Header starts with: "C_dist"
  block_lines <- extract_figure_block(all_lines, "Figure 5, inset", "C_dist")
  # Collapse lines into a single string.
  table_text <- paste(block_lines, collapse = "\n")
  # Read the text into a data frame (header=TRUE).
  df_inset <- read.table(text = table_text, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  # Assign specific column names.
  colnames(df_inset) <- c("C_distance", "NP1_1", "NP1_2", "NP1_3", "NP1_S", "NP1_4", "NP2", "D_value", "Extra") # Assumes 9 columns
  # Rename 'C_distance' to 'C_dist' for consistency.
  colnames(df_inset)[colnames(df_inset) == "C_distance"] <- "C_dist"
  # Convert key columns to numeric.
  df_inset$C_dist  <- as.numeric(df_inset$C_dist)
  df_inset$D_value <- as.numeric(df_inset$D_value)
  # Remove rows with NA in C_dist or D_value.
  df_inset <- df_inset[complete.cases(df_inset$C_dist, df_inset$D_value), ]
  # Return the cleaned data frame.
  df_inset
}


#' Parse and Process Data for Figure 6
#'
#' Extracts the data block for Figure 6, parses it, assigns meaningful column names,
#' filters out empty rows, and then transforms the data from a wide format
#' (multiple D and SE columns) to a long format suitable for plotting with ggplot2.
#' @param all_lines A character vector containing all lines from the uploaded file.
#' @return A data frame in long format containing the processed data for Figure 6.
parseDataFig6 <- function(all_lines) {
  # Extract the relevant lines for Figure 6.
  # Marker: "Figure 6", Header starts with: "Triplet"
  block_lines <- extract_figure_block(all_lines, "Figure 6", "Triplet")
  # Collapse lines into a single string.
  table_text <- paste(block_lines, collapse = "\n")
  # Define the expected column names for the wide format.
  new_names <- c("Triplet",
                 "D_P1P2_poly", "SE_P1P2_poly",
                 "D_P1P1_poly", "SE_P1P1_poly",
                 "D_P1P1S_poly", "SE_P1P1S_poly",
                 "D_P1P2_mono", "SE_P1P2_mono",
                 "D_P1P1_mono", "SE_P1P1_mono",
                 "D_P1P1S_mono", "SE_P1P1S_mono")
  # Read the data (no header in the block itself).
  df6 <- read.table(text = table_text, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, skip = 1) # Add skip = 1
  # Check if enough columns were parsed.
  if(ncol(df6) < length(new_names)){
    stop("Figure 6 block has fewer columns than expected.")
  }
  # Assign the column names to the first N columns.
  colnames(df6)[1:length(new_names)] <- new_names
  # Filter out any rows where the Triplet is NA or an empty string.
  df6 <- df6 %>% filter(!is.na(Triplet) & Triplet != "")
  
  # Transform data from wide to long format.
  df6_long <- df6 %>%
    # Pivot longer for D values: gather all columns starting with "D_" into two columns:
    # 'ComparisonType' (containing the original column name like "D_P1P2_poly")
    # 'Dvalue' (containing the corresponding D value).
    pivot_longer(
      cols = starts_with("D_"),
      names_to = "ComparisonType",
      values_to = "Dvalue"
    ) %>%
    # Pivot longer for SE values.
    pivot_longer(
      cols = starts_with("SE_"),
      names_to = "ComparisonTypeSE",
      values_to = "SEvalue"
    ) %>%
    # Filter to keep only rows where the D and SE values correspond to the same comparison
    filter(sub("SE_", "", ComparisonTypeSE) == sub("D_", "", ComparisonType)) %>%
    # Clean up the ComparisonType and ComparisonTypeSE columns by removing the "D_" and "SE_" prefixes.
    mutate(
      ComparisonType = sub("^D_", "", ComparisonType),
      ComparisonTypeSE = sub("^SE_", "", ComparisonTypeSE)
    ) %>%
    # Separate the combined ComparisonType (e.g., "P1P2_poly") into 'Scenario' (P1P2) and 'ArchaicsCategory' (poly).
    separate(ComparisonType, into = c("Scenario", "ArchaicsCategory"), sep = "_", extra = "merge", fill = "right") %>%
    # Create a more descriptive 'Comparison' column based on the 'Scenario'.
    mutate(
      Comparison = case_when(
        Scenario == "P1P2" ~ "YorFrc",     # Yoruba vs French
        Scenario == "P1P1" ~ "YorYor",     # Yoruba vs Yoruba
        Scenario == "P1P1S" ~ "YorSpike",  # Yoruba vs Spike-in
        TRUE ~ Scenario                   # Fallback
      ),
      # Determine Mutation Type (Transition/Transversion) based on whether the Triplet name ends with "_TV".
      MutationType = if_else(grepl("_TV$", Triplet), "TV", "TS")
    )
  
  # Convert Dvalue and SEvalue columns to numeric, removing header rows.
  df6_long <- df6_long %>%
    mutate(
      Dvalue = as.numeric(Dvalue),
      SEvalue = as.numeric(SEvalue)
    ) %>%
    filter(!is.na(Dvalue)) # Filter out rows where Dvalue is NA (header rows)
  
  # Return the processed long-format data frame.
  df6_long
}

#' Parse Data for Figure 11
#'
#' Extracts two separate data tables associated with Figure 11 (one starting with "slope",
#' one with "Het_dif") using a nested helper function and returns them as a list.
#'
#' @param all_lines A character vector containing all lines from the uploaded file.
#' @return A list containing two data frames: `slope` and `het`.
parseDataFig11 <- function(all_lines) {
  # Define a nested helper function specific to extracting tables within Figure 11 section.
  extract_table_block_11 <- function(lines, start_pattern = "^slope\\s+") {
    # Find the starting line based on the provided pattern (e.g., line starting with "slope ").
    start_idx <- grep(start_pattern, lines)
    if (!length(start_idx)) {
      stop(paste("No start pattern found:", start_pattern))
    }
    # Find potential stopping lines (next Figure marker or empty lines) *after* the start index.
    stop_candidates <- grep("^(Figure|\\s*$)", lines)
    stop_candidates <- stop_candidates[stop_candidates > start_idx[1]]
    if (!length(stop_candidates)) {
      # Handle case where no explicit stop line is found (might need adjustment based on file format).
      # Defaulting to end of file might be too broad. Consider adding a check or alternative stop pattern.
      warning("No explicit stopping line (Next Figure or empty line) found after start_pattern. Reading until end of block/file. Review data.")
      end_idx <- length(lines) # Or potentially stop earlier based on other heuristics if needed
    } else {
      end_idx <- stop_candidates[1] - 1
    }
    
    # Extract the lines belonging to this table.
    table_txt <- lines[start_idx[1]:end_idx]
    df <- read.table(text = table_txt, header = TRUE, fill = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
    df
  }
  
  # Extract the "slope" table using the nested helper.
  df_fig11_slope <- extract_table_block_11(all_lines, "^slope\\s+")
  # Extract the "heterozygosity difference" table using the nested helper.
  df_fig11_het    <- extract_table_block_11(all_lines, "^Het_dif\\s+")
  
  # Function to clean and convert a data frame
  clean_and_convert <- function(df) {
    for (col in names(df)) {
      # Remove any non-numeric characters (except ".")
      df[[col]] <- gsub("[^0-9.-]", "", df[[col]])
      # Replace empty strings with NA
      df[[col]][df[[col]] == ""] <- NA
      # Convert to numeric
      df[[col]] <- as.numeric(df[[col]])
    }
    return(df)
  }
  
  df_fig11_het <- clean_and_convert(df_fig11_het)
  df_fig11_slope <- clean_and_convert(df_fig11_slope)
  
  # Return both data frames in a named list.
  list(slope = df_fig11_slope, het = df_fig11_het)
}

#' Parse Data for Figure S3
#'
#' Extracts the data block for Figure S3, identifies the header, reads the data,
#' assigns expected column names, and returns the data frame. 
#' @param all_lines A character vector containing all lines from the uploaded file.
#' @return A data frame containing the parsed data for Figure S3.
parseDataFigS3 <- function(all_lines) {
  # Find the line containing the "Figure S3" marker.
  figS3_idx <- grep("Figure S3", all_lines, ignore.case = TRUE)
  if(length(figS3_idx) == 0){
    stop("No 'Figure S3' marker found in file")
  }
  # Find the header line (starting with "C_dist") *after* the Figure S3 marker.
  header_idx <- grep("^C_dist", all_lines)
  header_idx <- header_idx[header_idx > figS3_idx[1]][1] # Get the first one after the marker
  if(is.na(header_idx)){
    stop("No header line starting with 'C_dist' found after Figure S3")
  }
  # Extract all lines from the header onwards.
  block_lines <- all_lines[header_idx:length(all_lines)]
  # Filter these lines to keep only the header and lines starting with a number (or '-'),
  table_lines <- block_lines[grepl("^(C_dist)|(^-?[0-9])", block_lines)]
  # Collapse the filtered lines into a single string.
  table_text <- paste(table_lines, collapse="\n")
  # Read the data text into a data frame (header=TRUE).
  figS3_df <- read.table(text = table_text, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  # Define the expected column names.
  expected_names <- c("C_dist",
                      "N_P1_C_A_NA", "N_P2_C_A_NA",
                      "N_P1_C_A_NB", "N_P2_C_A_NB",
                      "N_P1_C_B_NA", "N_P2_C_B_NA",
                      "Dprime_C_A_NA", "Dprime_C_A_NB", "Dprime_C_B_NA")
  # Assign names if the parsed data frame has enough columns.
  if(ncol(figS3_df) >= length(expected_names)){
    colnames(figS3_df)[1:length(expected_names)] <- expected_names
  } else {
    warning("Figure S3 parsed data has fewer columns than expected. Column names may be incorrect.")
  }
  # Return the parsed data frame.
  figS3_df
}

#' Transform Figure S3 Data to Long Format
#'
#' Parses the Figure S3 data using `parseDataFigS3`, selects relevant D' columns,
#' pivots them into a long format, and recodes the 'Scenario' column for clarity.
#'
#' @param all_lines A character vector containing all lines from the uploaded file.
#' @return A data frame in long format containing the processed data for Figure S3 plotting.
figS3_long <- function(all_lines) {
  # Parse the raw Figure S3 data.
  df <- parseDataFigS3(all_lines)
  # Select the distance column and the D' columns, then pivot to long format.
  df_long <- df %>%
    select(C_dist, Dprime_C_A_NA, Dprime_C_A_NB, Dprime_C_B_NA) %>%
    pivot_longer(cols = -C_dist, names_to = "Scenario", values_to = "BranchAsym") # BranchAsym is the D' value
  
  # Recode the scenario names derived from column names into more descriptive labels.
  # The original names encode conditions like C(himpanzee)=A, N(eanderthal)=A etc.
  df_long$Scenario <- recode(df_long$Scenario,
                             Dprime_C_A_NA = "C=A, N=A", # Chimpanzee=A, Neanderthal=A
                             Dprime_C_A_NB = "C=A, N=B", # Chimpanzee=A, Neanderthal=B
                             Dprime_C_B_NA = "C=B, N=A") # Chimpanzee=B, Neanderthal=A
  # Return the long-format data frame.
  df_long
}


## -----------------------------------------------------------------------------
## UI (User Interface) Definition
## -----------------------------------------------------------------------------
# Defines the layout and appearance of the Shiny application in the web browser.
ui <- fluidPage(
  # Add custom CSS styles to the head of the HTML page.
  tags$head(
    tags$style(HTML("
      /* Card styling for sections (file input, buttons) */
      .section-card {
        background-color: white; /* White background */
        border-radius: 5px;      /* Rounded corners */
        padding: 15px;           /* Inner spacing */
        margin-bottom: 30px;     /* Vertical space between cards */
        box-shadow: 0 2px 5px rgba(0,0,0,0.1); /* Subtle shadow */
      }

      /* General button styling */
      .btn {
        font-size: 14px;         /* Text size */
        padding: 10px 15px;      /* Top/bottom and left/right padding */
        width: 100%;             /* Make buttons fill container width */
        margin-bottom: 12px;     /* Space below each button */
        border-radius: 4px;      /* Slightly rounded corners */
        transition: all 0.2s ease; /* Smooth transition for hover effects */
      }

      /* Button hover effect */
      .btn:hover {
        background-color: #3498db; /* Change background on hover */
        color: white;            /* Change text color on hover */
      }

      /* Section heading styling*/
      h4.section-heading {
        border-bottom: 2px solid #3498db; /* Blue bottom border */
        padding-bottom: 8px;      /* Space between text and border */
        margin-bottom: 20px;      /* Space below heading */
        color: #2c3e50;           /* Dark text color */
      }

      /* Styling for the file input container */
      .file-input-container {
        margin-bottom: 25px; /* Space below file input */
      }

      /* General styling for form controls (like file input text) */
      .form-control {
        height: 38px;       /* Control height */
        font-size: 14px;    /* Text size */
      }

      /* Styling for the plot output container */
      .plot-container {
        background-color: white; /* White background */
        border-radius: 5px;      /* Rounded corners */
        padding: 15px;           /* Inner spacing */
        box-shadow: 0 2px 5px rgba(0,0,0,0.1); /* Subtle shadow */
        height: 650px;           /* Fixed height for the plot area */
      }
    "))
  ),
  
  # Application title displayed at the top.
  div(class = "title-panel",
      h2("Unexpected D-tour: Combined Shiny App", style = "margin-top: 20px;")
  ),
  
  # Main row layout for the content.
  fluidRow(
    style = "margin-top: 30px;", # Add space below the title.
    # Left sidebar column (width 3 out of 12).
    column(3,
           # File upload section card.
           div(class = "section-card file-input-container",
               fileInput("datafile", "Upload Text File", # Input ID and label
                         accept = c(".txt"),             # Accept only .txt files
                         buttonLabel = "Browse...",      # Label for the browse button
                         placeholder = "No file selected") # Text shown before file selection
           ),
           
           # Main figures button section card.
           div(class = "section-card",
               h4("Main Figures", class = "section-heading"), # Section title
               # Action buttons for each main figure. Clicking triggers server logic.
               actionButton("fig4Button", "Plot Figure 4", class = "btn btn-default"),
               actionButton("fig5Button", "Plot Figure 5", class = "btn btn-default"),
               actionButton("fig6Button", "Plot Figure 6", class = "btn btn-default"),
               actionButton("fig7Button", "Plot Figure 7", class = "btn btn-default"),
               actionButton("fig8Button", "Plot Figure 8", class = "btn btn-default"),
               actionButton("fig9Button", "Plot Figure 9", class = "btn btn-default"),
               actionButton("fig10Button", "Plot Figure 10", class = "btn btn-default"),
               actionButton("fig11Button", "Plot Figure 11", class = "btn btn-default")
           ),
           
           # Supplementary figures button section card.
           div(class = "section-card",
               h4("Supplementary Figures", class = "section-heading"), # Section title
               # Action buttons for supplementary figures.
               actionButton("figS2Button", "Plot Figure S2", class = "btn btn-default"),
               actionButton("figS3Button", "Plot Figure S3", class = "btn btn-default")
           )
    ),
    
    # Right side main panel column (width 9 out of 12) for displaying the plot.
    column(9,
           div(class = "plot-container",
               # Placeholder for the plotly output, linked to 'plotOutput' in the server.
               plotlyOutput("plotOutput", height = "620px") # Height adjusted slightly less than container
           )
    )
  ),
  
  # Absolute panel for the dropdown instructions button (fixed position).
  absolutePanel(
    bottom = 30, right = 30, fixed = TRUE, # Position: 30px from bottom, 30px from right
    # Use dropdownButton from shinyWidgets for a pop-up style menu.
    dropdownButton(
      label = "Instructions",             # Button label (tooltip)
      icon = icon("info-circle"),       # Button icon (info circle)
      status = "primary",               # Button color style (e.g., blue)
      circle = TRUE,                    # Make the button circular
      width = "300px",                  # Width of the dropdown content area
      up = TRUE,                        # Dropdown opens upwards
      right = TRUE,                     # Align dropdown content to the right edge of the button
      # Content inside the dropdown panel.
      h5("Instructions", style = "margin-top: 0;"), # Title within dropdown
      p("1. Upload the combined text file from outputSource1 and outputSource2."), # Step 1
      p("2. Click any plot button to display the corresponding figure."),       # Step 2
      p("3. Interact with the plot using:"), # Interaction info
      tags$ul(                              # Unordered list for interaction details
        tags$li("Hover to see data values"),
        tags$li("Zoom in/out by clicking on +/- icons"),
        tags$li("Pan by clicking and dragging"),
        tags$li("Double-click to reset the view"),
        tags$li("Click the camera icon to save as image")
      ),
      p(tags$b("Note:"), "Figure 9 requires running the code twice for each population pair.") # Specific note for Fig 9
    )
  )
)


## -----------------------------------------------------------------------------
## Server Logic
## -----------------------------------------------------------------------------
# Defines the server-side logic of the Shiny application.
# This function takes 'input' and 'output' arguments provided by Shiny.
server <- function(input, output, session) {
  
  # Create a reactive value to store the ID of the figure currently requested for plotting.
  # Initialized to NULL. `reactiveVal` creates a reactive source.
  currentFigure <- reactiveVal(NULL)
  
  #----------------------------------------------------------------------------
  # Observers for Action Buttons
  #----------------------------------------------------------------------------
  # These 'observeEvent' blocks listen for clicks on the action buttons defined in the UI.
  # When a button is clicked, they update the 'currentFigure' reactive value
  # with the corresponding figure ID (e.g., "fig4", "figS3").
  
  observeEvent(input$fig4Button,  { currentFigure("fig4") })
  observeEvent(input$fig5Button,  { currentFigure("fig5") })
  observeEvent(input$fig6Button,  { currentFigure("fig6") })
  observeEvent(input$fig7Button,  { currentFigure("fig7") })
  observeEvent(input$fig8Button,  { currentFigure("fig8") })
  observeEvent(input$fig9Button,  { currentFigure("fig9") })
  observeEvent(input$fig10Button, { currentFigure("fig10") })
  observeEvent(input$fig11Button, { currentFigure("fig11") })
  observeEvent(input$figS2Button, { currentFigure("figS2") })
  observeEvent(input$figS3Button, { currentFigure("figS3") })
  
  #----------------------------------------------------------------------------
  # Reactive Expression: Read File Lines
  #----------------------------------------------------------------------------
  # This reactive expression reads the lines from the uploaded file.
  # It depends on `input$datafile`. If `input$datafile` changes (a new file is uploaded),
  # this expression automatically re-executes.
  fileLines <- reactive({
    # `req(input$datafile)` ensures that this code only runs *after* a file has been uploaded.
    # It prevents errors before the first file upload.
    req(input$datafile)
    # Read all lines from the temporary path where Shiny stores the uploaded file.
    readLines(input$datafile$datapath)
  })
  
  #----------------------------------------------------------------------------
  # Reactive Expression: Generate Plot Data and ggplot Object
  #----------------------------------------------------------------------------
  # This is the main reactive expression for generating the plots.
  # It depends on `currentFigure()` and `fileLines()`. It re-executes whenever
  # the user clicks a new figure button (changing `currentFigure`) or uploads
  # a new file (changing `fileLines`).
  plotFigures <- reactive({
    # Require that both `currentFigure` and `fileLines` are available.
    req(currentFigure(), fileLines())
    # Get the current content of the uploaded file.
    all_lines <- fileLines()
    # Get the ID of the figure to plot (e.g., "fig4").
    figChoice <- currentFigure()
    
    # Define a standard title format
    standard_title_layout <- function(title_text) {
      list(
        text = title_text,
        y = 0.95,       # Adjust vertical position as needed
        x = 0.5,       # Center horizontally
        font = list(size = 16), # Adjust font size as needed
        xanchor = "center",
        yanchor = "top"
      )
    }
    
    # Define standard margins
    standard_margins <- list(l = 50, r = 50, b = 50, t = 70, pad = 4)
    # Conditional logic to parse data and create the appropriate plot based on `figChoice`.
    # Each 'if' or 'else if' block corresponds to a specific figure button.
    
    # --- Figure 4 ---
    if (figChoice == "fig4") {
      df <- parseDataFig4(all_lines) # Parse data using the specific helper
      # Create the first panel (D_N vs C_dist)
      p4a <- ggplot(df, aes(x = C_dist, y = D_N)) +
        geom_point(color = "black", size = 2) +
        labs(x = "D'_C (C_dist)", y = "D'_N") +
        theme_minimal()
      # Create the second panel (D vs C_dist)
      p4b <- ggplot(df, aes(x = C_dist, y = D)) +
        geom_point(color = "black", size = 2) +
        labs(x = "D'_C (C_dist)", y = "D") +
        theme_minimal()
      # Combine the two plots side-by-side using plotly's subplot function.
      # Share X and Y axes for easier comparison.
      combined <- subplot(ggplotly(p4a), ggplotly(p4b), nrows = 1, shareX = TRUE, shareY = TRUE) %>%
        layout(
          title = standard_title_layout("Figure 4"),
          margin = standard_margins
        )
      return(combined)
    }
    
    # --- Figure 5 ---
    else if (figChoice == "fig5") {
      df_main <- parseDataFig5_main(all_lines)     # Parse main plot data
      df_inset <- parseDataFig5_inset(all_lines)   # Parse inset plot data
      # Create the main plot
      pMain <- ggplot(df_main, aes(x = C_dist, y = D_value)) +
        geom_point(color = "blue") +
        labs(x = "D'(C-Y1,C-X)", y = "D'(N-Y1,N-X)") +
        theme_minimal()
      # Create the inset plot
      pInset <- ggplot(df_inset, aes(x = C_dist, y = D_value)) +
        geom_point(color = "red") +
        labs(x = "D'(C-Y1,C-X)", y = "D'(N-Y1,N-X)") + # Consistent labels
        theme_minimal()
      # Combine plots side-by-side using subplot.
      combined <- subplot(ggplotly(pMain), ggplotly(pInset), nrows = 1) %>%
        layout(
          title = standard_title_layout("Figure 5"),
          margin = standard_margins
        )
      return(combined)
      
      # --- Figure 6 ---
    } else if (figChoice == "fig6") {
      df6_long <- parseDataFig6(all_lines) # Parse and process data into long format
      # Filter data for the two panels (polymorphic vs monomorphic sites)
      df_poly <- df6_long %>% filter(ArchaicsCategory == "poly")
      df_mono <- df6_long %>% filter(ArchaicsCategory == "mono")
      
      # Define a helper function to create one panel of the Figure 6 plot
      plot_figure6_panel <- function(data) {
        # Define custom colors for different comparison/mutation types
        custom_colors <- c(
          "YorFrc_TV"   = "red",       # Yoruba-French, Transversion
          "YorFrc_TS"   = "gold",      # Yoruba-French, Transition
          "YorYor_TV"   = "black",     # Yoruba-Yoruba, Transversion
          "YorYor_TS"   = "gray60",    # Yoruba-Yoruba, Transition
          "YorSpike_TV" = "darkblue",  # Yoruba-Spike, Transversion
          "YorSpike_TS" = "lightblue"  # Yoruba-Spike, Transition
        )
        # Create the ggplot object
        ggplot(data, aes(x = Triplet, y = Dvalue,
                         # Map color and group aesthetics to the interaction of Comparison and MutationType
                         color = interaction(Comparison, MutationType, sep = "_"),
                         group = interaction(Comparison, MutationType, sep = "_"))) +
          geom_point() +        
          geom_line() +         
          geom_errorbar(aes(ymin = Dvalue - SEvalue, ymax = Dvalue + SEvalue), # Add error bars
                        width = 0.3, alpha = 0.5) + # Set error bar width and transparency
          scale_color_manual(values = custom_colors) + # Apply the custom color scheme
          labs(x = "Triplet", y = "D") +             
          theme_bw() +                               
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
      }
      
      # Create the two panels using the helper function
      p1 <- plot_figure6_panel(df_poly)
      p2 <- plot_figure6_panel(df_mono)
      
      # Combine the panels using subplot
      combined <- subplot(ggplotly(p1), ggplotly(p2), nrows = 1, shareX = FALSE, shareY = FALSE) %>%
        layout(
          title = standard_title_layout("Figure 6"),
          annotations = list(
            list(x = 0.25, y = 1.02, text = "Panel (a): Polymorphic Archaics",
                 showarrow = FALSE, xref = "paper", yref = "paper", xanchor = 'center', font = list(size = 14)),
            list(x = 0.75, y = 1.02, text = "Panel (b): Monomorphic Archaics",
                 showarrow = FALSE, xref = "paper", yref = "paper", xanchor = 'center', font = list(size = 14))
          ),
          legend = list(orientation = "h", x = 0.1, y = -0.2),
          margin = standard_margins
        )
      return(combined)
      
      # --- Figure 11 ---
    } else if (figChoice == "fig11") {
      fig11_data <- parseDataFig11(all_lines) # Parse the two tables into a list
      df_het <- fig11_data$het     # Extract heterozygosity difference data
      df_slope <- fig11_data$slope # Extract slope data
      
      # Clean dataframes right before plotting
      clean_dataframe <- function(df) {
        df %>%
          mutate_all(~ifelse(is.infinite(.) | is.na(.), NA, .))
      }
      
      df_het <- clean_dataframe(df_het)
      df_slope <- clean_dataframe(df_slope)
      
      # Create the first panel (D vs Heterozygosity Difference)
      p_het <- ggplot(df_het, aes(x = Het_dif)) +
        geom_point(aes(y = D_ALL,   color = "D_ALL"),
                   data = . %>% filter(is.finite(Het_dif) & is.finite(D_ALL))) +
        geom_smooth(aes(y = D_ALL,  color = "D_ALL"), method = "lm", se = FALSE,
                    data = . %>% filter(is.finite(Het_dif) & is.finite(D_ALL))) +
        geom_point(aes(y = D_MajorA, color = "D_MajorA"),
                   data = . %>% filter(is.finite(Het_dif) & is.finite(D_MajorA))) +
        geom_point(aes(y = D_MajorB, color = "D_MajorB"),
                   data = . %>% filter(is.finite(Het_dif) & is.finite(D_MajorB))) +
        labs(
          title = "D vs Heterozygosity Difference", # Panel title
          x = "Heterozygosity Difference (Bantu - Japanese)",
          y = "D(...)"
        ) +
        scale_color_manual( # Define colors for legend/mapping
          name = "D Type",  # Legend title
          values = c("D_ALL" = "black", "D_MajorA" = "blue", "D_MajorB" = "red")
        ) +
        theme_minimal()
      
      # Create the second panel (D vs Slope)
      p_slope <- ggplot(df_slope, aes(x = slope)) +
        geom_point(aes(y = D_ALL,   color = "D_ALL"),
                   data = . %>% filter(is.finite(slope) & is.finite(D_ALL))) +
        geom_smooth(aes(y = D_ALL,  color = "D_ALL"), method = "lm", se = FALSE,
                    data = . %>% filter(is.finite(slope) & is.finite(D_ALL))) +
        geom_point(aes(y = D_MajorA, color = "D_MajorA"),
                   data = . %>% filter(is.finite(slope) & is.finite(D_MajorA))) +
        geom_point(aes(y = D_MajorB, color = "D_MajorB"),
                   data = . %>% filter(is.finite(slope) & is.finite(D_MajorB))) +
        labs(
          title = "D vs Slope (Het vs Dist)", # Panel title
          x = "Slope (Het vs. Dist. from Africa)",
          y = "D(...)" # Y-axis label shared implicitly by subplot
        ) +
        scale_color_manual( # Define colors (consistent with p_het)
          name = "D Type", # Legend title
          values = c("D_ALL" = "black", "D_MajorA" = "blue", "D_MajorB" = "red")
        ) +
        theme_minimal()

      combined <- subplot(ggplotly(p_het), ggplotly(p_slope), nrows = 1, shareY = TRUE, titleX = TRUE) %>%
        layout(
          title = standard_title_layout("Figure 11"),
          margin = standard_margins
        )
      return(combined)
      
      # --- Figure S3 ---
    } else if (figChoice == "figS3") {
      # Parse and transform data to long format using the specific helper.
      # Also renames columns and recodes 'Scenario' for better plot labels.
      df_long <- figS3_long(all_lines) %>%
        mutate(Scenario = recode(Scenario, 
                                 "C=A, N=A" = "Altai = chimpanzee = A",
                                 "C=A, N=B" = "Altai = A, chimpanzee = B",
                                 "C=B, N=A" = "Altai = B, chimpanzee = A")) %>%
        rename( # Rename columns for plot axes/legend labels
          HetDiff = C_dist,
          BranchAsymmetry = BranchAsym,
          SiteClass = Scenario
        )
      
      # Create the ggplot object
      p <- ggplot(df_long, aes(x = HetDiff, y = BranchAsymmetry, color = SiteClass)) +
        geom_point() +
        ggtitle("Figure S3") +
        labs(x = "Het Diff (Yoruba - French)", y = "Branch Asymmetry (D')") +
        theme_minimal() +
        scale_color_manual( # Define colors for the different site classes
          name = "Site Class", # Legend title
          values = c(
            "Altai = chimpanzee = A"    = "blue",
            "Altai = A, chimpanzee = B" = "black",
            "Altai = B, chimpanzee = A" = "red"
          ))

      return(ggplotly(p) %>%
               layout(
                 title = standard_title_layout("Figure S3"),
                 margin = standard_margins
               ))
      
      # --- Figure 7 ---
    } else if (figChoice == "fig7") {
      # Extract the relevant block for Figure 7.
      block <- extract_figure_block(all_lines, "Figure 7", "Triplet")
      # Parse the data using the general parser.
      df <- parseDataWithFix(block) 
      colnames(df)[1] <- "Triplet" # Ensure first column is named 'Triplet'
      # Pivot data to long format for plotting multiple lines.
      df_long <- pivot_longer(df, cols = -Triplet, names_to = "Type", values_to = "Value")
      # Recode the 'Type' column (original column names) for better legend labels.
      df_long$Type <- recode(df_long$Type,
                             "TS"   = "AABAA_TS",
                             "TS.1" = "PAAAA_TS",
                             "TS.2" = "PBBBA_TS",
                             "TV"   = "AABAA_TV",
                             "TV.1" = "PAAAA_TV",
                             "TV.2" = "PBBBA_TV"
      )
      # Define a color map for the different types.
      color_map <- c(
        "AABAA_TS" = "red",
        "PAAAA_TS" = "purple",
        "PBBBA_TS" = "darkblue",
        "AABAA_TV" = "yellow",
        "PAAAA_TV" = "pink",
        "PBBBA_TV" = "lightblue"
      )
      
      # Create the ggplot object.
      p <- ggplot(df_long, aes(x = Triplet, y = Value, color = Type, group = Type)) +
        geom_line(linewidth = 1) + 
        geom_point() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_color_manual(values = color_map, name = "Category")
      
      return(ggplotly(p) %>%
               layout(
                 title = standard_title_layout("Figure 7"),
                 margin = standard_margins
               ))
      
      # --- Figure 8 ---
    } else if (figChoice == "fig8") {
      # Extract the block for Figure 8.
      block <- extract_figure_block(all_lines, "Figure 8", "Triplet")
      # Parse data, expecting 7 columns.
      df <- parseDataWithFix(block)
      colnames(df)[1] <- "Triplet"
      colnames(df)[2:7] <- c("propn_ABBA", "propn_BABA", "diff_ABBA",
                             "propn_PBBBA", "propn_PAAAA", "diff_PBBBA")
      # Create scatter plot comparing the two difference metrics.
      p <- ggplot(df, aes(x = diff_ABBA, y = diff_PBBBA)) +
        geom_point(size = 3) + # Add points
        geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) + # Add linear regression line
        xlab("Relative Contribution to ABBAs (propn_ABBA - propn_BABA)") + 
        ylab("Relative Mutation Rate (propn_PBBBA - propn_PAAAA)") + 
        theme_minimal()
      
      return(ggplotly(p) %>%
               layout(
                 title = standard_title_layout("Figure 8"),
                 margin = standard_margins
               ))
      
      # --- Figure 9 ---
    } else if (figChoice == "fig9") {
      # Extract the block for Figure 9.
      block <- extract_figure_block(all_lines, "Figure 9", "Triplet")
      df <- parseDataWithFix(block)
      colnames(df)[1] <- "Triplet"
      colnames(df)[2:5] <- c("TV_D_major_A", "TV_D_major_B", "TS_D_major_A", "TS_D_major_B")
      df$Triplet <- as.character(df$Triplet) # Ensure Triplet is character for potential factor conversion later
      
      # Separate data for Transversions (TV) and Transitions (TS) and pivot longer.
      df_tv <- df %>%
        select(Triplet, TV_D_major_A, TV_D_major_B) %>%
        pivot_longer(cols = c("TV_D_major_A", "TV_D_major_B"), names_to = "LineID", values_to = "D")
      df_ts <- df %>%
        select(Triplet, TS_D_major_A, TS_D_major_B) %>%
        pivot_longer(cols = c("TS_D_major_A", "TS_D_major_B"), names_to = "LineID", values_to = "D")
      
      # Add columns indicating Major allele (A or B) and Mutation type (TV or TS).
      df_tv <- df_tv %>% mutate(Major = ifelse(grepl("A$", LineID), "A", "B"), # Check last character of LineID
                                Mutation = "TV")
      df_ts <- df_ts %>% mutate(Major = ifelse(grepl("A$", LineID), "A", "B"),
                                Mutation = "TS")
      
      # Combine the TV and TS data frames.
      df_combined <- bind_rows(df_tv, df_ts) %>%
        mutate(LineID2 = paste(Major, Mutation, sep = "_")) # Create a combined ID for grouping/coloring
      
      # Define a specific order for the Triplets on the x-axis.
      triplets_order <- c("AAA", "AAC", "AAG", "AAT", "CAA", "CAC", "CAG", "CAT",
                          "GAA", "GAC", "GAG", "GAT", "TAA", "TAC", "TAG", "TAT",
                          "ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT",
                          "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT")
      # Convert Triplet column to a factor with the specified order.
      df_combined$Triplet <- factor(df_combined$Triplet, levels = triplets_order, ordered = TRUE)
      
      # Define colors for the four combinations of Major allele and Mutation type.
      color_map <- c("B_TV" = "purple", "B_TS" = "darkblue",
                     "A_TV" = "red",    "A_TS" = "darkgreen")
      
      # Create the ggplot object.
      p <- ggplot(df_combined, aes(x = Triplet, y = D, group = LineID2, color = LineID2)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        scale_color_manual(values = color_map, name = "Category (Major_Mutation)") + 
        theme_bw(base_size = 13) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
              legend.title = element_blank()) + 
        labs(title = "Figure 9", x = "Mutating Triplet", y = "D value")
      
      return(ggplotly(p) %>%
               layout(
                 title = standard_title_layout("Figure 9"),
                 margin = standard_margins
               ))
      
      # --- Figure 10 ---
    } else if (figChoice == "fig10") {
      # Extract the block for Figure 10.
      block <- extract_figure_block(all_lines, "Figure 10", "Triplet") # Assumes header is "Triplet..."
      # Figure out how many lines to skip before the data starts. Often 1 if header line included in block.
      skipN <- 1 # Adjust if the actual header is different or missing in the block
      df_10_noheader <- parseFigure10Block(block, skipN = skipN)
      
      # Attempt to convert columns 2 through 14 (assumed numeric data) to numeric.
      # Using suppressWarnings to avoid messages if conversion fails for some columns (e.g., if text remains).
      suppressWarnings({
        df_10_noheader[, 2:min(14, ncol(df_10_noheader))] <- lapply(df_10_noheader[, 2:min(14, ncol(df_10_noheader))], as.numeric)
      })
      
      # Check if enough columns were parsed (at least 13 expected).
      stopifnot(ncol(df_10_noheader) >= 13) # Stop if data structure seems wrong
      colnames(df_10_noheader)[1] <- "Triplet"
      colnames(df_10_noheader)[2:13] <- c("AP_TV", "AP_TS", "PA_TV", "PA_TS", "BP_TV", "BP_TS",
                                          "PB_TV", "PB_TS", "PP_TV", "PP_TS", "ALL_TV", "ALL_TS")
      if (ncol(df_10_noheader) >= 14) {
        colnames(df_10_noheader)[14] <- "TOTAL"
      }
      # Ensure Triplet is character and filter out any potential empty/NA rows.
      df_10_noheader$Triplet <- as.character(df_10_noheader$Triplet)
      df_10 <- df_10_noheader %>% filter(!is.na(Triplet) & Triplet != "")
      
      # Pivot the data to a long format suitable for faceting.
      df_long <- df_10 %>%
        pivot_longer(
          # Select all the D* value columns.
          cols = c(AP_TV, AP_TS, PA_TV, PA_TS, BP_TV, BP_TS,
                   PB_TV, PB_TS, PP_TV, PP_TS, ALL_TV, ALL_TS),
          # Specify how to split the column names:
          # names_to: Names of the new columns to create ('Scenario', 'MutationType').
          # names_pattern: A regex to capture parts of the old column names.
          #   "(.*)" captures the Scenario part (e.g., "AP", "PA", "ALL").
          #   "_(TV|TS)$" captures the MutationType part ("TV" or "TS") at the end.
          names_to = c("Scenario", "MutationType"),
          names_pattern = "(.*)_(TV|TS)$",
          # values_to: Name of the new column for the cell values ('Dstar').
          values_to = "Dstar"
        ) %>%
        # Recode the captured 'Scenario' abbreviations into more readable labels.
        mutate(Scenario = recode(Scenario,
                                 "AP" = "A/P", "PA" = "P/A", "BP" = "B/P",
                                 "PB" = "P/B", "PP" = "P/P", "ALL" = "ALL"))
      
      # Separate the data for individual scenarios (A/P, P/A, etc.) and the overall 'ALL' scenario.
      df_scenario <- df_long %>% filter(Scenario != "ALL")
      # Extract the 'ALL' scenario data to be used as a reference line in plots.
      df_overall  <- df_long %>% filter(Scenario == "ALL") %>%
        select(Triplet, MutationType, OverallDstar = Dstar) # Select relevant columns and rename Dstar
      
      # Join the scenario-specific data with the overall D* data.
      df_final <- df_scenario %>% left_join(df_overall, by = c("Triplet", "MutationType"))
      
      df_final$Scenario <- factor(df_final$Scenario,
                                  levels = c("A/P", "P/A", "B/P", "P/B", "P/P"),
                                  ordered = TRUE)
      df_final$MutationType <- factor(df_final$MutationType,
                                      levels = c("TV", "TS"), ordered = TRUE)
      # Convert Triplet to a factor using the order from the original data frame to maintain consistency.
      tripletLevels <- unique(df_10$Triplet)
      df_final$Triplet <- factor(df_final$Triplet, levels = tripletLevels, ordered = TRUE)
      # Use `as.numeric(Triplet)` for x-axis to treat it as continuous for lines, but use labels later.
      p <- ggplot(df_final, aes(x = as.numeric(Triplet))) +
        # Add the overall D* line (black, thinner). Group=1 ensures a single line per facet.
        geom_line(aes(y = OverallDstar, group = 1), color = "black", linewidth = 0.5) +
        # Add the scenario-specific D* line (red, thicker).
        geom_line(aes(y = Dstar, group = 1), color = "red", linewidth = 0.8) +
        # Add a horizontal line at y=0 for reference.
        geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
        # Create a grid of plots faceted by Scenario (rows) and MutationType (columns).
        # `scales = "free_y"` allows y-axis scales to differ between rows.
        # `switch = "y"` moves the y-axis facet labels (Scenario names) to the left.
        facet_grid(Scenario ~ MutationType, scales = "free_y", switch = "y") +
        # Set the x-axis breaks and labels to correspond to the Triplet levels.
        scale_x_continuous(breaks = seq_along(levels(df_final$Triplet)),
                           labels = levels(df_final$Triplet)) +
        labs(x = "Triplet", y = "D*", title = "Figure 10") + 
        theme_minimal(base_size = 12) + # Minimal theme
        theme(axis.text.x = element_text(angle = 90, hjust = 1), 
              strip.placement = "outside") # Ensure facet labels are outside axes

      return(ggplotly(p) %>%
               layout(
                 title = standard_title_layout("Figure 10"),
                 margin = standard_margins
               ))
      
      # --- Figure S2 ---
    } else if (figChoice == "figS2") {
      # Extract the block for Figure S2.
      block <- extract_figure_block(all_lines, "Figure S2", "Triplet")
      # Parse data using the general parser.
      df <- parseDataWithFix(block)
      colnames(df)[1] <- "Triplet" # Ensure first column is Triplet
      # Pivot to long format.
      df_long <- df %>% pivot_longer(cols = -Triplet, names_to = "Type", values_to = "Value")
      
      # Create the ggplot object.
      p <- ggplot(df_long, aes(x = Triplet, y = Value, color = Type, group = Type)) +
        geom_line(linewidth = 1) +
        geom_point() +
        ggtitle("Figure S2") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_color_manual( 
          name = "Category", 
          values = c( "AABAA" = "red", "AAABA" = "blue", "PAAAA" = "black" ) 
        )

      return(ggplotly(p) %>%
               layout(
                 title = standard_title_layout("Figure S2"),
                 margin = standard_margins
               ))
    }
    
    # Add a fallback or error message if figChoice doesn't match any known figure
    else {
      # Return NULL or a placeholder plot/message if the figure choice is not handled
      # For example:
      # return(plotly_empty(type = "scatter", mode = "text") %>% layout(title = paste("Plotting logic for", figChoice, "not yet implemented.")))
      # Or simply stop with an error:
      stop(paste("Plotting logic for figure choice '", figChoice, "' is not defined.", sep=""))
    }
    
    
  }) 
  
  #----------------------------------------------------------------------------
  # Render Plot Output
  #----------------------------------------------------------------------------
  # This defines the output element 'plotOutput' that was created in the UI.
  # It calls the `plotFigures()` reactive expression to get the plot object.
  output$plotOutput <- renderPlotly({
    # Require necessary inputs before attempting to render.
    req(input$datafile, currentFigure())
    # Execute the reactive expression to generate and return the plot.
    plotFigures()
  }) # End of renderPlotly
  
} # End of server function

#------------------------------------------------------------------------------
# Run the Shiny Application
#------------------------------------------------------------------------------
# This command starts the Shiny app, combining the UI and Server definitions.
shinyApp(ui, server)