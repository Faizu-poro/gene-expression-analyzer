# app.R

# Load required packages
required_packages <- c("shiny", "plotly", "reactable", "shinyFeedback", "dplyr", "tidyr", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.rstudio.com/")
  }
}

# Load preprocessCore from Bioconductor for quantile normalization
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "http://cran.rstudio.com/")
  }
  BiocManager::install("preprocessCore")
}

library(shiny)
library(plotly)
library(reactable)
library(shinyFeedback)
library(dplyr)
library(tidyr)
library(preprocessCore)
library(ggplot2)

# Increase upload size
options(shiny.maxRequestSize = 30 * 1024^2) 

# Custom GCT file parser
read_gct <- function(file_path) {
  conn <- file(file_path, "r")
  on.exit(close(conn))
  
  # Read version and dimensions
  version <- readLines(conn, n = 1)  # Line 1: Version
  dimensions <- readLines(conn, n = 1)  # Line 2: Dimensions
  dims <- strsplit(dimensions, "\t")[[1]]
  
  # Extract values
  n_data_rows <- as.numeric(dims[1])  # Number of data rows (genes)
  n_cols <- as.numeric(dims[2])  # Number of columns
  n_metadata_rows <- ifelse(length(dims) >= 4, as.numeric(dims[4]), 0)  # Metadata rows
  
  # Read column headers (Line 3)
  col_headers <- strsplit(readLines(conn, n = 1), "\t")[[1]]
  
  # Read metadata rows (if any)
  metadata <- NULL
  if (n_metadata_rows > 0) {
    metadata <- readLines(conn, n = n_metadata_rows)  # Read metadata rows
  }
  
  # Read expression data (only n_data_rows rows)
  expr_data <- read.delim(
    conn,
    header = FALSE,
    nrows = n_data_rows,  # Read only the specified number of data rows
    col.names = col_headers,
    check.names = FALSE
  )
  
  # Convert expression columns to numeric
  expr_matrix <- as.data.frame(lapply(expr_data[, -1], function(x) as.numeric(as.character(x))))
  rownames(expr_matrix) <- expr_data[, 1]  # Use first column as row names (gene IDs)
  
  return(list(
    metadata = metadata,  # Metadata lines
    expression = expr_matrix  # Expression matrix
  ))
}

# Normalization function
normalize_data <- function(data, method) {
  mat <- data$expression
  switch(method,
         "log" = log2(mat + 1),
         "zscore" = scale(mat, center = TRUE, scale = TRUE),
         "quantile" = {
           norm_mat <- preprocessCore::normalize.quantiles(as.matrix(mat))
           colnames(norm_mat) <- colnames(mat)
           rownames(norm_mat) <- rownames(mat)
           norm_mat
         },
         # If "none", return original matrix
         mat)
}

# UI
ui <- fluidPage(
  useShinyFeedback(),
  titlePanel("Gene Expression Analyzer"),
  sidebarLayout(
    sidebarPanel(
      fileInput("gct_file", "Upload GCT File", accept = ".gct"), # define filename and text to appear for file upload
      selectInput("norm_method", "Normalization Method:",
                  choices = c("Log2" = "log",
                              "Z-score" = "zscore",
                              "Quantile" = "quantile")),
      uiOutput("norm_explanation"),  # render normalization explanation
      actionButton("analyze", "Run Analysis", class = "btn-primary"), #define action button with name
      hr(), # horizontal rule to seperate from analysis info
      h4("Analysis Info"), # header for analysis details 
      verbatimTextOutput("data_info") # render from server using gctfile and normalization method selected
    ),
    mainPanel(
      # Define tabs and their titles
      tabsetPanel(
        tabPanel("PCA Plot",
                 plotlyOutput("pca_plot"),
                 uiOutput("download_pca_button")),
        tabPanel("Metadata",
                 reactableOutput("metadata_table")),
        tabPanel("Expression Data",
                 reactableOutput("expression_table"))
      )
    )
  )
)

# Server
server <- function(input, output) {
  
  # normalization explanation
  output$norm_explanation <- renderUI({
    req(input$norm_method) # need norm_method to not be empty
    explanation <- switch(input$norm_method,
                          "log" = paste(
                            "Log2 Transformation:", 
                            "Converts data using a logarithmic scale (base 2).",
                            "Reduces skewness from extreme high/low values.", 
                            "Adds 1 to avoid log(0) errors.",
                            "Useful for visualizing gene expression trends.",
                            sep = "\n"
                          ),
                          "zscore" = paste(
                            "Z-Score Normalization:", 
                            "Centers data: Subtracts the mean (average expression of each gene becomes 0).",
                            "Scales data: Divides by the standard deviation (variation becomes 1).", 
                            "Allows comparison of genes across different scales.",
                            "Use to highlight relative expression differences.",
                            sep = "\n"
                          ),
                          "quantile" = paste(
                            "Quantile Normalization:", 
                            "Forces all samples to have the same statistical distribution.",
                            "Removes technical variations (e.g., batch effects).", 
                            "Preserves biological differences.",
                            "Commonly used in microarray/RNA-seq data.",
                            "Computationally intensive but highly effective.",
                            sep = "\n"
                          )
    )
    helpText(explanation) # return help with respect to norm_method selected by user
  })
  
  # Reactive expression since we need it to change with error on reading the GCT file 
  gct_data <- reactive({ 
    # Make sure file is uploadedd
    req(input$gct_file)
    # try block to catch exceptions
    tryCatch({
      # Check file extension and return error if invalid file
      if (tools::file_ext(input$gct_file$name) != "gct") { # Check if file extension is .gct
        showFeedbackDanger("gct_file", "Please upload a .gct file") # Show error message
        return(NULL)
      }
      
      # Read the GCT file
      data <- read_gct(input$gct_file$datapath) # Get temp path and use custom function to parse gct file and extract expression data
      
      # Return data for try block 
      data
      
    }, error = function(e) {
      # Handle specific exceptions
      if (grepl("invalid GCT format", e$message, ignore.case = TRUE)) {
        showFeedbackDanger("gct_file", "Invalid GCT file format")
      } else {
        showFeedbackDanger("gct_file", "Error reading file")
      }
      return(NULL)
    })
  })
  
  # Run normalization when "Run Analysis" is clicked
  norm_data <- eventReactive(input$analyze, {
    req(gct_data(), input$norm_method) # Gct file and user selected norm_method required
    withProgress(message = "Normalizing data...", {
      normalize_data(gct_data(), input$norm_method) # call normalize function
    })
  })
  
  # Reactive expression to perform PCA on normalized data
  pca_result <- reactive({
    req(norm_data())
    withProgress(message = "Running PCA...", {
      nd <- norm_data()
      # Remove genes (rows) with zero variance
      non_constant <- apply(nd, 1, function(x) {
        sd_val <- sd(x, na.rm = TRUE)
        !is.na(sd_val) && sd_val > 0
      })
      nd_filtered <- nd[non_constant, , drop = FALSE]
      prcomp(t(nd_filtered), scale. = TRUE)
    })
  })
  
  # Create PCA plot object using Plotly
  pca_plot_obj <- reactive({
    pca <- pca_result()
    var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
    plot_ly(as.data.frame(pca$x),
            x = ~PC1, y = ~PC2,
            text = ~rownames(pca$x),
            type = "scatter", mode = "markers",
            hoverinfo = "text+x+y") %>%
      layout(
        xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")),
        yaxis = list(title = paste0("PC2 (", var_exp[2], "%)"))
      )
  })
  
  # Render the interactive PCA plot
  output$pca_plot <- renderPlotly({
    pca_plot_obj()
  })
  
  # Render metadata as a structured table
  output$metadata_table <- renderReactable({
    req(gct_data())
    metadata <- gct_data()$metadata
    if (!is.null(metadata)) {
      # Parse metadata into a structured table
      metadata_list <- strsplit(metadata, "\t")
      metadata_df <- do.call(rbind, lapply(metadata_list, function(row) {
        data.frame(
          Field = row[1],
          Values = paste(row[-1], collapse = ", ")
        )
      }))
      reactable(metadata_df, filterable = TRUE, searchable = TRUE)
    } else {
      reactable(data.frame(Message = "No metadata found."))
    }
  })
  
  # Render expression data as a data table using reactable
  output$expression_table <- renderReactable({
    req(gct_data())
    expr_data <- gct_data()$expression
    reactable(expr_data, filterable = TRUE, searchable = TRUE)
  })
  
  # Display data information
  output$data_info <- renderPrint({
    req(gct_data())
    cat("Samples:", ncol(gct_data()$expression), "\n")
    cat("Genes:", nrow(gct_data()$expression), "\n")
    cat("Normalization:", input$norm_method, "\n")
  })
  
  # download button appears when PCA is created
  output$download_pca_button <- renderUI({
    req(pca_result())  # Only show button if PCA results exist
    downloadButton("dl_pca", "Download PCA Plot")
  })
  
  output$dl_pca <- downloadHandler(
    filename = function() { "pca_plot.png" },
    content = function(file) {
      req(pca_result())
      pca <- pca_result()
      var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
      pc_labels <- colnames(pca$x)[1:2]
      
      # Create static ggplot
      p <- ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) +
        geom_point(color = "blue") +
        labs(
          x = paste0(pc_labels[1], " (", var_exp[1], "%)"),
          y = paste0(pc_labesl[2], " (", var_exp[2], "%)"),
          title = "PCA Plot (PC1 vs PC2)"
        ) +
        theme_minimal()
      
      # Save plot
      ggsave(file, p, width = 10, height = 6, dpi = 300)
    }
  )
}

# Run the Shiny application
shinyApp(ui = ui, server = server)