# Gene Expression Analyzer

A Shiny application for analyzing gene expression data from GCT files. Supports data normalization, PCA analysis, and interactive visualization.

## Features
- Upload and parse GCT files.
- Normalization methods: Log2, Z-score, Quantile.
- Interactive PCA plot (PC1 vs. PC2).
- Display metadata and expression data in tabular format.
- Download PCA plots.

## Installation

### Prerequisites
- R (≥ 4.0)
- RStudio (recommended)

### Installation  
The app will auto-install missing packages on first run. For manual setup:  
    ```r
    # Run these commands in R once:  
    install.packages(c("shiny", "plotly", "reactable", "shinyFeedback", "dplyr", "tidyr", "ggplot2))  
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")  
    BiocManager::install("preprocessCore")  

## Steps
1. **Clone the repository**:
    ```bash
    git clone https://github.com/Faizu-poro/gene-expression-analyzer.git
    cd gene-expression-analyzer

2. **Run the app**
    a. Linux:
        ```bash
        R -e "shiny::runApp('app.R')"

    b. Rstudio:
        i. Open Rstudio
        ii. Set the current working directory to the folder containing app.R:
        iii. Click the Run App button or 
            ```r
            shiny::runApp("app.R")
    
    c. R:
        If you are running R outside of Rstudio, open R, set working directory to the cloned repository and execute:
            ```r
            shiny::runApp("app.R")    

### Usage
1. **Upload a GCT File**
    - Click "Upload GCT File" and select your file
    - Supported extenions: .gct

2. **Select Normalization Method**
    - Choose from Log2, Z-Score or Quantile

3. **Run Analysis**
    - Click "Run Analysis" to generate PCA plots and tables

4. **Explore Results**
    - PCA Plot: Interactive plot of PC1 vs PC2
    - Metadata: Structured metadata table
    - Expression Data: Numeric expression matrix

5. Downloading the Results

- Use the **Download PCA Plot** button to save the PCA plot as a PNG.

### Configuration Guide
1. File Upload Limit
    - Purpose: Adjust maximum file size (default: 30MB)
    - Modify in app.R:
        ```r
        options(shiny.maxRequestSize = 30 * 1024^2)  # Change 30 to desired MB

2. Metadata Parsing
    - Purpose: Customize metadata parsing logic
    - Modify read_gct() function in app.R:
        ```r
        # Example: Skip first 57 metadata rows as provided in line 2 cell 4
        # Here dims is elements in line 2 where dimensions of the gct file are normally observed.
        # Fourth element in the sample dataset corresponds to number of rows of metadata information after headers. 
        # Adjust number in gct file as necessary or assign number of metadata lines here to n_metadata_rows
        n_metadata_rows <- ifelse(length(dims) >= 4, as.numeric(dims[4]), 0)

3. Normalization Parameters
    - Purpose: Tweak normalization methods
    - Modify normalize_data() in app.R:
        ```r
        # Example: Change log2 offset. Default offset set to 1 to avoid 0 values.
        "log" = log2(mat + 1)  # Replace 1 with another value

## Application Screenshots

Screenshots of the app (e.g., file upload screen, analysis view, plots) can be found in the `screenshots/` folder.

## Troubleshooting

- **File Format Issues:**  
  Ensure your GCT file conforms to the expected format. The first two lines should contain version/dimension info, followed by metadata and numeric data.
- **Performance:**  
  Larger files may take longer to process. Progress indicators will show the status.
- **Installation Issues:**  
  If packages fail to install, make sure you have write permissions or try running R/RStudio as an administrator.