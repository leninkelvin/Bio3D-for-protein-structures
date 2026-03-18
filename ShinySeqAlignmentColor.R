# Load necessary libraries
library(shiny)
library(Biostrings)
library(msaR)  # For colored MSA visualization
library(DT)
library(shinyjs)
library(shinythemes)
library(seqinr)    # For additional file formats

# Define UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  useShinyjs(),
  tags$head(
    tags$style(HTML("
      .msaR {
        height: 600px !important;
      }
      .well {
        background-color: #f8f9fa;
      }
      .content-wrapper {
        padding: 20px;
      }
      #msaROutput {
        border: 1px solid #ddd;
        border-radius: 5px;
        padding: 10px;
        background-color: white;
      }
    "))
  ),
  titlePanel("MSA Visualization Tool"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("alignmentFile", "Upload Alignment File",
                accept = c(".fasta", ".fa", ".aln", ".clustal", ".txt")),
      
      selectInput("fileFormat", "File Format",
                  choices = c("FASTA" = "fasta", 
                              "Clustal" = "clustal", 
                              "Stockholm" = "stockholm",
                              "MSF" = "msf",
                              "Auto-detect" = "auto"),
                  selected = "auto"),
      
      hr(),
      
      selectInput("colorScheme", "Color Scheme",
                  choices = c("Clustal" = "clustal", 
                              "Taylor" = "taylor", 
                              "Zappo" = "zappo", 
                              "Helix Propensity" = "helix", 
                              "Strand Propensity" = "strand", 
                              "Turn Propensity" = "turn", 
                              "Buried Index" = "buried",
                              "Hydrophobicity" = "hydro",
                              "Nucleotide" = "nucleotide"),
                  selected = "clustal"),
      
      checkboxInput("showConsensus", "Show Consensus Sequence", TRUE),
      checkboxInput("showConservation", "Show Conservation", TRUE),
      checkboxInput("showLogo", "Show Sequence Logo", TRUE),
      
      hr(),
      
      sliderInput("labelWidth", "Label Width", 
                  min = 80, max = 300, value = 150, step = 10),
      
      sliderInput("fontSize", "Font Size", 
                  min = 8, max = 16, value = 10, step = 1),
      
      hr(),
      
      downloadButton("downloadImage", "Download Visualization (PNG)"),
      downloadButton("downloadMSAR", "Download Interactive HTML")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Visualization", 
                 div(class = "content-wrapper",
                     uiOutput("msaROutput")
                 )),
        tabPanel("Sequence Information", 
                 div(class = "content-wrapper",
                     dataTableOutput("sequenceInfo")
                 )),
        tabPanel("Color Scheme Information", 
                 div(class = "content-wrapper",
                     h3("Color Scheme Details"),
                     uiOutput("colorSchemeInfo"),
                     plotOutput("colorLegend", height = "300px")
                 )),
        tabPanel("Alignment Statistics", 
                 div(class = "content-wrapper",
                     verbatimTextOutput("consensusOutput"),
                     dataTableOutput("alignmentStats")
                 ))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive value to store the alignment
  alignment <- reactiveVal(NULL)
  
  # Reactive value to store the consensus sequence
  consensus <- reactiveVal(NULL)
  
  # Description of color schemes
  colorSchemeDescriptions <- list(
    clustal = "The Clustal color scheme highlights physicochemical properties. Red: small/hydrophobic, Blue: acidic, Magenta: basic, Green: hydroxyl/amine/basic.",
    taylor = "Taylor's scheme assigns each amino acid a unique color based on distinct properties.",
    zappo = "Zappo highlights physicochemical properties. Pink: aliphatic, Orange: aromatic, Green: hydrophilic, Red: acidic, Blue: basic.",
    helix = "Colors based on alpha-helix propensity. Red indicates high propensity, blue indicates low propensity.",
    strand = "Colors based on beta-strand propensity. Red indicates high propensity, blue indicates low propensity.",
    turn = "Colors based on turn propensity. Red indicates high propensity, blue indicates low propensity.",
    buried = "Colors based on buried index. Red indicates high tendency to be buried in protein core.",
    hydro = "Colors based on hydrophobicity. Blue: hydrophilic, red: hydrophobic.",
    nucleotide = "For nucleotide sequences. Red: A, Green: T, Blue: C, Yellow: G."
  )
  
  # Function to read alignment file in different formats
  readAlignmentFile <- function(file, format) {
    if(format == "auto") {
      # Try to guess format based on file extension
      ext <- tools::file_ext(file)
      if(ext %in% c("fa", "fasta")) {
        format <- "fasta"
      } else if(ext %in% c("aln", "clustal")) {
        format <- "clustal"
      } else if(ext == "msf") {
        format <- "msf"
      } else if(ext == "stockholm") {
        format <- "stockholm"
      } else {
        # Default to FASTA if can't determine
        format <- "fasta"
      }
    }
    
    # Read alignment based on format
    if(format == "fasta") {
      return(readAAMultipleAlignment(file))
    } else if(format == "clustal") {
      # Use seqinr for Clustal format
      clustal <- read.alignment(file, format = "clustal")
      # Convert to AAMultipleAlignment
      seqs <- AAStringSet(clustal$seq)
      names(seqs) <- clustal$nam
      return(as(seqs, "AAMultipleAlignment"))
    } else if(format == "msf") {
      # Use seqinr for MSF format
      msf <- read.alignment(file, format = "msf")
      seqs <- AAStringSet(msf$seq)
      names(seqs) <- msf$nam
      return(as(seqs, "AAMultipleAlignment"))
    } else if(format == "stockholm") {
      # For Stockholm format, try to convert
      # This is simplified; may need enhancement for complex Stockholm files
      lines <- readLines(file)
      seqData <- list()
      
      current_name <- NULL
      current_seq <- ""
      
      for(line in lines) {
        if(grepl("^#", line) || line == "//") {
          # Skip comments and end marker
          next
        } else if(grepl("^\\s*$", line)) {
          # Skip empty lines
          next
        } else {
          # Sequence data line
          parts <- strsplit(trimws(line), "\\s+")[[1]]
          if(length(parts) >= 2) {
            name <- parts[1]
            seq <- parts[2]
            
            if(name %in% names(seqData)) {
              seqData[[name]] <- paste0(seqData[[name]], seq)
            } else {
              seqData[[name]] <- seq
            }
          }
        }
      }
      
      # Convert to AAStringSet
      seqs <- AAStringSet(unlist(seqData))
      names(seqs) <- names(seqData)
      return(as(seqs, "AAMultipleAlignment"))
    }
  }
  
  # Read uploaded alignment file
  observeEvent(input$alignmentFile, {
    req(input$alignmentFile)
    
    # Try to read sequences
    tryCatch({
      aligned <- readAlignmentFile(input$alignmentFile$datapath, input$fileFormat)
      alignment(aligned)
      
      # Calculate consensus sequence
      cons_seq <- msaConsensusSequence(aligned)
      consensus(cons_seq)
      
      # Show sequence information
      output$sequenceInfo <- renderDataTable({
        seqs <- as(aligned, "AAStringSet")
        seqInfo <- data.frame(
          Name = names(seqs),
          Length = width(seqs),
          GapsCount = sapply(as.character(seqs), function(s) sum(strsplit(s, "")[[1]] == "-")),
          GapsPercentage = round(sapply(as.character(seqs), function(s) 
            100 * sum(strsplit(s, "")[[1]] == "-") / nchar(s)), 2)
        )
        return(seqInfo)
      })
      
      # Show consensus sequence
      output$consensusOutput <- renderPrint({
        cat("Consensus sequence:\n")
        cat(cons_seq)
      })
      
      # Calculate and display alignment statistics
      output$alignmentStats <- renderDataTable({
        seqs <- as(aligned, "AAStringSet")
        align_width <- width(seqs)[1]
        
        # Count identical columns (all characters the same)
        seq_strings <- as.character(seqs)
        identical_count <- 0
        similar_count <- 0
        
        for(i in 1:align_width) {
          column <- sapply(seq_strings, function(s) substr(s, i, i))
          unique_chars <- unique(column[column != "-"])
          if(length(unique_chars) == 1 && length(column[column != "-"]) > 0) {
            identical_count <- identical_count + 1
          } else if(length(unique_chars) <= 2) {
            similar_count <- similar_count + 1
          }
        }
        
        # Calculate statistics
        stats <- data.frame(
          Metric = c("Number of Sequences", 
                     "Alignment Length", 
                     "Identical Positions",
                     "Identical Percentage",
                     "Similar Positions",
                     "Similar Percentage",
                     "Average Gaps Percentage"),
          Value = c(length(seqs), 
                    align_width,
                    identical_count,
                    round(100 * identical_count / align_width, 2),
                    similar_count,
                    round(100 * similar_count / align_width, 2),
                    round(mean(sapply(seq_strings, function(s) 
                      100 * sum(strsplit(s, "")[[1]] == "-") / nchar(s))), 2))
        )
        
        return(stats)
      })
      
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
    })
  })
  
  # Display colored alignment with msaR
  output$msaROutput <- renderUI({
    req(alignment())
    
    msaR(alignment(), 
         menu = TRUE, 
         overviewbox = TRUE,
         colorscheme = input$colorScheme,
         conservation = input$showConservation,
         seqlogo = input$showLogo,
         labelNameLength = input$labelWidth,
         labelid = TRUE)
  })
  
  # Display color scheme information
  output$colorSchemeInfo <- renderUI({
    scheme <- input$colorScheme
    HTML(paste("<p>", colorSchemeDescriptions[[scheme]], "</p>"))
  })
  
  # Generate color legend based on selected scheme
  output$colorLegend <- renderPlot({
    scheme <- input$colorScheme
    
    # Different amino acid categories based on scheme
    if(scheme == "clustal") {
      categories <- list(
        "Small/hydrophobic" = c("A", "V", "L", "I", "M", "F", "W"),
        "Acidic" = c("D", "E"),
        "Basic" = c("R", "K"),
        "Hydroxyl/amine/basic" = c("S", "T", "N", "Q"),
        "Glycine" = "G",
        "Proline" = "P",
        "Tyrosine" = "Y",
        "Cysteine" = "C"
      )
      colors <- c("red", "blue", "magenta", "green", "orange", "yellow", "cyan", "pink")
    } else if(scheme %in% c("zappo", "taylor", "helix", "strand", "turn", "buried", "hydro")) {
      # Simplified for demonstration - actual schemes have more detailed mappings
      amino_acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                       "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
      categories <- as.list(amino_acids)
      names(categories) <- amino_acids
      
      # Generate a color gradient for demonstration
      if(scheme == "hydro") {
        colors <- colorRampPalette(c("blue", "white", "red"))(20)
      } else {
        colors <- rainbow(20)
      }
    } else if(scheme == "nucleotide") {
      categories <- list(
        "Adenine (A)" = "A",
        "Thymine (T)" = "T",
        "Cytosine (C)" = "C",
        "Guanine (G)" = "G"
      )
      colors <- c("red", "green", "blue", "yellow")
    }
    
    # Plot legend
    par(mar = c(5, 10, 4, 2))
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, length(categories)), 
         axes = FALSE, xlab = "", ylab = "")
    
    for(i in 1:length(categories)) {
      rect(0, i-0.5, 0.2, i+0.5, col = colors[i], border = "black")
      text(0.25, i, paste(names(categories)[i], ": ", 
                          paste(categories[[i]], collapse = ", ")), 
           pos = 4, cex = 1.2)
    }
    title("Color Scheme Legend")
  })
  
  # Download handlers
  output$downloadImage <- downloadHandler(
    filename = function() {
      paste("msa-visualization-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # This is a placeholder - in a real app you would use webshot or similar
      # to capture the visualization as an image
      showNotification("PNG export feature requires additional setup with webshot package.", 
                       type = "warning")
      # For demonstration purposes, create a simple placeholder image
      png(file, width = 800, height = 600)
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
           main = "MSA Visualization (Placeholder)")
      text(1, 1, "The actual MSA visualization would be exported here")
      dev.off()
    }
  )
  
  output$downloadMSAR <- downloadHandler(
    filename = function() {
      paste("msa-visualization-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      # This is a placeholder - in a real app you would use saveWidget or similar
      # to save the msaR visualization as a standalone HTML file
      showNotification("HTML export feature requires additional setup with htmlwidgets package.", 
                       type = "warning")
      # For demonstration purposes, create a simple HTML file
      writeLines("<html><body><h1>MSA Visualization</h1><p>This is a placeholder for the interactive MSA visualization.</p></body></html>", file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)