library(shiny)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File"),
      tags$p("The SNP positions MUST be from Genome 37"),
      tags$p("The CSV file should contain at least the following columns:"),
      tags$ul(
        tags$li("chr"),
        tags$li("start"),
        tags$li("end")
      ),
      tags$p("Example format:"),
      tableOutput("exampleTable"),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      fluidRow(
        column(8),  # Empty column to push the button to the right
        column(4, align = "right", downloadButton("downloadTable", "Download Table (CSV)")),
        column(4, align = "right", downloadButton("downloadGeneData", "Download Gene Data"))
      ),
      tableOutput("finalTable")
    )
  )
)

server <- function(input, output) {
  
  finalData <- reactiveVal()
  
  # Example table data
  exampleData <- data.frame(
    chr = c(1, 1),
    start = as.integer(c(110365045, 172864224)),
    end = as.integer(c(110365045, 172864224)),
    stringsAsFactors = FALSE
  )
  
  output$exampleTable <- renderTable({
    exampleData
  }, digits = 0)  # Ensure no decimals are shown
  
  observeEvent(input$submit, {
    inFile <- input$file
    
    if (!is.null(inFile)) {
      file_path <<- inFile$datapath
      
      showModal(modalDialog(
        title = "Processing",
        "File received. Loading results... Please wait.",
        footer = NULL
      ))
      
      source("C:\\Users\\bryan\\Desktop\\USC_research\\Huaiyu\\Snps-to-panther\\geneIDs_AnnoQ_updated.R")
      
      df <- read.csv(file_path)
      
      # Prepare columns for the annotation dataframe
      annotations_df <- data.frame(
        chr = character(),
        start = integer(),
        end = integer(),
        Ensembl_Genes = character(),
        RefSeq_Genes = character(),
        Enhancers = character(),
        rs_dbSNP151 = character(),
        stringsAsFactors = FALSE
      )
      
      # Loop through each row and process annotations
      for (i in 1:nrow(df)) {
        chr <- df$chr[i]
        start <- df$start[i]
        end <- df$end[i]
        
        annotation <- process_annotations(chr, start, end)
        annotations_df <- rbind(annotations_df, annotation)
      }
      
      final_df <- match_annotations_to_df(df, annotations_df)
      
      # Select and rename columns for display
      final_df <- final_df[, c("chr", "start", "end", "rs_dbSNP151", "Ensembl_Genes", "RefSeq_Genes", "Enhancers")]
      colnames(final_df) <- c("Chromosome", "Start", "End", "rsID", "Ensembl Gene", "RefSeq Gene", "Enhancers")
      
      finalData(final_df)
      
      removeModal()
    }
  })
  
  output$finalTable <- renderTable({
    finalData()
  })
  
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("Table-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(finalData(), file, row.names = FALSE)
    }
  )
  
  output$downloadGeneData <- downloadHandler(
    filename = function() {
      paste("Gene_Data-", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      gene_data <- finalData()[, c("Ensembl Gene", "RefSeq Gene", "Enhancers")]
      # Convert the data frame columns to a single vector
      all_values <- unlist(gene_data)
      # Remove both NA and empty strings
      all_values <- all_values[all_values != "" & !is.na(all_values)]
      # Collapsing all values into a single string separated by commas
      collapsed_values <- paste(all_values, collapse = ",")
      writeLines(collapsed_values, file)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
