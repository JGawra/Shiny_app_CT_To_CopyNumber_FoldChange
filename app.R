library(shiny)
library(readr)
library(dplyr)
library(tidyr)
library(openxlsx)

ui <- fluidPage(
  titlePanel("qPCR Fold Change and Copy Number Calculator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File with Raw Ct Values", accept = ".csv"),
      uiOutput("housekeeping_select"),
      uiOutput("reference_group_select"),
      downloadButton("download_results", "Download All Results (Excel)"),
      actionButton("show_about", "About / Instructions", icon = icon("info-circle"), class = "btn-info", style = "margin-bottom: 10px; width: 100%;"),      
      
    ),
    mainPanel(
      h4("Housekeeping Gene Stability (SD)"),
      tableOutput("hk_stability"),
      
      h4("ΔCt (Housekeeping - Target)"),
      tableOutput("delta_ct"),
      
      h4("Copy Number Table (2^ΔCt * 1000)"),
      tableOutput("copy_number"),
      

      h4("Fold Change (Copy number each sample by average of control for each gene)"),
      tableOutput("copy_number_fold_change"),
      
      h4("Copy Number Table (2^ΔCt * 1000)"),
      tableOutput("copy_number"),
      downloadButton("download_copy_number_csv", "Download Copy Number (CSV)"),

      h4("Fold Change (Copy number each sample by average of control for each gene)"),
      tableOutput("copy_number_fold_change"),
      downloadButton("download_fold_change_csv", "Download Fold Change (CSV)"),
      
      br(),
      div(
        HTML(
          paste(
            "Developed by Janan Gawra",
            "<a href='https://www.linkedin.com/in/janangawra/' target='_blank'>",
            "<i class='fab fa-linkedin'></i> LinkedIn Profile</a>",
            sep = " | "
          )
        ),
        style = "text-align: right; font-size: 0.9em; color: #666; margin-top: 20px;"
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Load and validate input data
  raw_data <- reactive({
    req(input$file)
    df <- read_csv(input$file$datapath, show_col_types = FALSE)
    validate(
      need(all(c("Replicate", "Treatment") %in% colnames(df)),
           "The file must include columns named 'Replicate' and 'Treatment'")
    )
    df
  })
  
  observeEvent(raw_data(), {
    df <- raw_data()
    gene_cols <- setdiff(colnames(df), c("Replicate", "Treatment"))
    updateSelectInput(session, "housekeeping", choices = gene_cols)
    updateSelectInput(session, "reference_group", choices = unique(df$Treatment))
  })
  
  output$housekeeping_select <- renderUI({
    req(raw_data())
    gene_cols <- setdiff(colnames(raw_data()), c("Replicate", "Treatment"))
    selectInput("housekeeping", "Select Housekeeping Gene", choices = gene_cols)
  })
  
  output$reference_group_select <- renderUI({
    req(raw_data())
    selectInput("reference_group", "Select Reference Group", choices = unique(raw_data()$Treatment))
  })
  output$copy_number_fold_change <- renderTable({ head(relative_copy_number_fold_change(), 10) })
  
  # Calculate HK gene stability
  hk_stability_df <- reactive({
    df <- raw_data()
    gene_cols <- setdiff(colnames(df), c("Replicate", "Treatment"))
    df %>%
      summarise(across(all_of(gene_cols), sd, na.rm = TRUE)) %>%
      pivot_longer(cols = everything(), names_to = "Gene", values_to = "SD") %>%
      arrange(SD)
  })
  
  output$hk_stability <- renderTable({ hk_stability_df() })
  
  # ΔCt calculation
  delta_ct_table <- reactive({
    df <- raw_data()
    req(input$housekeeping)
    hk <- df[[input$housekeeping]]
    gene_cols <- setdiff(colnames(df), c("Replicate", "Treatment", input$housekeeping))
    
    df %>%
      transmute(Replicate, Treatment,
                across(all_of(gene_cols), ~ hk - .x, .names = "{.col}"))
  })
  
  # Copy number table
  copy_number_table <- reactive({
    df <- delta_ct_table()
    gene_cols <- setdiff(colnames(df), c("Replicate", "Treatment"))
    df %>%
      mutate(across(all_of(gene_cols), ~ 2^.x * 1000))
  })
  
  # Relative Copy Number Fold Change (copy number divided by control mean)
  relative_copy_number_fold_change <- reactive({
    req(input$reference_group)
    df <- copy_number_table()
    gene_cols <- setdiff(colnames(df), c("Replicate", "Treatment"))
    
    # Calculate mean copy number for each gene in reference group
    control_means <- df %>%
      filter(Treatment == input$reference_group) %>%
      summarise(across(all_of(gene_cols), mean, na.rm = TRUE))
    
    # Divide each gene's copy number by the control mean
    df_out <- df
    for (gene in gene_cols) {
      df_out[[gene]] <- df[[gene]] / control_means[[gene]]
    }
    df_out
  })
  output$download_copy_number_csv <- downloadHandler(
    filename = function() {
      paste0("Copy_Number_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(copy_number_table(), file, row.names = FALSE)
    }
  )
  
  output$download_fold_change_csv <- downloadHandler(
    filename = function() {
      paste0("Fold_Change_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(relative_copy_number_fold_change(), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$show_about, {
    showModal(modalDialog(
      title = "About & Instructions",
      easyClose = TRUE,
      footer = NULL,
      size = "l",
      HTML(
        "<b>About:</b><br>
      This Shiny app calculates qPCR Fold Change and Copy Number.<br><br>
      <b>How to Use:</b>
      <ol>
        <li>Prepare your data in CSV format (see example below).</li>
        <li>Upload the CSV file using the file input on the left.</li>
        <li>Select your Housekeeping gene and Reference (Control) group.</li>
        <li>View results in the main panel, or download as Excel/CSV files.</li>
      </ol>
      <b>Input File Format:</b><br>
      <ul>
        <li>The first two columns must be: <b>Replicate</b> and <b>Treatment</b>.</li>
        <li>Each additional column must be a gene (header = gene name; values = Ct values).</li>
      </ul>
      <b>Example:</b><br>
      <pre>
Replicate,Treatment,GeneA,GeneB,GAPDH
R1,Control,24.1,25.3,21.2
R2,Control,23.8,25.1,21.0
R1,Treated,22.7,24.2,20.5
R2,Treated,22.5,24.1,20.3
      </pre>
      <b>Steps Performed:</b>
      <ol>
        <li><b>HK Stability:</b> Calculates SD of Ct for each gene.</li>
        <li><b>ΔCt:</b> For each target gene, ΔCt = Ct(Housekeeping) - Ct(Target).</li>
        <li><b>Copy Number:</b> 2^ΔCt * 1000 for each gene/sample.</li>
        <li><b>Fold Change:</b> For each gene/sample, Copy Number / Average Copy Number of Reference group.</li>
      </ol>
      "
      )
    ))
  })
  
  
  
  
  output$delta_ct <- renderTable({ head(delta_ct_table(), 10) })
  output$copy_number <- renderTable({ head(copy_number_table(), 10) })
  output$copy_number_fold_change <- renderTable({ head(relative_copy_number_fold_change(), 10) })
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("qPCR_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      wb <- createWorkbook()
      addWorksheet(wb, "Raw Input Data")
      writeData(wb, "Raw Input Data", raw_data())
      addWorksheet(wb, "HK Stability")
      writeData(wb, "HK Stability", hk_stability_df())
      addWorksheet(wb, "Delta Ct")
      writeData(wb, "Delta Ct", delta_ct_table())
      addWorksheet(wb, "Copy Number")
      writeData(wb, "Copy Number", copy_number_table())
      addWorksheet(wb, "Fold Change")
      writeData(wb, "Fold Change", relative_copy_number_fold_change())
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
}

shinyApp(ui, server)

