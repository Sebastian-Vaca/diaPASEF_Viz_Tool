## app.R ##
library(shinydashboard)
library(shiny)
library(shinyjs)
library(shinycssloaders)
source("functions.R")

ui <- dashboardPage(
  dashboardHeader(title = "dia-PASEF Viz Tools",titleWidth = 250),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem("Inputs and Parameters", tabName = "Parameters", icon = icon("gears")),
      menuItem("Summaries", tabName = "Summaries", icon = icon("clipboard-check")),
      menuItem("Metrics", tabName = "Metrics", icon = icon("dashboard"),startExpanded = F,
               menuSubItem('Counts', tabName = "metrics_counts", icon = icon('export', lib = 'glyphicon')),
               menuSubItem('CVs', tabName = "metrics_cvs", icon = icon('bolt')),
               menuSubItem('QC metrics', tabName = "metrics_QCs", icon = icon("bolt")),
               menuSubItem('Dynamic range', tabName = "metrics_dynamicRange", icon = icon("bolt")),
               menuSubItem('Correlations/PCA', tabName = "metrics_correlations", icon = icon("bolt")),
               menuSubItem('IM vs m/z heatmap', tabName = "metrics_IM_mz_map", icon = icon("bolt"))
               ),
      menuItem("Upset plots", tabName = "Upset_plots", icon = icon("stats", lib= "glyphicon")),
      menuItem("PTMs", tabName = "PTMs", icon = icon("flask")),
      menuItem("Dilution series", tabName = "DilSeries", icon = icon("vial")),
      menuItem("Targeted Peptide/Protein Tools", tabName = "PRM", icon = icon("bullseye", lib = "font-awesome")),
      menuItem("GCT Builder", tabName = "GCT", icon = icon("table")),

      tags$hr(),
      checkboxInput("reassign_checkbox_input", "Reassign/Reorder groups", T),
      checkboxInput("remove_checkbox_input", "Remove selected runs", F),
      
      tags$hr()
      

    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "Parameters",
              shinyjs::useShinyjs(),
              titlePanel("dia-PASEF Visualization Tools (v23.11.1)"),
              radioButtons("radio_input_software",
                           label = h3("Which software?"),
                           choices = list("Spectronaut" = "Spectronaut",
                                          "DIANN" = "DIANN",
                                          "timsDIANN" = "timsDIANN"), 
                           selected = 1),
              
              fluidRow(column(10, textOutput("value"))),
              
              # Input: Select a file ----
              hidden(textOutput("Input_file_text"),
              tags$style("#Input_file_text {font-size:25px;
               color:black;
               display:block; }")),
              
              hidden(
                fileInput("file1", "Choose input file",
                          multiple = TRUE,
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv",
                                     "text/tsv",
                                     ".xls",
                                     ".tsv"))),
              
              # Input: Checkbox if file has header ----
              # h3("Metadata"),
              hidden(textOutput("Metadata_file_input_text"),
                     tags$style("#Metadata_file_input_text {font-size:25px;
               color:black;
               display:block; }")),
              
              
              # Metadata
              hidden(textOutput("Metadata_file_input_text_2"),
                     tags$style("#Metadata_file_input_text_2 {font-size:14px;
               color:black;
               display:block; }")),
              hidden(
                downloadButton(outputId = "download_metadata_template", label = "Download Metadata Template")),
              
              # Input: Select a file ----
              hidden(fileInput("groups_id_csv", "Add metadata file",
                        multiple = TRUE,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv"))),
              br(),
              
              actionButton(inputId = "DataFiltering_ActionButton", align = "left",label = "Update Filtering Criteria", style="color: #FFFFFF; background-color: #0071BC"),
              
              
              # hidden(tabBox(id = "DataFiltering_tabBox",
              #                title = tagList(shiny::icon("gear"), "Filtering parameters"), width =9,
              #   tabPanel("q-values", 
              #            checkboxInput("DataFiltering_EG_qValue_checkbox", "Elution Group q-value", value = T),
              #            textInput("DataFiltering_EG_qValue", "Max q-value:", value = 0.01),
              #            checkboxInput("DataFiltering_PG_qValue_checkbox", "Protein Group q-value", value = T),
              #            textInput("DataFiltering_PG_qValue", "Max q-value:", value = 0.01)),
              #   tabPanel("Peptide Length", 
              #            checkboxInput("DataFiltering_PeptideLength_checkbox", "Peptide Length", value = F),
              #            textInput("DataFiltering_PeptideLength", "min peptideLength:", value = 6),
              #            textInput("DataFiltering_PeptideLength", "min peptideLength:", value = 50)),
              #   tabPanel("Proteotypicity", 
              #            checkboxInput("DataFiltering_Proteotypicity_checkbox", "Keep only Proteotypic peptides", value = F)),
              # actionButton(inputId = "DataFiltering_ActionButton", align = "left",label = "Update Filtering Criteria", style="color: #FFFFFF; background-color: #0071BC"))
              # )
              
              uiOutput("DataFiltering_tabBox")
              
              
              # img(src = "diaQuito2.png", height = 150, width = 150),
              # img(src = "shiny.PNG", height = 50, width = 50)
              
      ),
      # Summaries
      tabItem(tabName = "Summaries",
              tabsetPanel(type = "tabs",
                          
                          tabPanel("Metadata", tableOutput("metadata_table")),
                          tabPanel("Column_check", tableOutput("column_check")),
                          tabPanel("Reports",
                                   useShinyjs(),
                                   h4("Cofficients of variation"),
                                   h6("Download CV data (precursors)"),
                                   uiOutput("downloadRender_CVs_data"),
                                   # downloadButton(outputId = "Download_CVs_data", label = "Download CV data (precursors)"),
                                   h6("Download CV data (protein groups)"),
                                   uiOutput("downloadRender_CVs_proteins_data"),
                                   # downloadButton(outputId = "Download_CVs_proteins_data", label = "Download CV data (protein groups)"),
                                   br(),
                                   h4("Peptide and proteins counts"),
                                   h6("Download count data (per MS run)"),
                                   uiOutput("downloadRender_counts_per_MSRun_data"),
                                   # downloadButton(outputId = "Download_counts_per_MSRun_data", label = "Download count data (per MS run)"),
                                   h6("Download count data (per condition)"),
                                   uiOutput("downloadRender_counts_per_condition"),
                                   # downloadButton(outputId = "Download_counts_per_condition", label = "Download count data (per condition)")
                                   ),
                          tabPanel("Table_head", tableOutput("contents"))
              )
      ),
      #Metrics
      tabItem(tabName = "Metrics",
              tabsetPanel(type = "tabs")
              ),
    
      tabItem(tabName = "metrics_counts",
              tabsetPanel(type = "tabs",
                          tabPanel("Counts per MS Run",
                                   selectInput(inputId = "countsPerRun_plot_type", label ="What do you want to plot? [Leave blank to see all]",
                                               choices = c("Precursors","Peptides","Proteins","Proteins_with2pepts"), multiple = T),
                                   plotOutput("Counts_per_MSRun")),
                          tabPanel("Unique counts per Condition", plotOutput("Counts")),
                          tabPanel("Mean counts per Condition",
                                   selectInput(inputId = "CountsErrorBars_plot_type", label ="What do you want to plot? [Leave blank to see all]",
                                               choices = c("Precursors","Peptides","Proteins","Proteins_with2pepts"), multiple = T),
                                   plotOutput("CountsErrorBars")),
                          tabPanel("Data Completeness", plotOutput("data_completeness")),
                          tabPanel("Peptides per Protein",
                                   selectInput(inputId = "PeptidesperProtein_item", label ="What do you want to plot?",
                                               choices = c("Precursors","Peptides"), multiple = F),
                                   selectInput(inputId = "PeptidesperProtein_plotType", label ="How do you want to plot it?",
                                               choices = c("histogram","boxplot"), multiple = F),
                                   selectizeInput(inputId = "PeptidesperProtein_item_conditions", label ="Select Conditions [Leave blank to see all]",
                                                  choices = NULL, multiple = T, options = list(maxItems = 10)),
                                   checkboxInput("PeptidesperProtein_outliers_checkbox", "Remove outliers (95th percentile)", value = T),
                                   actionButton(inputId = "PeptidesperProtein_actionbutton",label = "Update plot", style="color: #FFFFFF; background-color: #0071BC"),
                                   hidden(div(id= "PeptidesperProtein_hider", plotOutput("PeptidesperProtein_plot") %>% withSpinner()) )) 
                          )
      ),
      tabItem(tabName = "metrics_cvs",
              tabsetPanel(type = "tabs",
                          tabPanel("CV distribution",
                                   selectizeInput(inputId = "CVs_conditions", label ="Select Conditions [Leave blank to see all]",
                                                  choices = "Null", multiple = T, options = list(maxItems = 5)),
                                   selectInput(inputId = "CVs_level", label ="Select Protein/Peptide level",
                                               choices = c("Protein", "Precursor"), multiple = F),
                                   selectInput(inputId = "CVs_plot_type", label ="Select plot type",
                                               choices = c("freqpoly", "histogram","boxplot"), multiple = F),
                                   checkboxInput("CVs_checkbox", "Remove outliers (95th percentile)", F),
                                   actionButton(inputId = "CVs_actionbutton",label = "Update plot", style="color: #FFFFFF; background-color: #0071BC"),
                                   hidden(div(id= "CVs_Plot_hider", plotOutput("CVs") %>% withSpinner() ))
                                   ),
                          tabPanel("CV cutoff: Precursor level",
                                   fluidRow(
                                     selectizeInput(inputId = "CVs_cutoff_precursors_conditions", label ="Select Conditions [Leave blank to see all]",
                                                    choices = "Null", multiple = T, options = list(maxItems = 10)),
                                     selectInput(inputId = "CVs_cutoff_precursors_type", label ="What would you like to plot? [Leave blank to see all]",
                                                 choices = c("Identified","Quantified","CV < 40","CV < 20","CV < 10"), multiple = T),
                                     actionButton(inputId = "CVs_cutoff_precursors_actionbutton",label = "Update plot", style="color: #FFFFFF; background-color: #0071BC"),
                                     hidden(div(id= "CVs_cutoff_precursors_Plot_hider", plotOutput("CVs_cutoff_Precursors") %>% withSpinner())))
                                   ),
                          tabPanel("CV cutoff: Protein Group level",
                                   fluidRow(
                                     selectizeInput(inputId = "CVs_cutoff_proteins_conditions", label ="Select Conditions [Leave blank to see all]",
                                                    choices = "Null", multiple = T, options = list(maxItems = 10)),
                                     selectInput(inputId = "CVs_cutoff_proteins_type", label ="What would you like to plot? [Leave blank to see all]",
                                                 choices = c("Identified","Quantified","CV < 40","CV < 20","CV < 10"), multiple = T),
                                     actionButton(inputId = "CVs_cutoff_proteins_actionbutton",label = "Update plot", style="color: #FFFFFF; background-color: #0071BC"),
                                     hidden(div(id= "CVs_cutoff_proteins_Plot_hider", plotOutput("CVs_cutoff_proteins") %>% withSpinner()))))
                          )
              ),
      
      tabItem(tabName = "metrics_dynamicRange",
              tabsetPanel(type = "tabs",
                          tabPanel("Dynamic range",
                                   fluidRow(
                                     selectInput("waterfall_condition", "Select condition", choices = "Null"),
                                     selectInput("waterfall_quant_value_to_use", "Select: Peptide or Protein quant", choices = c("PG.MS2Quantity", "FG.MS2Quantity")),
                                     selectInput("waterfall_mean_or_median", "Select: mean or median", choices = c("mean", "median")),
                                     plotOutput("watefall_plot_simplified")
                                   )),
                          )
      ),
      
      tabItem(tabName = "metrics_QCs",
              tabsetPanel(type = "tabs",
                          tabPanel("Metric distributions",
                                   fluidRow(
                                     selectizeInput(inputId = "BoxHist_conditions", label ="Select Conditions [Leave blank to see all]",
                                                    choices = "Null", multiple = T, options = list(maxItems = 5)),
                                     selectizeInput(inputId = "BoxHist_metric_to_plot", label ="Select Metric",
                                                    choices = NULL, multiple = F),
                                     selectInput(inputId = "BoxHist_plot_type", label ="Select plot type",
                                                 choices = c("freqpoly", "histogram","boxplot"), multiple = F),
                                     checkboxInput("BoxHist_checkbox", "Remove outliers (5th and 95th percentile)", T),
                                     actionButton(inputId = "BoxHist_actionbutton",label = "Update plot", style="color: #FFFFFF; background-color: #0071BC"),
                                     hidden(div(id= "BoxHist_Plot_hider", plotOutput("BoxHist_Plot") %>% withSpinner() ))
                                   )),
                          tabPanel("Peptide lenghts",
                                   fluidRow(
                                     selectizeInput(inputId = "PeptideLength_conditions", label ="Select Conditions [Leave blank to see all]",
                                                    choices = "Null", multiple = T, options = list(maxItems = 10)),
                                     selectInput(inputId = "PeptideLength_plot_type", label ="Select plot type",
                                                 choices = c("freqpoly", "histogram"), multiple = F),
                                     actionButton(inputId = "PeptideLength_actionbutton",label = "Update plot", style="color: #FFFFFF; background-color: #0071BC"),
                                     hidden(div(id= "PeptideLength_Plot_hider", plotOutput("PeptideLength_Plot") %>% withSpinner() ))
                                   )),
                          tabPanel("Charge States",
                                   fluidRow(
                                     selectizeInput(inputId = "ChargeStates_conditions", label ="Select Conditions [Leave blank to see all]",
                                                    choices = "Null", multiple = T, options = list(maxItems = 10)),
                                     actionButton(inputId = "ChargeStates_actionbutton",label = "Update plot", style="color: #FFFFFF; background-color: #0071BC"),
                                     hidden(div(id= "ChargeStates_Plot_hider", plotOutput("ChargeStates_Plot") %>% withSpinner() ))
                                   )),
                          tabPanel("TICs of ID peptides", plotlyOutput("plot_TICS", height = "800px", width = "800px")),
                          
                          )
      ),
      
      tabItem(tabName = "metrics_correlations",
              tabsetPanel(type = "tabs",
                          tabPanel("Correlations",
                                   fluidRow(
                                     selectizeInput(inputId = "CorrPlot_conditions", label ="Select Conditions [Leave blank to see all]",
                                                    choices = "Null", multiple = T, options = list(maxItems = 5)),
                                     selectInput(inputId = "CorrPlot_metric_to_plot", label ="Select Metric",
                                                 choices = c("PG.MS2Quantity","FG.MS2Quantity", "EG.IonMobility", "EG.ApexRT"), multiple = F),
                                     checkboxInput(inputId = "CorrPlot_fillUpper_checkbox",label = "Color upper triangle"),
                                     actionButton(inputId = "CorrPlot_actionbutton",label = "Update plot", style="color: #FFFFFF; background-color: #0071BC"),
                                     hidden(div(id ="Corr_Plot_hider", plotOutput("Corr_Plot") %>% withSpinner() ))
                                   )),
                          tabPanel("PCA",
                                   fluidRow(
                                     actionButton(inputId = "PCA_actionbutton",label = "Update PCA plot", style="color: #FFFFFF; background-color: #0071BC"),
                                     hidden(div(id ="PCA_hider_scree", plotlyOutput("PCA_Plot_scree", height = "200px", width = "500px") %>% withSpinner() )),
                                     hidden(div(id ="PCA_hider_pca", plotlyOutput("PCA_Plot_pca", height = "400px", width = "500px") %>% withSpinner() ))
                                     ))
                          )
      ),
      tabItem(tabName = "metrics_IM_mz_map",
              tabsetPanel(type = "tabs",
                          tabPanel("IM vs m/z heatmap",
                                   fluidRow(
                                     selectInput("im_mz_map_condition", "Select condition", choices = "Null"),
                                     selectInput("im_mz_map_replicate", "Select replicate", choices = "Null"),
                                     numericInput("im_mz_map_upper_limit_counts", "Set upper limit for the counts", value = NULL),
                                     numericInput("im_mz_map_binwidth_mz", "Set width for m/z (Da)", value = 50),
                                     numericInput("im_mz_map_binwidth_im", "Set width for ion-mobility (1/k0)", value = 0.01),
                                     sliderInput("im_mz_map_sliderIM", label = h3("Slider Range"), min = 0.6, 
                                                 max = 1.7, value = c(0.7, 1.3),step = 0.05),
                                     sliderInput("im_mz_map_sliderMZ", label = h3("Slider Range"), min = 200, 
                                                 max = 1600, value = c(350, 1200),step = 25),
                                     actionButton(inputId = "im_mz_map_actionbutton",label = "Refresh plot", style="color: #FFFFFF; background-color: #0071BC"),
                                     hidden(div(id="im_mz_map_plot_hider",plotOutput("im_mz_map_plot") %>% withSpinner(type = 8, color = "#0071BC", size = 2)))
                                   ))
                          )
      ),
      # PTMs
      tabItem(tabName = "PTMs",
              tabsetPanel(type = "tabs",
                          tabPanel("PTMS: modifications found",
                                   selectInput("ptm_counts_selectinput", "Select modification", choices = "Null"),
                                   tableOutput("Modifications_found_table")),
                          tabPanel("PTMs: Mean counts per replicate",
                                   selectInput(inputId = "PTMcounts_plot_type", label ="What do you want to plot? [Leave blank to see all]",
                                               choices = c("Precursors","Peptides","Proteins","Proteins_with2pepts"), multiple = T),
                                   plotOutput("PTM_plot_Counts")),
                          tabPanel("PTMs: Enrichment yield",
                                   selectInput(inputId = "PTMenrich_plot_type", label ="What do you want to plot?",
                                               choices = c("Precursors","Peptides","Proteins","Proteins_with2pepts [Leave blank to see all]"), multiple = T),
                                   plotOutput("PTM_EnrichYield")),
                          tabPanel("PTMs: CVs", plotOutput("PTM_CVs"))
              )
      ),
      # Upset Plots
      tabItem(tabName = "Upset_plots",
              tabsetPanel(type = "tabs",
                          tabPanel("Proteins with condition filter",
                                   selectizeInput(inputId = "upset_Prots_withConditionsFilter_conditions", label ="Select Conditions [Leave blank to see all]",
                                                  choices = NULL, multiple = T, options = list(maxItems = 10)),
                                   actionButton(inputId = "upset_plot_prots_actionbutton",label = "Refresh plot", style="color: #FFFFFF; background-color: #0071BC"),
                                   h6("Download upset matrix"),
                                   hidden(div(id="upset_plot_prots_hider_1",
                                              uiOutput("upset_plot_prots_downloadRender") %>% 
                                                withSpinner(type = 2, color.background =  "#0071BC"))),
                                   hidden(div(id="upset_plot_prots_hider_2",
                                              plotOutput("upset_plot_prots_withConditionsFilter") %>%
                                                withSpinner() ))
                                   ),
                          tabPanel("Precursors with condition filter",
                                   selectizeInput(inputId = "upset_Precs_withConditionsFilter_conditions", label ="Select Conditions [Leave blank to see all]",
                                                  choices = NULL, multiple = T, options = list(maxItems = 10)),
                                   actionButton(inputId = "upset_plot_precs_actionbutton",label = "Refresh plot", style="color: #FFFFFF; background-color: #0071BC"),
                                   h6("Download upset matrix"),
                                   hidden(div(id="upset_plot_precs_hider_1",
                                              uiOutput("upset_plot_precs_downloadRender") %>% 
                                                withSpinner(type = 2, color.background =  "#0071BC"))),
                                   hidden(div(id="upset_plot_precs_hider_2",
                                              plotOutput("upset_plot_precs_withConditionsFilter") %>%
                                                withSpinner() ))
                                   )
                                    
              )
      ),
      ## Dilution series
      tabItem(tabName = "DilSeries",
              tabsetPanel(type = "tabs",
                          tabPanel("Dilution series: Parameters",
                                   # Horizontal line ----
                                   tags$hr(),
                                   
                                   h3("Dilution series experiment"),
                                   
                                   selectInput("dil_series_Conditions", "Select normalizing condition", choices = "Null"),
                                   verbatimTextOutput("normalizing_cond_value2"),
                                   fluidRow(column(numericInput("plot_height", "Height",value = 800),width = 2),
                                            column(numericInput("plot_width", "Width",value = 800),width = 2)),
                                   plotOutput("calcurve_normalized2"))
                          
              )
      ),
      ## PRM
      tabItem(tabName = "PRM",
              tabsetPanel(type = "tabs",
                          tabPanel("Load list of Peptide and Protein Targets", 
                                   h3("Load CSV files containing the peptide/Protein targets"),
                                   fileInput("Targets_input", "Choose Peptide targets input file",
                                             accept = c("text/csv",
                                                        "text/comma-separated-values,text/plain",
                                                        ".csv",
                                                        "text/tsv",
                                                        ".xls",
                                                        ".tsv")),
                                   h3("Download list of all observed Protein and Peptides"),
                                   downloadButton(outputId = "Download_list_of_Observed_Proteins_n_Peptides", label = "Download Protein/Peptide list")
                                   ),
                          tabPanel("Targets: Mean counts per replicate", plotOutput("plot_counts_PeptideTargets")),
                          tabPanel("Targets: Mean counts per replicate", plotOutput("plot_counts_Protein_Targets")),
                          tabPanel("Targets: Calibration curves",
                                   h3("Select normalizing condition"),
                                   selectInput("Calcurve_NormCondition_selection", "Select normalizing condition", choices = "Null"),
                                   selectInput("Calcurve_ModPept_selection", "Select Modified peptide", choices = "Null"),
                                   plotOutput("cal_curve_norm_to_single_condition_singlePept")),
                          tabPanel("Peptide Targets: Intensity accross runs", plotlyOutput("plotly_Peptide_Intensity_accross_runs")),
                          tabPanel("Peptide Targets: HeatMap log2(ratios)", plotlyOutput("plot_HeatMap_Peptide_Targets_ratios_by_condition")),
                          tabPanel("Peptide Targets: HeatMap log10(intensity)", plotlyOutput("plot_HeatMap_Peptide_Targets")),
                          
                          tabPanel("Protein Targets: Intensity accross runs", plotlyOutput("plotly_Protein_Intensity_accross_runs")),
                          tabPanel("Protein Targets: HeatMap log2(ratios)", plotlyOutput("plot_HeatMap_Protein_Targets_ratios_by_condition")),
                          tabPanel("Protein Targets: HeatMap log10(intensity)", plotlyOutput("plot_HeatMap_Protein_Targets")),
                          tabPanel("Dynamic range",
                                   fluidRow(
                                     selectInput("waterfall_condition", "Select condition", choices = "Null"),
                                     selectInput("waterfall_quant_value_to_use", "Select: Peptide or Protein quant", choices = c("PG.MS2Quantity", "FG.MS2Quantity")),
                                     selectInput("waterfall_mean_or_median", "Select: mean or median", choices = c("mean", "median")),
                                   checkboxInput(inputId = "waterfall_checkbox_input",label = "Highlight targets"),
                                   # selectizeInput(inputId = "waterfall_targets", label ="Select targets ProteinGroup",
                                   #                choices = "Null", multiple = T, options = list(maxItems = 10)),
                                   actionButton(inputId = "waterfall_actionbutton",label = "Plot!", style="color: #FFFFFF; background-color: #0071BC"),
                                   plotOutput("watefall_plot")
                                   )),
                          tabPanel("Create PRM Method", tableOutput("format_prm_method"))
                          
              )
      ),
      ## GCT
      tabItem(tabName = "GCT",
              tabsetPanel(type = "tabs",
                          tabPanel("GCT Builder",
                                   h3("GCT tools: More features coming soon!"),
                                   downloadButton(outputId = "download_gct_PGLevel", label = "Download GCT file (Protein Level)")
                          )
              )
      )
    ),
    
    ## Appearance
    tags$head(
      tags$style(HTML(".skin-blue .main-sidebar {background-color:  #004E82;}
                      .skin-blue .main-header .logo {background-color: #0071BC;}
                      .skin-blue .main-header .navbar {background-color: #0071BC;}"
                      
                      )))
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  options(shiny.maxRequestSize = 15000*1024^2)
  
  ### Parameters
  
  observeEvent(input$radio_input_software, {
    toggle("file1",condition = T)
  })
  
  observeEvent(input$radio_input_software, {
    toggle("Input_file_text",condition = T)
  })
  observeEvent(input$radio_input_software, {
    output$Input_file_text <- renderText({"Input file"})
  })
  
  observeEvent(input$file1, {
    toggle("Metadata_file_input_text",condition = T)
  })
  observeEvent(input$file1, {
    output$Metadata_file_input_text <- renderText({"Metadata file"})
  })
  
  observeEvent(input$file1, {
    toggle("Metadata_file_input_text_2",condition = T)
  })
  observeEvent(input$file1, {
    output$Metadata_file_input_text_2 <- renderText({"If you don't have already a Metadata file you can download it below after the upload is complete."})
  })
  
  
  observeEvent(input$file1, {
    toggle("download_metadata_template",condition = T)
  })
  observeEvent(input$file1, {
    toggle("groups_id_csv",condition = T)
  })
  
  # observeEvent(input$file1, {
  #   toggleElement("DataFiltering_tabBox", condition = T)
  # })

  output$DataFiltering_tabBox <- renderUI({
    req(input$groups_id_csv)
    
    tabBox(id = "DataFiltering_tabBox",
           title = tagList(shiny::icon("gear"), "Filtering parameters"), width =9,
           tabPanel("q-values",
                    checkboxInput("DataFiltering_EG_qValue_checkbox", "Elution Group q-value", value = T),
                    textInput("DataFiltering_EG_qValue", "Max q-value:", value = 0.01),
                    checkboxInput("DataFiltering_PG_qValue_checkbox", "Protein Group q-value", value = T),
                    textInput("DataFiltering_PG_qValue", "Max q-value:", value = 0.01)),
           tabPanel("Peptide Length",
                    checkboxInput("DataFiltering_PeptideLength_checkbox", "Peptide Length", value = F),
                    textInput("DataFiltering_PeptideLength_MinValue", "Min peptideLength:", value = 6),
                    textInput("DataFiltering_PeptideLength_MaxValue", "Max peptideLength:", value = 50)),
           tabPanel("Proteotypicity",
                    checkboxInput("DataFiltering_Proteotypicity_checkbox", "Keep only Proteotypic peptides", value = F)))
  })

  
  
  
  ## UI
  output$value <- renderText({
    a <-  ifelse(is.null(input$radio_input_software), "Please choose the software",
                 switch(input$radio_input_software, 
                        "DIANN" = "DIANN: requires the output .tsv file and a metadata file that can be downloaded below.",
                        "timsDIANN"= "timsDIANN: requires the output .tsv file and a metadata file that can be downloaded below.",
                        "Spectronaut" = "Spectronaut: requires a report created witht the schema that can be downloaded on our GitHub page.")
    )
    a
  })
  
  ########
  #####
  ## Global datasets
  ####
  
  
  data_input_unfiltered <- reactive({
    
    
    req(input$file1)
    req(input$radio_input_software)
    
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    df <- load_data_2(raw_input_file_path = input$file1$datapath,
                      metadata_filepath = input$groups_id_csv$datapath,
                      software = software_used,
                      needs_reassign = input$reassign_checkbox_input)
    df
  })
  
  data_column_verif <- reactive({
    req(data_input_unfiltered())
    
    load_data_input_verification_afterLoad(dt = data_input_unfiltered())
  })
  
  output$column_check <- renderTable({
    return(data_column_verif())
  })
  
  
  
  data_input <- reactiveVal(NULL)
  
  observeEvent(input$DataFiltering_ActionButton, {
    if (!is.null(data_input_unfiltered())) {
    
      req(data_input_unfiltered())
      df = data_input_unfiltered()
      
      # df <- df %>% remove_data_based_on_metadata(remove_selected_runs = input$remove_checkbox_input)
      
      df <- load_data_filter_data(dt_unfiltered = df,
                                  filter_EG_qValue = input$DataFiltering_EG_qValue_checkbox,
                                  EG_qValue_cutoff = as.numeric(input$DataFiltering_EG_qValue),
                                  filter_PG_qValue = input$DataFiltering_PG_qValue_checkbox,
                                  PG_qValue_cutoff = as.numeric(input$DataFiltering_PG_qValue),
                                  filter_PeptideLenght = input$DataFiltering_PeptideLength_checkbox,
                                  PeptideLenght_min_cutoff = as.numeric(input$DataFiltering_PeptideLength_MinValue),
                                  PeptideLenght_max_cutoff = as.numeric(input$DataFiltering_PeptideLength_MaxValue),
                                  filter_isProteotypic = input$DataFiltering_Proteotypicity_checkbox)
      
      data_input(df)
      
      }
    
  })

  
  metadata_input <- reactive({
    req(input$groups_id_csv)
    group_dt <- fread(input$groups_id_csv$datapath)
    group_dt
  })
  
  metadata_template <- reactive({
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    metadata_template2 <- create_metadata_template(
      input_filepath = input$file1$datapath,
      software = software_used)
    
    metadata_template2
  })
  
  Targets_input <- reactive({
    req(input$Targets_input)
    PeptideTargets_dt <- load_Targets(input$Targets_input$datapath)
    PeptideTargets_dt
  })
  
  ##### Summaries and table outputs

  output$downloadRender_CVs_data <- renderUI({
    req(CVs_data())
    downloadButton("Download_CVs_data")
  })
  
  output$Download_CVs_data <- downloadHandler(
    
    filename = function() {
      paste("CVs_precursor_level-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      CVs_data = CVs_data()
      
      fwrite(CVs_data, file, sep = ",")
    })
  
  output$downloadRender_CVs_proteins_data <- renderUI({
    req(CVs_data())
    downloadButton("Download_CVs_proteins_data")
  })
  
  output$Download_CVs_proteins_data <- downloadHandler(
    
    filename = function() {
      paste("CVs_protein_level-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      CVs_data = CVs_proteins_data()
      
      fwrite(CVs_data, file, sep = ",")
    })
  
  output$downloadRender_counts_per_MSRun_data <- renderUI({
    req(CVs_data())
    downloadButton("Download_counts_per_MSRun_data")
  })
  output$Download_counts_per_MSRun_data <- downloadHandler(
    
    filename = function() {
      paste("counts_per_MSRun-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      counts = counts_per_MSRun_data()
      
      fwrite(counts, file, sep = ",")
    })
  
  output$downloadRender_counts_per_condition <- renderUI({
    req(CVs_data())
    downloadButton("Download_counts_per_condition")
  })
  
  output$Download_counts_per_condition <- downloadHandler(
    
    filename = function() {
      paste("counts_per_MSRun-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      counts = counts_data()
      
      fwrite(counts, file, sep = ",")
    })
  
  
  ##CVs
  CVs_data <- reactive({
    
    CVs <- Calculate_CVs(data_input())
    
    CVs
  })
  
  CVs_proteins_data <- reactive({
    
    CVs <- Calculate_protein_CVs(data_input())
    
    CVs
  })
  counts_per_MSRun_data<- reactive({
    
    counts <- create_counts_table_per_MSRun(data_input())
    
    counts
  })
  counts_data <- reactive({
    
    counts <- create_counts_table(data_input())
    
    counts
  })
  
  
  #####
  ## Plots
  ####
  
  output$download_metadata_template <- downloadHandler(
    
    filename = function() { 
      paste("Metadata_template-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      fwrite(metadata_template(), file, sep = ",")
    })
  
  output$contents <- renderTable({
    
    df <- data_input()
    
    return(head(df))
  })
  
  
  
  output$contents <- renderTable({
    
    return(head(data_input()))
  })
  
  output$metadata_table <- renderTable({
    return(metadata_input())
  })
  
  ### Upset Plots
  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "upset_Prots_withConditionsFilter_conditions", choices = unique(metadata_input()[,"R.Condition"]))
    updateSelectInput(inputId = "upset_Precs_withConditionsFilter_conditions", choices = unique(metadata_input()[,"R.Condition"]))
    
  })
  #### Upset Prots

  Upset_protein_matrix <- eventReactive(input$upset_plot_prots_actionbutton, {
    
    Upset_protein_matrix <- upset_plot_from_spectronaut_prots_matrix(
      data_input(),
      input$upset_Prots_withConditionsFilter_conditions)
    
  })
  
  observeEvent(input$upset_plot_prots_actionbutton, {
    show("upset_plot_prots_hider_1")
    toggle(id = "upset_plot_prots_downloadRender", condition = T)
    show("upset_plot_prots_hider_2")
    toggle(id = "upset_plot_prots_withConditionsFilter", condition = T)
    
    output$upset_plot_prots_withConditionsFilter <- renderPlot({
      
      isolate({
        G = upset_plot_from_spectronaut_prots_withConditionsFilter(data_input(),
                                                                   input$upset_Prots_withConditionsFilter_conditions)
      })
      
      G
    })
  })

  observeEvent(input$upset_plot_prots_actionbutton, {
    output$upset_plot_prots_downloadRender <- renderUI({
      isolate({
        req(Upset_protein_matrix())
        downloadButton("upset_plot_prots_Download")
      })
    })
  })
  
  output$upset_plot_prots_Download <- downloadHandler(
    
    filename = function() {
      paste("Upset_proteins-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      counts = Upset_protein_matrix()
      
      fwrite(counts, file, sep = ",")
    })
  

  #### Upset Precursors
  Upset_precursors_matrix <- eventReactive(input$upset_plot_precs_actionbutton, {

      Upset_precursors_matrix <- upset_plot_from_spectronaut_precs_matrix(
        data_input(),
        input$upset_Precs_withConditionsFilter_conditions)

    })
  
  observeEvent(input$upset_plot_precs_actionbutton, {
    show("upset_plot_precs_hider_1")
    toggle(id = "upset_plot_precs_downloadRender", condition = T)
    show("upset_plot_precs_hider_2")
    toggle(id = "upset_plot_precs_withConditionsFilter", condition = T)
    
    output$upset_plot_precs_withConditionsFilter <- renderPlot({
    
      isolate({
      G = upset_plot_from_spectronaut_precs_withConditionsFilter(data_input(),
                                                             input$upset_Precs_withConditionsFilter_conditions)
      })
      
      G
    })
  })
  
  observeEvent(input$upset_plot_precs_actionbutton, {
    output$upset_plot_precs_downloadRender <- renderUI({
      isolate({
        req(Upset_precursors_matrix())
       downloadButton("upset_plot_precs_Download")
      })
    })
  })
  
output$upset_plot_precs_Download <- downloadHandler(

    filename = function() {
      paste("Upset_precursors-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      counts = Upset_precursors_matrix()
      
      fwrite(counts, file, sep = ",")
})
    
  

  
  
  ## Counts
  
  output$Counts_per_MSRun <- renderPlot({
    # 
    # plot_count_per_MSRun_pepts_prots_precs_reordered(counts_per_MSRun_data(),
    #                                        data_input(),
    #                                        input$reassign_checkbox_input )
    plot_count_per_MSRun_pepts_prots_precs_reordered_filter(counts_per_MSRun_data(),
                                                            data_input(),
                                                            plot_type = input$countsPerRun_plot_type,
                                                            input$reassign_checkbox_input )
  })
  
  
  output$Counts <- renderPlot({
    
    plot_count_pepts_prots_precs_reordered(counts_data(),
                                           data_input(),
                                           input$reassign_checkbox_input )
  })
  
  output$CountsErrorBars <- renderPlot({
    
    # plot_count_pepts_prots_precs_withErrorBars_reordered(counts_data(),
    #                                                      data_input(),
                                                         # input$reassign_checkbox_input )
    plot_count_pepts_prots_precs_withErrorBars_reordered_filter(counts_data(),
                                                                data_input(),
                                                                plot_type = input$CountsErrorBars_plot_type,
                                                                input$reassign_checkbox_input )
    
  })
  
  output$plot_TICS <- renderPlotly({
    
    plot_TICS(data_input())
    
  })
  
  # Peptides per Protein
  
  Precursors_per_Protein_dt <- reactive({
    req(input$file1)
    req(input$groups_id_csv)
    
    Precursors_per_Protein_dt <- Precursors_per_Protein(data_input())
    Precursors_per_Protein_dt
  })
  
  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "PeptidesperProtein_item_conditions", choices = unique(metadata_input()[,"R.Condition"]))
  })
  
  observeEvent(input$PeptidesperProtein_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("PeptidesperProtein_hider")
    toggle(id = "PeptidesperProtein_plot", condition = T)
    
    output$PeptidesperProtein_plot<- renderPlot({
      
      isolate({
        G = plot_number_of_peptides(Precursors_per_Protein_dt = Precursors_per_Protein_dt(),
                                    item = input$PeptidesperProtein_item,
                                    conditions = input$PeptidesperProtein_item_conditions,
                                    remove_outliers = input$PeptidesperProtein_outliers_checkbox,
                                    plot_type = input$PeptidesperProtein_plotType)
          
      })
      
      G
      
    })
  })
  
  output$watefall_plot_simplified <- renderPlot({
      
      
      G = plot_and_table_waterfall(dt = data_input(),
                                   condition = input$waterfall_condition, 
                                   mean_or_median = input$waterfall_mean_or_median,
                                   quant_value_to_use = input$waterfall_quant_value_to_use,
                                   highlight_targets = F,
                                   list_of_targetProteinGroup = NULL
      )
      
      G
      
    },width = 400,height = 600)
  
  ### im vs. mz map
  im_mz_map_replicates_available <- reactive({
    
    req(input$im_mz_map_condition)
    selected_condition = as.character(input$im_mz_map_condition)

    dt = metadata_input() %>%
      filter(.data[["R.Condition"]] %in% selected_condition)

    unique(dt$R.Replicate)
    

  })
  
  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "im_mz_map_condition", choices = unique(metadata_input()[,"R.Condition"]))
  })
  
  observeEvent(input$im_mz_map_condition, {
    updateSelectInput(inputId = "im_mz_map_replicate", choices = im_mz_map_replicates_available())
    
  })
  
  observeEvent(input$im_mz_map_actionbutton, {
    
    show("im_mz_map_plot_hider")
    toggle(id = "im_mz_map_plot", condition = T)
    
    output$im_mz_map_plot <- renderPlot({
      
      isolate({
        G = plot_mz_IM_map(dt = data_input(),
                           condition = input$im_mz_map_condition,
                           replicate = input$im_mz_map_replicate,
                           binwidth_mz = input$im_mz_map_binwidth_mz,
                           binwidth_im = input$im_mz_map_binwidth_im,
                           count_upper_limit = input$im_mz_map_upper_limit_counts,
                           mz_limits = input$im_mz_map_sliderMZ,
                           im_limits = input$im_mz_map_sliderIM
        )
      })
      
      G
      
    },width = 400,height = 400)
    })
  
  output$RT_distribution_plot<- renderPlot({
    
    
    G = plot_median_RT_2(dt = data_input(),
                         input$RT_distribution_checkbox)
    
    G
    
  })
  
  ## Corr plots
  
  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "CorrPlot_conditions", choices = unique(metadata_input()[,"R.Condition"]))
  })
  
  list_matrices_for_each_metric_2 <- reactiveValues()
  
  observeEvent(input$CorrPlot_actionbutton,{ #updates list of matrix tables with each metric
    
    req(input$file1)
    req(input$groups_id_csv)
    
    if(!input$CorrPlot_metric_to_plot %in% names(list_matrices_for_each_metric_2)){
      
      software_used <-  ifelse(
        is.null(input$radio_input_software), "null",
        switch(input$radio_input_software,
               "DIANN" = "DIANN",
               "timsDIANN" = "timsDIANN",
               "Spectronaut" = "Spectronaut"))
      
      isolate({
        list_matrices_for_each_metric_2 = create_list_matrices_for_each_metric_2(
          list_matrices_for_each_metric_2,
          dt = data_input(),
          metadata = metadata_input(),
          software = software_used,
          metric_to_plot = input$CorrPlot_metric_to_plot)
        })
      }
    })
  
  observeEvent(input$CorrPlot_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("Corr_Plot_hider")
    toggle(id = "Corr_Plot", condition = T)
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    output$Corr_Plot<- renderPlot({
      
      isolate({
        # G = Boxplot_OneMetric(list_matrices_for_each_metric,
        #                       metric_to_plot = input$BoxHist_metric_to_plot,
        #                       conditions = input$BoxHist_conditions,
        #                       plot_type = input$BoxHist_plot_type,
        #                       remove_outliers = input$BoxHist_checkbox)
        
        G = plot_correlations_fromList(list_matrices_for_each_metric_2,
                                       conditions = input$CorrPlot_conditions,
                                       metadata = metadata_input(),
                                       metric_to_plot = input$CorrPlot_metric_to_plot)
        
      })
      
      G
      
    })
  })
  
  #### PCA
  list_matrices_for_each_metric_pca <- reactiveValues()

  observeEvent(input$PCA_actionbutton,{ #updates list of matrix tables with each metric

    req(input$file1)
    req(input$groups_id_csv)

    if(!any(names(list_matrices_for_each_metric_pca) %in% c("PG.MS2Quantity"))){

      software_used <-  ifelse(
        is.null(input$radio_input_software), "null",
        switch(input$radio_input_software,
               "DIANN" = "DIANN",
               "timsDIANN" = "timsDIANN",
               "Spectronaut" = "Spectronaut"))

      isolate({

        list_matrices_for_each_metric_pca = create_list_matrices_for_each_metric_2(
          list_matrices_for_each_metric_pca,
          dt = data_input(),
          metadata = metadata_input(),
          software = software_used,
          metric_to_plot = "PG.MS2Quantity")


      })
    }
  })
  
  observeEvent(input$PCA_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("PCA_hider_pca")
    toggle(id = "PCA_Plot_pca", condition = T)
    
    
    output$PCA_Plot_pca<- renderPlotly({
      
      isolate({

        pca_data_2 <- create_PCA_data(list_matrices_for_each_metric_pca, metadata = metadata_input())
        
        
        G = plot_PCA(pca.data = pca_data_2,
                     plot_type = "pca",
                     metadata = metadata_input())
        
      })
      
      G
      
    })
  })
  
  observeEvent(input$PCA_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("PCA_hider_scree")
    toggle(id = "PCA_Plot_scree", condition = T)
    
    
    output$PCA_Plot_scree<- renderPlotly({
      
      isolate({
        
        pca_data_2 <- create_PCA_data(list_matrices_for_each_metric_pca, metadata = metadata_input())
        
        
        G = plot_PCA(pca.data = pca_data_2,
                     plot_type = "scree",
                     metadata = metadata_input())
        
      })
      
      G
      
    })
  })
  
  ## QC metrics plots
  
  observeEvent(input$groups_id_csv, {
    req(input$groups_id_csv)
    updateSelectInput(inputId = "BoxHist_conditions", choices = unique(metadata_input()[,"R.Condition"]))
  })
  
  observe({
    req(input$file1)
    req(input$groups_id_csv)
    
    columns_available = names(data_input())[names(data_input()) %in% c(
      "EG.Cscore", "EG.Qvalue", "FG.PrecMz","EG.PeakWidth", "EG.FWHM","PG.MS2Quantity","FG.MS2Quantity",
      "EG.ApexRT","EG.IonMobility", "FG.MS2RawQuantity","EG.Cscore","FG.PrecMz","EG.FWHM", "EG.PeakWidth")]
    
    updateSelectInput(inputId = "BoxHist_metric_to_plot", choices = columns_available)
  })
  
  
  list_matrices_for_each_metric <- reactiveValues()
  
  observeEvent(input$BoxHist_actionbutton,{#updates list of long format tables with each metric
    
    req(input$file1)
    req(input$groups_id_csv)
    
    if(!input$BoxHist_metric_to_plot %in% names(list_matrices_for_each_metric)){
      
      software_used <-  ifelse(
        is.null(input$radio_input_software), "null",
        switch(input$radio_input_software,
               "DIANN" = "DIANN",
               "timsDIANN" = "timsDIANN",
               "Spectronaut" = "Spectronaut"))
      
      isolate({
        list_matrices_for_each_metric = create_list_matrices_for_each_metric(
          list_matrices_for_each_metric,
          dt = data_input(),
          metadata = metadata_input(),
          software = software_used,
          metric_to_plot = input$BoxHist_metric_to_plot)
        })
      
    }
  })
  
  observeEvent(input$BoxHist_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("BoxHist_Plot_hider")
    toggle(id = "BoxHist_Plot", condition = T)
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    output$BoxHist_Plot<- renderPlot({
      
      isolate({
        G = Boxplot_OneMetric(list_matrices_for_each_metric,
                              metric_to_plot = input$BoxHist_metric_to_plot,
                              conditions = input$BoxHist_conditions,
                              plot_type = input$BoxHist_plot_type,
                              remove_outliers = input$BoxHist_checkbox)
        })
      
      G
      
    })
  })
  
  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "PeptideLength_conditions", choices = unique(metadata_input()[,"R.Condition"]))
    updateSelectInput(inputId = "ChargeStates_conditions", choices = unique(metadata_input()[,"R.Condition"]))
    
  })
  
  observeEvent(input$PeptideLength_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("PeptideLength_Plot_hider")
    toggle(id = "PeptideLength_Plot", condition = T)
    
    output$PeptideLength_Plot<- renderPlot({
      
      isolate({
        G = plot_peptide_length(dt = data_input(),
                               conditions = input$PeptideLength_conditions,
                               type = input$PeptideLength_plot_type)
      })
      
      G
      
    })
  })
  
  observeEvent(input$ChargeStates_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("ChargeStates_Plot_hider")
    toggle(id = "ChargeStates_Plot", condition = T)
    
    output$ChargeStates_Plot<- renderPlot({
      
      isolate({
        G = plot_charge_states(dt = data_input(),
                               conditions = input$ChargeStates_conditions)
      })
      
      G
      
    })
  })
  
  #### PTMS

  Modifications_found <- reactive({
    req(data_input())
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    mod = data.frame(modifications = find_modifications(data_input(), software = software_used))
    return(mod$modifications)
  })
  
  output$Modifications_found_table <- renderTable({
    req(data_input())
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    mod = data.frame(modifications = find_modifications(data_input(), software = software_used))
    return(mod)
  })
  
  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "ptm_counts_selectinput", choices = unique(Modifications_found()))
  })
  
  ptm_plotlist_counts <- reactive({
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    ptm_plotlist_counts <- PTM_counts_list_O_plots_filter(data_input(),
                                                          reorder =input$reassign_checkbox_input,
                                                          plot_type = input$PTMcounts_plot_type,
                                                          software= software_used)
    
    ptm_plotlist_counts
  })
  ptm_plotlist_enrichment <- reactive({
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    ptm_plotlist_enrich <- PTM_enrichment_yield_list_O_plots_filter(dt = data_input(),
                                                             counts = counts_data(),
                                                             software= software_used,
                                                             plot_type = input$PTMenrich_plot_type,
                                                             reorder = input$reassign_checkbox_input
                                                             )
    
    ptm_plotlist_enrich
  })
  ptm_plotlist_CVs <- reactive({
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    ptm_plotlist_CVs <- PTM_CVs_list_O_plots(CVs = CVs_data(),
                                             dt = data_input(),
                                             reorder = input$reassign_checkbox_input,
                                             software= software_used)
    
    ptm_plotlist_CVs
  })
  
  
  output$PTM_plot_Counts <- renderPlot({
    req(data_input())
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    PTM_plot_counts(ptm_plotlist_counts(),
                    input$ptm_counts_selectinput)
  })
  
  output$PTM_EnrichYield <- renderPlot({
    req(data_input())
    
    PTM_plot_enrichmentYield(ptm_plotlist_enrichment(),
                     input$ptm_counts_selectinput)
  })
  
  output$PTM_CVs <- renderPlot({
    req(data_input())
    
    PTM_plot_CVs(ptm_plotlist_CVs(),
                             input$ptm_counts_selectinput)
  })
  

  
  
  
  ###CVs
  observeEvent(input$groups_id_csv, {
    req(input$groups_id_csv)
    updateSelectInput(inputId = "CVs_conditions", choices = unique(metadata_input()[,"R.Condition"]))
  })
  
  observeEvent(input$CVs_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("CVs_Plot_hider")
    toggle(id = "CVs_Plot", condition = T)
    
    software_used <-  ifelse(
      is.null(input$radio_input_software), "null",
      switch(input$radio_input_software,
             "DIANN" = "DIANN",
             "timsDIANN" = "timsDIANN",
             "Spectronaut" = "Spectronaut"))
    
    CVs_level_2 <-  ifelse(
      is.null(input$CVs_level), "null",
      switch(input$CVs_level,
             "Protein" = "Protein",
             "Precursor" = "Precursor"))
    
    output$CVs<- renderPlot({
      
      isolate({
        
        if (CVs_level_2 == "Protein") {
          CVs_2 = CVs_proteins_data()
        } else {
          CVs_2 = CVs_data()
        }
        
        title = ifelse(CVs_level_2 == "Protein",
                       yes = "CVs(%) - Protein level",
                       no = "CVs(%) - Precursor level" )
        
        G = plot_CVs_reordered_2(CVs_2,
                                 metadata_input(),
                                 reorder = input$reassign_checkbox_input,
                                 conditions = input$CVs_conditions,
                                 plot_type = input$CVs_plot_type,
                                 remove_outliers = input$CVs_checkbox)+
          ggtitle(title)
        
        
        
      })
      
      G
      
    })
  })

  
  output$CVs_prots <- renderPlot({
    plot_CVs_reordered(CVs_proteins_data(), data_input(), input$reassign_checkbox_input)
  },
  height = 800,
  width = 800)
  
  
  
  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "CVs_cutoff_precursors_conditions", choices = unique(metadata_input()[,"R.Condition"]))
    updateSelectInput(inputId = "CVs_cutoff_proteins_conditions", choices =unique(metadata_input()[,"R.Condition"]))
  })
  
  observeEvent(input$CVs_cutoff_precursors_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("CVs_cutoff_precursors_Plot_hider")
    toggle(id = "CVs_cutoff_Precursors", condition = T)
    
    
    output$CVs_cutoff_Precursors<- renderPlot({
      
      isolate({
        
        G = plot_CVs_counts_Precursor(CVs_data(), data_input(),
                                      input$reassign_checkbox_input,
                                      input$CVs_cutoff_precursors_conditions,
                                      type = input$CVs_cutoff_precursors_type)
        
        
        
      })
      
      G
      
    })
  })
  
  observeEvent(input$CVs_cutoff_proteins_actionbutton, {
    req(input$file1)
    req(input$groups_id_csv)
    
    show("CVs_cutoff_proteins_Plot_hider")
    toggle(id = "CVs_cutoff_proteins", condition = T)
    
    
    output$CVs_cutoff_proteins<- renderPlot({
      
      isolate({
        
        G = plot_CVs_counts_Protein(CVs_proteins_data(), data_input(),
                                      input$reassign_checkbox_input,
                                      input$CVs_cutoff_proteins_conditions,
                                      type = input$CVs_cutoff_proteins_type)
      })
      
      G
      
    })
  })
  

  # output$CVs_cutoff_Proteins <- renderPlot({
  #   plot_CVs_counts_Protein(CVs_proteins_data(), data_input(),
  #                           input$reassign_checkbox_input,
  #                           input$CVs_cutoff_proteins_conditions)
  # })
  
  output$data_completeness <- renderPlot({
    plot_data_completeness(data_input())
  })
  
  output$intensity_boxplot_plot <- renderPlot({
    plot_intensity_boxplot(data_input(), input$reassign_checkbox_input)
  })
  ### QCs
  
  
  ### Dilution series
  
  normalizing_cond_reactive2 <- reactiveValues(data = NULL)
  
  observeEvent(input$dil_series_Conditions, {
    normalizing_cond_reactive2$normalizing_cond_text_input2 <- input$dil_series_Conditions
  })
  
  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "dil_series_Conditions", choices = unique(metadata_input()[,"R.Condition"]))
  }) 

  output$normalizing_cond_value2 <- renderText({
    paste0("normalizing to condition ", normalizing_cond_reactive2$normalizing_cond_text_input2) })
  

  output$calcurve_normalized2 <- renderPlot({
    cal_curve_norm_to_single_condition(dt = data_input(),
                                       normalizing_condition = normalizing_cond_reactive2$normalizing_cond_text_input2,
                                       reorder = input$reassign_checkbox_input)
  },
  width = function() input$plot_width,
  height = function() input$plot_height)
  
  ### PRM tools
  output$Download_list_of_Observed_Proteins_n_Peptides <- downloadHandler(
    
    filename = function() {
      paste("list_of_Observed_Proteins_n_Peptides-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      list_of_Observed_Proteins_n_Peptides = List_of_Observed_Proteins_n_Peptides(data_input())
      
      fwrite(list_of_Observed_Proteins_n_Peptides, file, sep = ",")
    })
  
  output$plot_counts_PeptideTargets <- renderPlot({
    plot_counts_Targets(dt = data_input(),
                        Targets = Targets_input(),
                        reorder = input$reassign_checkbox_input,
                        pepts_or_prots = "pepts")
  })
  
  observeEvent(input$Targets_input, {
    updateSelectInput(inputId = "Calcurve_ModPept_selection", choices = unique(Targets_input()[,"EG.ModifiedPeptide"]))
    
  })
  
  
  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "Calcurve_NormCondition_selection", choices = unique(metadata_input()[,"R.Condition"]))
  }) 
  output$cal_curve_norm_to_single_condition_singlePept <- renderPlot({
    cal_curve_norm_to_single_condition_singlePept(dt = data_input(),
                                                  normalizing_condition = input$Calcurve_NormCondition_selection,
                                                  reorder = input$reassign_checkbox_input,
                                                  PepSeq = input$Calcurve_ModPept_selection)
  })

  output$plotly_Peptide_Intensity_accross_runs <- renderPlotly({
    
    plotly_Peptide_Intensity_accross_runs(data_input(),Targets_input(),input$reassign_checkbox_input)
    
  })
  
  
  
  output$plot_HeatMap_Peptide_Targets <- renderPlotly({
    
    plot_HeatMap_Peptide_Targets(data_input(), Targets_input())
    
  })
  
  output$plot_HeatMap_Peptide_Targets_ratios_by_condition <- renderPlotly({
    
    plot_HeatMap_Peptide_Targets_ratios_by_condition(data_input(),
                                                    Targets_input(),
                                                    input$Calcurve_NormCondition_selection)
    
  })
  
  output$plot_counts_Protein_Targets <- renderPlot({
    plot_counts_Targets(dt = data_input(),
                        Targets = Targets_input(),
                        reorder = input$reassign_checkbox_input,
                        pepts_or_prots = "prots")
  })
  
  output$plotly_Protein_Intensity_accross_runs <- renderPlotly({
    
    plotly_Protein_Intensity_accross_runs(data_input(),Targets_input(),input$reassign_checkbox_input)
    
  })
  
  
  
  output$plot_HeatMap_Protein_Targets <- renderPlotly({
    
    plot_HeatMap_Protein_Targets(data_input(), Targets_input())
    
  })
  
  output$plot_HeatMap_Protein_Targets_ratios_by_condition <- renderPlotly({
    
    plot_HeatMap_Protein_Targets_ratios_by_condition(data_input(),
                                                     Targets_input(),
                                                     input$Calcurve_NormCondition_selection)
    
  })
  
  ###waterfall plot
  

  observeEvent(input$groups_id_csv, {
    updateSelectInput(inputId = "waterfall_condition", choices = unique(metadata_input()[,"R.Condition"]))
     # updateSelectInput(inputId = "waterfall_targets", choices = unique(data_input()[,"PG.ProteinGroups"]))
    
  })
  
  
  list_of_prot_targets <- reactive({
    req(input$Targets_input)
    Targets_input <- fread(input$Targets_input$datapath)
    unique(Targets_input$PG.ProteinGroups)
    
  })
  
  observeEvent(input$waterfall_actionbutton, {
    output$watefall_plot <- renderPlot({
      

      G = plot_and_table_waterfall(dt = data_input(),
                                   condition = input$waterfall_condition, 
                                   mean_or_median = input$waterfall_mean_or_median,
                                   quant_value_to_use = input$waterfall_quant_value_to_use,
                                   highlight_targets = input$waterfall_checkbox_input,
                                   # list_of_targetProteinGroup = input$waterfall_targets
                                   list_of_targetProteinGroup =list_of_prot_targets()
                                   )
      
      G
      
    },width = 400,height = 600)
  })

  
  
  output$format_prm_method <- renderTable({
    
    PRM_table <- format_diaData_to_PRMmethod(data_input(),
                                             Targets = Targets_input(),
                                             Isolation_Width = 3,
                                             RT_Range_seconds=300,
                                             IM_range = 0.07,
                                             experiment_name = "Targeted_peptides")
    
    return(PRM_table)
  })
  
  ## GCT tools
  
  GCT_file <- reactive({
    
    GCT_file <- create_gct_PG.Level(data_input(), metadata = metadata_input())
    
    GCT_file
  })
  
  output$download_gct_PGLevel <- downloadHandler(
    
    filename = "gct_PGLevel.gct",
    
    content = function(file) {
      cmapR::write_gct(ds = GCT_file(),
                       ofile= file,
                       precision = 3, appenddim = F, ver = 3)},
    contentType = "gct")
}

# Create Shiny app ----
shinyApp(ui, server)