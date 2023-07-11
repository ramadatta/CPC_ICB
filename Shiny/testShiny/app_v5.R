# Function to process files
process_files <- function(ct_file, ko_file, day) {
  # Read the input file
  mat1 <- read.csv(ct_file, row.names = 1, check.names = FALSE, header = TRUE)
  mat1 <- as.matrix(mat1)
  mat2 <- read.csv(ko_file, row.names = 1, check.names = FALSE, header = TRUE)
  mat2 <- as.matrix(mat2)
  
  ## do not consider autocrine feedback
  diag(mat1) = 0
  diag(mat2) = 0
  
  # Difference of the conditions
  mat <- mat2 - mat1
  #View(mat)
  
  
  ## Custom thresholding. We don't want to show all of the arrows.
  thresh = 0.05
  #mat[mat < thresh] = 0
  mat[abs(mat) < thresh] = 0
  
  #reorder cols and rows to make plot visually more appealing
  matrix <- mat[c('AT2',
                  'Airway Epithelium',
                  'Alveolar Macrophages',
                  'Capillary cells',
                  'Dendritic cells',
                  'Fibroblasts and Protector cells',
                  'Granulocytes',
                  'Innate Lymphoid cells',
                  'Interstitial Macrophages',
                  'non classical Monocytes',
                  'classical Monocytes',
                  'Monocytes and M2 Macrophages',
                  'NKT cells',
                  'Plasma B cells',
                  'T cells',
                  'Vascular Endothelial cells'),
                c('AT2',
                  'Airway Epithelium',
                  'Alveolar Macrophages',
                  'Capillary cells',
                  'Dendritic cells',
                  'Fibroblasts and Protector cells',
                  'Granulocytes',
                  'Innate Lymphoid cells',
                  'Interstitial Macrophages',
                  'non classical Monocytes',
                  'classical Monocytes',
                  'Monocytes and M2 Macrophages',
                  'NKT cells',
                  'Plasma B cells',
                  'T cells',
                  'Vascular Endothelial cells')]
  
  ## make igraph object
  test.gr <- graph_from_adjacency_matrix(matrix, mode="directed", weighted=T)
  #View(test.gr)
  #class(test.gr)
  
  
  ## convert to VisNetwork-list
  test.visn <- toVisNetworkData(test.gr)
  
  ## highlight differences
  #test.visn$edges$weight <- test.visn$edges$weight * 20
  
  ## copy column "weight" to new column "value" in list "edges"
  test.visn$edges$value <- test.visn$edges$weight
  
  # make sure some values are negative
  #head(test.visn$edges)
  
  
  
  ## add colors to edge
  
  # Initialize color vector with NA
  edge_colors <- rep(NA, length(test.visn$edges$value))
  
  # Assign colors for positive, negative, and zero weights
  positive_indices <- which(test.visn$edges$value > 0)
  negative_indices <- which(test.visn$edges$value < 0)
  
  edge_colors[positive_indices] <- color.scale(test.visn$edges$value[positive_indices], extremes = c(pal_material('red')(10)[1], pal_material('red')(10)[10]))
  edge_colors[negative_indices] <- color.scale(test.visn$edges$value[negative_indices], extremes = c(pal_material('light-blue')(10)[10], pal_material('light-blue')(10)[1]))
  
  test.visn$edges$color <- edge_colors
  
  ## So far so good, can almost replicate niklas figure
  
  # feature 1 ---------------------------------------------------------------
  
  # when a node is selected - show top 5 target genes and ligands from KO (to confirm?) in a text box- we can tinker below code to show text 
  
  # feature 3b with real ligands --------------------------------------------
  
  # Shell command
  
  # for d in nichenet_ligand_act*.csv; do echo -en "$d\t"; head -6 $d | tail -n+2 | awk -F, '{print $2}' | sed 's/"//g' | tr '\n' ',' | sed 's/,$//g'; echo ""; done | sed -e 's/nichenet_ligand_act_//g' -e 's/_/ /g' -e 's/.csv//g' >top5_ligands_KO.txt
  # Run makeLigand_TargetGene_List_v2.R
  
  # Create dataframe *************Need modification***********
  targets_cells_top5ligands_df <- data.table::fread("Top5_Ligands_TargetGenes_6d.txt",sep="\t", header = TRUE) 
  
  # Print the dataframe
  #head(targets_cells_top5ligands_df)
  
  cells_top5ligands_df <- targets_cells_top5ligands_df %>% 
    group_by(Cell) %>% 
    summarize(Ligand = paste(Ligand, collapse = ", "))
  
  # feature 5 ---------------------------------------------------------------
  
  # when a edge is selected - show top 5 target genes and ligands from KO (to confirm?) for both nodes in a text box - need to explore
  # I can get the information displayed when hovered over the edge - but it would be nice if the nodes are selected and two boxes appear 
  # USeful links
  
  #https://stackoverflow.com/questions/39758792/display-information-in-network-edges-using-visnetwork?rq=4
  #https://stackoverflow.com/questions/64669941/adding-additional-information-to-a-visnetwork
  #https://stackoverflow.com/questions/39905061/get-node-and-edge-data-from-visnetwork-graph?rq=3
  #https://stackoverflow.com/questions/52519042/update-a-datatable-based-on-visnetwork-drop-down-node-selection-in-a-shiny-app
  
  # head(test.visn$edges)
  # head(targets_cells_top5ligands_df,10)
  # head(cells_top5ligands_df)
  
  collapsed_df <- targets_cells_top5ligands_df %>%
    group_by(Cell) %>%
    reframe(EdgeInfo = paste("<b>Cell:</b>", Cell, "<br><br><b>Ligand: Target</b><br><br>", 
                             paste0("--", Ligand, ": ", TargetGenes, collapse = "<br>"), sep = "\n")) %>% unique() %>% as.data.frame()
  
  # class(collapsed_df)
  # 
  # head(collapsed_df)
  
  merged_df <- left_join(test.visn$edges, collapsed_df, by = c("from" = "Cell"))
  
  
  merged_df <- left_join(merged_df, collapsed_df, by = c("to" = "Cell")) 
  
  #tibble(merged_df)
  
  merged_df <- merged_df %>% 
    unite(title, EdgeInfo.x, EdgeInfo.y, sep = "<br><br>") 
  
  #tibble(merged_df)
  
  # Assign the updated "title" column back to test.visn$edges
  test.visn$edges$title <- merged_df$title
  tibble(test.visn$edges)
  
  return(list(nodes = test.visn$nodes, edges = test.visn$edges, targets_cells_top5ligands_df = targets_cells_top5ligands_df))
}

ui <- fluidPage(
  
  # Generate Title Panel at the top of the app
  titlePanel("MHV68 Project Network Visualization Dashboard"),
  
  # Add a select input for time points
  # sidebarLayout(
  #   sidebarPanel(
  #     # selectInput(
  #     #   "time_point",
  #     #   "Select Time Point:",
  #     #   choices = c("3d", "6d", "15d", "28d", "45d", "100d"),
  #     #   selected = "6d"
  #     # )
  #         selectInput(
  #           inputId = "time_point",
  #           label = "Select a Condition",
  #           choices = c(
  #             "Difference CT_uninfected_CT_3d, KO_uninfected_KO_3d" = "3d",
  #             "Difference CT_uninfected_CT_6d, KO_uninfected_KO_6d" = "6d",
  #             "Difference CT_uninfected_CT_15d, KO_uninfected_KO_15d" = "15d",
  #             "Difference CT_uninfected_CT_28d, KO_uninfected_KO_28d" = "28d",
  #             "Difference CT_uninfected_CT_45d, KO_uninfected_KO_45d" = "45d",
  #             "Difference CT_uninfected_CT_100d, KO_uninfected_KO_100d" = "100d"),
  #           selected = '100d')
  #   ),
    
  #  mainPanel(
      
      fluidRow(         selectInput(
        inputId = "time_point",
        label = "Select a Condition",
        choices = c(
          "Difference CT_uninfected_CT_3d, KO_uninfected_KO_3d" = "3d",
          "Difference CT_uninfected_CT_6d, KO_uninfected_KO_6d" = "6d",
          "Difference CT_uninfected_CT_15d, KO_uninfected_KO_15d" = "15d",
          "Difference CT_uninfected_CT_28d, KO_uninfected_KO_28d" = "28d",
          "Difference CT_uninfected_CT_45d, KO_uninfected_KO_45d" = "45d",
          "Difference CT_uninfected_CT_100d, KO_uninfected_KO_100d" = "100d"),
        selected = '100d', width = "50%")),
      fluidRow(
            br(),
            br(),
            br(),
            br(),
        column(width = 5, DTOutput('tbl')),
        column(width = 7, visNetworkOutput("network"))
      ),
      
      fluidRow(
        column(width = 6, offset = 1, 
               downloadButton("download", "Download Selected Nodes Data")
        )
      ),
      
      fluidRow(
        column(width = 7, offset = 6,
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               "Select multiple nodes by clicking and holding for a second. The datatable on the left will update accordingly."
        )
      )
   # ) #main panel
 # ) #sidebarLayout
) #end of fluidPage

server <- function (input, output, session){
  
  setwd("/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny")
  
  
  # Reactive function to process the selected files
  processed_files <- reactive({
    selected_day_option <- input$time_point
    
    if (selected_day_option %in% c("3d", "6d", "15d", "28d", "45d", "100d")) {
      # Construct the file names based on the selected time point
      ct_file <- paste0("2021_01_24_CT_uninfected_CT_", selected_day_option, "-CT_qual.csv")
      ko_file <- paste0("2021_01_24_KO_uninfected_KO_", selected_day_option, "-KO_qual.csv")
      
      # Call the function to process the files
      process_files(ct_file, ko_file, selected_day_option)
      

    }
  })
  
  output$network <- renderVisNetwork({
    
    if (!is.null(processed_files()$nodes)) {
      visNetwork(nodes = processed_files()$nodes, edges = processed_files()$edges, height = "600px") %>%
      visOptions(highlightNearest=TRUE, 
                 nodesIdSelection = TRUE) %>%
      #allow for long click to select additional nodes
      visInteraction(multiselect = TRUE) %>%
      visIgraphLayout(layout = "layout_in_circle") %>% 
      visEdges(arrows = 'to', smooth = TRUE, color = list(highlight = "#28A99E", hover = "#28A99E")) %>%
      visNodes(size = 20, color = list(background = "#28A99E", border = "black", highlight = "yellow",
                                       hover = "yellow"),
               shapeProperties = list(useImageSize=FALSE, 
                                      interpolation=FALSE)) %>%
      visOptions(highlightNearest = list(enabled = F, degree = 0, 
                                         hover = T, 
                                         hideColor = 'rgba(200,200,200,0.2)',
                                         #hideColor = 'rgba(200,200,200,0)',
                                         labelOnly=TRUE),
                 nodesIdSelection=TRUE) %>%
      visInteraction(hover=TRUE, zoomView = TRUE,
                     navigationButtons = TRUE,
                     tooltipStyle = 'position: fixed;visibility:hidden;padding: 5px;
                font-family: sans-serif;font-size:12px;
                font-color:#000000;background-color: #e3fafa;
                -moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;
                 border: 0px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.5);
                 max-width:300px;overflow-wrap: normal') %>%
      addFontAwesome() %>%
      visLayout(randomSeed = 02143, improvedLayout=TRUE) %>% 
      
      #Use visEvents to turn set input$current_node_selection to list of selected nodes
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_selection', nodes.nodes);
                ;}")
    
  }
  })
  
  # # render data table restricted to selected nodes
  # output$tbl <- renderDT(
  #   processed_files()$targets_cells_top5ligands_df %>% 
  #     filter(Cell %in% input$current_node_selection),
  #   options = list(lengthChange = FALSE)
  # )
  # 
  # # download button
  # output$download <- downloadHandler(
  #   filename = function() {
  #     paste("selected_nodes_data", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv", sep = "")
  #   },
  #   content = function(file) {
  #     selected_nodes_data <- processed_files()$targets_cells_top5ligands_df %>% 
  #       filter(Cell %in% input$current_node_selection)
  #     write.csv(selected_nodes_data, file, row.names = FALSE)
  #   }
  # )
  
  # render data table restricted to selected nodes
  output$tbl <- renderDT({
    selected_nodes <- input$current_node_selection
    if (!is.null(selected_nodes)) {
      filtered_df <- processed_files()$targets_cells_top5ligands_df %>% filter(Cell %in% selected_nodes)
      print(filtered_df)
      datatable(filtered_df, options = list(lengthChange = FALSE))
    } else {
      datatable(NULL)
    }
  })
  
  # download button
  output$download <- downloadHandler(
    filename = function() {
      paste("selected_nodes_data", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv", sep = "")
    },
    content = function(file) {
      selected_nodes_data <- processed_files()$targets_cells_top5ligands_df %>% filter(Cell %in% input$current_node_selection)
      write.csv(selected_nodes_data, file, row.names = FALSE)
    }
  )
  
  
  
  
}


shinyApp(ui, server)