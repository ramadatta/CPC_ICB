#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

server <- function (input, output, session){
  
  source('/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R',local = TRUE)  
  #source('ui.R')
  
  test.visn <- list(nodes = NULL, edges = NULL)
  
  observeEvent(input$df_test, {
    selected_option <- input$df_test
    

    if (!is.null(selected_option)) {
      # if (selected_option == "CT_uninfected_CT_3d_vs_KO_uninfected_KO_3d") {
      #   # Process CT_uninfected_CT_3d and KO_uninfected_KO_3d files
      #   ct_file <- "2021_01_24_CT_uninfected_CT_3d-CT_qual.csv"
      #   ko_file <- "2021_01_24_KO_uninfected_KO_3d-KO_qual.csv"
      #   day <- "3d"
      #   
      #   # Load the '/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R' script
      #   source("/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R")
      #   # Call the function with the file names
      #   result <- process_files(ct_file, ko_file, day)
      #   test.visn$nodes <- result$nodes
      #   test.visn$edges <- result$edges
      # 
      # } else
      
    if (selected_option == "CT_uninfected_CT_6d_vs_KO_uninfected_KO_6d") {
      # Process CT_uninfected_CT_6d and KO_uninfected_KO_6d files
      ct_file <- "2021_01_24_CT_uninfected_CT_6d-CT_qual.csv"
      ko_file <- "2021_01_24_KO_uninfected_KO_6d-KO_qual.csv"
      day <- "6d"
      
      # Load the '/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R' script
      source("/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R")
      
      # Call the function with the file names
      result <- process_files(ct_file, ko_file, day)
      test.visn$nodes <- result$nodes
      test.visn$edges <- result$edges
      targets_cells_top5ligands_df <- result$targets_cells_top5ligands_df
      
    } else if (selected_option == "CT_uninfected_CT_15d_vs_KO_uninfected_KO_15d") {
      # Process CT_uninfected_CT_15d and KO_uninfected_KO_15d files
      ct_file <- "2021_01_24_CT_uninfected_CT_15d-CT_qual.csv"
      ko_file <- "2021_01_24_KO_uninfected_KO_15d-KO_qual.csv"
      day <- "15d"
      
      # Load the '/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R' script
      source("/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R")
      
      # Call the function with the file names
      result <- process_files(ct_file, ko_file, day)
      test.visn$nodes <- result$nodes
      test.visn$edges <- result$edges
      targets_cells_top5ligands_df <- result$targets_cells_top5ligands_df

      
    } else if (selected_option == "CT_uninfected_CT_28d_vs_KO_uninfected_KO_28d") {
      # Process CT_uninfected_CT_28d and KO_uninfected_KO_28d files
      ct_file <- "2021_01_24_CT_uninfected_CT_28d-CT_qual.csv"
      ko_file <- "2021_01_24_KO_uninfected_KO_28d-KO_qual.csv"
      day <- "28d"
      
      # Load the '/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R' script
      source("/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R")
      
      # Call the function with the file names
      result <- process_files(ct_file, ko_file, day)
      test.visn$nodes <- result$nodes
      test.visn$edges <- result$edges
      targets_cells_top5ligands_df <- result$targets_cells_top5ligands_df
    
      
    } else if (selected_option == "CT_uninfected_CT_45d_vs_KO_uninfected_KO_45d") {
      # Process CT_uninfected_CT_45d and KO_uninfected_KO_45d files
      ct_file <- "2021_01_24_CT_uninfected_CT_45d-CT_qual.csv"
      ko_file <- "2021_01_24_KO_uninfected_KO_45d-KO_qual.csv"
      day <- "45d"
      
      # Load the '/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R' script
      source("/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R")
      
      # Call the function with the file names
      result <- process_files(ct_file, ko_file, day)
      test.visn$nodes <- result$nodes
      test.visn$edges <- result$edges
      targets_cells_top5ligands_df <- result$targets_cells_top5ligands_df
      
    } else if (selected_option == "CT_uninfected_CT_100d_vs_KO_uninfected_KO_100d") {
      
      # Process CT_uninfected_CT_100d and KO_uninfected_KO_100d files
      ct_file <- "2021_01_24_CT_uninfected_CT_100d-CT_qual.csv"
      ko_file <- "2021_01_24_KO_uninfected_KO_100d-KO_qual.csv"
      day <- "100d"
      
      # Load the '/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R' script
      source("/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny/testShiny/MHV68/basic.R")
      
      # Call the function with the file names
      result <- process_files(ct_file, ko_file, day)
      test.visn$nodes <- result$nodes
      test.visn$edges <- result$edges
      targets_cells_top5ligands_df <- result$targets_cells_top5ligands_df
      
    }
    
    
  }
    })
  
  
  
  output$network <- visNetwork::renderVisNetwork({
    visNetwork::visNetwork(nodes = test.visn$nodes,test.visn$edges %>% 
                             visNetwork::visOptions(highlightNearest=TRUE, 
                 nodesIdSelection = TRUE) %>%
      #allow for long click to select additional nodes
        visNetwork::visInteraction(multiselect = TRUE) %>%
        visNetwork::visIgraphLayout(layout = "layout_in_circle") 
      %>%
        visNetwork::visEdges(arrows = 'to', smooth = TRUE, color = list(highlight = "#28A99E", hover = "#28A99E")) %>%
        visNetwork::visNodes(size = 20, color = list(background = "#28A99E", border = "black", highlight = "yellow",
                                       hover = "yellow"),
               shapeProperties = list(useImageSize=FALSE,
                                      interpolation=FALSE)) %>%
        visNetwork::visOptions(highlightNearest = list(enabled = F, degree = 0,
                                         hover = T,
                                         hideColor = 'rgba(200,200,200,0.2)',
                                         #hideColor = 'rgba(200,200,200,0)',
                                         labelOnly=TRUE),
                 nodesIdSelection=TRUE) %>%
        visNetwork::visInteraction(hover=TRUE, zoomView = TRUE,
                     navigationButtons = TRUE,
                     tooltipStyle = 'position: fixed;visibility:hidden;padding: 5px;
                font-family: sans-serif;font-size:12px;
                font-color:#000000;background-color: #e3fafa;
                -moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;
                 border: 0px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.5);
                 max-width:300px;overflow-wrap: normal') %>%
      addFontAwesome() %>%
        visNetwork::visLayout(randomSeed = 02143, improvedLayout=TRUE) %>%

      #Use visEvents to turn set input$current_node_selection to list of selected nodes
        visNetwork::visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_selection', nodes.nodes);
                ;}")
      )
      # 
  })
  
  #render data table restricted to selected nodes
  # output$tbl <- renderDT(
  #   test.visn$edges %>% 
  #     filter((to %in% input$current_node_selection)|(from %in% input$current_node_selection)),
  #   options = list(lengthChange = FALSE)
  # )
  
  # render data table restricted to selected nodes
  # output$tbl <- renderDT(
  #   targets_cells_top5ligands_df %>% 
  #     filter(Cell %in% input$current_node_selection),
  #   options = list(lengthChange = FALSE)
  # )
  
  # download button
  output$download <- downloadHandler(
    filename = function() {
      paste("selected_nodes_data", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv", sep = "")
    },
    content = function(file) {
      selected_nodes_data <- targets_cells_top5ligands_df %>% 
        filter(Cell %in% input$current_node_selection)
      write.csv(selected_nodes_data, file, row.names = FALSE)
    }
  )
  
 
  
  
}

#shinyApp(ui, server)