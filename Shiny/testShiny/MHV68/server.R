server <- function (input, output, session){
  
  setwd("/Users/ramadatta/Analysis/Herberts_Lab/2_MHV68_Mice_Analysis/from_Lukas_HD/Datta_Redo/NN/run_3_uninfected_same_condition_results/Shiny")
  
  
  # Reactive function to process the selected files
  processed_files <- reactive({
    selected_option <- input$time_point
    
    if (selected_option %in% c("3d", "6d", "15d", "28d", "45d", "100d")) {
      # Construct the file names based on the selected time point
      ct_file <- paste0("2021_01_24_CT_uninfected_CT_", selected_option, "-CT_qual.csv")
      ko_file <- paste0("2021_01_24_KO_uninfected_KO_", selected_option, "-KO_qual.csv")
      
      # Call the function to process the files
      process_files(ct_file, ko_file, selected_option)
    }
  })
  
  output$network <- renderVisNetwork({
    #visNetwork(nodes = test.visn$nodes,test.visn$edges) %>% 
    visNetwork(nodes = processed_files()$nodes, edges = processed_files()$edges) %>%
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
    
  })
  
}