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