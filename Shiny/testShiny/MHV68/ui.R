#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(igraph)
library(visNetwork)
library(dplyr)
library(shiny)
library(shinythemes)
library(DT)
library(plotrix)
library(ggsci)
library(tidyverse)


ui <- fluidPage(
  
  # Generate Title Panel at the top of the app
  titlePanel("MHV68 Project Network Visualization Dashboard"),
  
  # Add a select input for time points
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "time_point",
        "Select Time Point:",
        choices = c("3d", "6d", "15d", "28d", "45d", "100d"),
        selected = "6d"
      )
    ),
    
    mainPanel(
      fluidRow(
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
    )
  )
) #end of fluidPage

# 
# ui <- shinyUI(fluidPage(
#   
#   # Generate Title Panel at the top of the app
#   titlePanel("MHV68 Project Network Visualization Dashboard"),
#   # fluidRow(
#   #   selectInput(
#   #     inputId = "df_test",
#   #     label = "Select a Condition",
#   #     choices = c(
#   #      # "Difference CT_uninfected_CT_3d, KO_uninfected_KO_3d" = "CT_uninfected_CT_3d_vs_KO_uninfected_KO_3d",
#   #       "Difference CT_uninfected_CT_6d, KO_uninfected_KO_6d" = "CT_uninfected_CT_6d_vs_KO_uninfected_KO_6d",
#   #       "Difference CT_uninfected_CT_15d, KO_uninfected_KO_15d" = "CT_uninfected_CT_15d_vs_KO_uninfected_KO_15d",
#   #       "Difference CT_uninfected_CT_28d, KO_uninfected_KO_28d" = "CT_uninfected_CT_28d_vs_KO_uninfected_KO_28d",
#   #       "Difference CT_uninfected_CT_45d, KO_uninfected_KO_45d" = "CT_uninfected_CT_45d_vs_KO_uninfected_KO_45d",
#   #       "Difference CT_uninfected_CT_100d, KO_uninfected_KO_100d" = "CT_uninfected_CT_100d_vs_KO_uninfected_KO_100d"
#   #     ),
#   #     selected = 'Difference CT_uninfected_CT_100d, KO_uninfected_KO_100d',
#   #     width = "50%"
#   #   ),
#   #  # DT::dataTableOutput("test_table")
#   # ),
#   
#   fluidRow(
#     br(),
#     br(),
#     br(),
#     br(),
#     column(width = 5, DTOutput('tbl')),
#     column(width = 7, visNetworkOutput("network")) #note that column widths in a fluidRow should sum to 12
#   ),
#   
#   fluidRow(
#     column(width = 6, offset = 1, 
#            downloadButton("download", "Download Selected Nodes Data")
#            #downloadButton("download", paste0("Download Selected Nodes Data_", format(Sys.time(), "%Y%m%d%H%M%S")))
#            
#     ),
#     fluidRow(
#     column(width = 7, offset = 6,
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            br(),
#            "Select multiple nodes by clicking and holding for a second. The datatable on the left will update accordingly."
#     )
#   )
#   
# ) #end of fluidPage
# )
# )