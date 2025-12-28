#!/usr/bin/env Rscript
# Author: C. Ryan Miller
# Created: 2025-12-28 02:20 CST
# Commit: added as source_shiny_app.R (renamed from kinome_shiny_app.R)

library(shiny)
library(DT)
library(data.table)

# Load data (adjust filename if your pipeline output differs)
dt <- fread("human_kinome_with_group_family_GOlast_custom2.csv")

ui <- fluidPage(
  titlePanel("Human Source Annotation Explorer"),
  sidebarLayout(
    sidebarPanel(
      helpText("Filter, search, and explore source (kinase) annotations."),
      selectInput("group_filter", "Group:", choices = c("All", sort(unique(dt$Group))), selected = "All"),
      selectInput("family_filter", "Family:", choices = c("All", sort(unique(dt$Family))), selected = "All"),
      downloadButton("downloadData", "Download Filtered Table")
    ),
    mainPanel(
      DTOutput("sourceTable")
    )
  )
)

server <- function(input, output, session) {
  filtered_dt <- reactive({
    d <- dt
    if (!is.null(input$group_filter) && input$group_filter != "All") {
      d <- d[d$Group == input$group_filter, ]
    }
    if (!is.null(input$family_filter) && input$family_filter != "All") {
      d <- d[d$Family == input$family_filter, ]
    }
    d
  })

  output$sourceTable <- renderDT({
    datatable(
      filtered_dt(),
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE),
      extensions = c('Buttons')
    )
  }, server = TRUE)

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("filtered_source_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_dt(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
