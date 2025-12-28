# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

# kinome_shiny_app.R
# Shiny app for interactive viewing of human kinase annotations

library(shiny)
library(DT)
library(data.table)

# Load data
dt <- fread("human_kinome_with_group_family_GOlast_custom2.csv")

ui <- fluidPage(
  titlePanel("Human Kinome Annotation Explorer"),
  sidebarLayout(
    sidebarPanel(
      helpText("Filter, search, and explore kinase annotations. Use the table below to sort, filter, and download results."),
      selectInput("group_filter", "Group:", choices = c("All", sort(unique(dt$Group))), selected = "All"),
      selectInput("family_filter", "Family:", choices = c("All", sort(unique(dt$Family))), selected = "All"),
      downloadButton("downloadData", "Download Filtered Table")
    ),
    mainPanel(
      DTOutput("kinomeTable")
    )
  )
)

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

  output$kinomeTable <- renderDT({
    datatable(
      filtered_dt(),
      filter = "top",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      ),
      extensions = c('Buttons')
    )
  }, server = TRUE)

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("filtered_kinome_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_dt(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
