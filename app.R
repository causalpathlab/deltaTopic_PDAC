library(shiny, quietly = TRUE)
library(ggpubr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(Matrix, quietly = TRUE)
library(data.table, quietly = TRUE)
library(magrittr, quietly = TRUE)
SaveFolderPath <- "models/TotalDeltaETM_allgenes_ep1000_nlv32_bs512_combinebyadd_lr0.01_train_size1"
my_gene_list_df <- fread(paste0(SaveFolderPath, paste0("/","delta","_topK_genes.csv.gz"))) %>% as.data.frame()

K <- 10
my_gene_list <- c()
for(j in seq_len(ncol(my_gene_list_df))){
    my_gene_list <- c(my_gene_list, as.character(my_gene_list_df[1:K, j]))
}
my_gene_list <- as.factor(my_gene_list) %>% unique()

topics_df <- fread(paste0(SaveFolderPath, "/topics.csv.gz"))
topics_df <- topics_df[, -1]
topics_df$topics <- colnames(topics_df)[apply(topics_df, 1, which.max)]

U <- readRDS("figures/fig1b/data_top_u.rds")
S <- readRDS("figures/fig1b/data_top_s.rds")

ui <- pageWithSidebar(
  headerPanel('Too tired to think of a title'),
  sidebarPanel(
    checkboxGroupInput("topics", "topics to show:",
                        unique(topics_df$topics),
                        selected = unique(topics_df$topics)[1:5]),
    selectInput('my_gene_1', 'Gene 1', my_gene_list, selected = 'B2M'),
    selectInput('my_gene_2', 'Gene 2', my_gene_list, selected = 'MALAT1'),
    selectInput('my_gene_3', 'Gene 3', my_gene_list, selected = 'RPL10')
  ),
  fluidRow(
    column(width = 7,
    tableOutput('table1')
    ),
    column(width = 7,
    plotOutput('cell_count_per_topic')
    ),
    column(width = 12,
    plotOutput('plot1')
    ),
    column(width = 12,
    plotOutput('plot2')
    ),
    column(width = 12,
    plotOutput('plot3')
    )
  )
)

server <-  function(input, output, session) {

  df_1 <- reactive({
    df <- data.frame(S = S[, input$my_gene_1], U = U[, input$my_gene_1])
    df <- cbind(df, topics_df)
    df <- df[df$topics %in% input$topics, ]
    df
  })

  df_2 <- reactive({
    df <- data.frame(S = S[, input$my_gene_2], U = U[, input$my_gene_2])
    df <- cbind(df, topics_df)
    df <- df[df$topics %in% input$topics, ]
    df
  })

  df_3 <- reactive({
    df <- data.frame(S = S[, input$my_gene_3], U = U[, input$my_gene_3])
    df <- cbind(df, topics_df)
    df <- df[df$topics %in% input$topics, ]
    df
  })

  output$table1 <- renderTable({
    head(my_gene_list_df[, unique(topics_df$topics)], 5)
  })

  output$cell_count_per_topic <- renderPlot({
    ggplot(topics_df, aes(x = topics)) +
    geom_bar()+
    geom_text(stat='count', aes(label=..count..), vjust=-1)
  })

  output$plot1 <- renderPlot({
    df <- df_1()
    p1 <- ggplot(df, aes(x = S, y = U,color = topics)) +
      geom_point(alpha = 0.1) + 
      facet_grid(~topics,scale = "fixed") +
      geom_abline(intercept = 0, slope = 1, colour = "red")+ 
      ggtitle(input$my_gene_1) 
    p2 <- ggplot(df, aes(x = S, y = U, color = topics)) +
      geom_point(alpha = 0.5) +
      #geom_jitter() + 
      geom_abline(intercept = 0, slope = 1, colour = "red") 
      #ggtitle(input$my_gene)
      ggarrange(p1, p2, 
            labels = c("A", "B"),
            ncol = 1, nrow = 2)
  })

  output$plot2 <- renderPlot({
    df <- df_2()
    p1 <- ggplot(df, aes(x = S, y = U,color = topics)) +
      geom_point(alpha = 0.1) + 
      facet_grid(~topics,scale = "fixed") +
      geom_abline(intercept = 0, slope = 1, colour = "red") + 
      ggtitle(input$my_gene_2)
    p2 <- ggplot(df, aes(x = S, y = U, color = topics)) +
      geom_point(alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, colour = "red") 
      #ggtitle(input$my_gene)
      ggarrange(p1, p2, 
              labels = c("A", "B"),
              ncol = 1, nrow = 2)
  })

  output$plot3 <- renderPlot({
    df <- df_3()
    p1 <- ggplot(df, aes(x = S, y = U,color = topics)) +
      geom_point(alpha = 0.1) + 
      facet_grid(~topics,scale = "fixed") +
      geom_abline(intercept = 0, slope = 1, colour = "red") + 
      ggtitle(input$my_gene_3)
    p2 <- ggplot(df, aes(x = S, y = U, color = topics)) +
      geom_point(alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, colour = "red") 
      #ggtitle(input$my_gene)
      ggarrange(p1, p2, 
              labels = c("A", "B"),
              ncol = 1, nrow = 2)
  })
}
shinyApp(ui, server)