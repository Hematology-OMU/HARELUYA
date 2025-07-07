library(shiny)
library(randomForestSRC)
library(survival)
library(RColorBrewer)

#
# const
#

APP_TITLE <- "Calculator of the personalized OS curves for 10 allo-HSCT procedures using Random Survival Forest"
APP_NAME <- "HARELUYA"
CHART_DESC <- c(
  "Allo-HSCT procedures classify into the following ten kinds.",
  "R-MRD; RIC/MRD/non-PTCY, M-MRD; MAC/MRD/non-PTCY,",
  "R-MUD; RIC/MUD/non-PTCY, M-MUD; MAC/MUD/non-PTCY,",
  "R-UCB; RIC/UCB/non-PTCY, M-UCB; MAC/UCB/non-PTCY,",
  "R-Haplo-CY; RIC/Haplo/PTCY, M-Haplo-CY; MAC/Haplo/PTCY,",
  "R-Haplo-nCY; RIC/Haplo/non-PTCY, M-Haplo-nCY; MAC/Haplo/non-PTCY.",
  paste(
    "Abbreviation: Haplo; HLA haploidentical related donor,",
    "MAC; myeloablative conditioning, MRD; HLA matched related donor,",
    "MUD; HLA matched unrelated donor, PTCY; post-transplant cyclophosphamide,",
    "RIC; reduced-intensity conditioning, UCB; unrelated cord blood."
  )
)

Tx_pattern.names <- c(
  'R-MRD',
  'R-MUD',
  'R-UCB',
  'R-Haplo-CY',
  'R-Haplo-nonCY',
  'M-MRD',
  'M-MUD',
  'M-UCB',
  'M-Haplo-CY',
  'M-Haplo-nonCY'
)

INIT_DATA_PATH <- "www/data/data.csv"

#
# function
#

updateInitData <- function (preds) {
  df <- data.frame(preds[[1]]$time.interest)
  df <- Reduce(function (ldf, pred) {
    return (cbind(ldf, data.frame(c(pred$survival))))
  }, preds, df)
  names(df) <- c("x", Tx_pattern.names)
  write.csv(df, INIT_DATA_PATH, row.names=F)
}

#
# global
#

modelRFSRC <- readRDS("./modelRFSRC_TotalCase.obj")

ui <- navbarPage(
  APP_TITLE,
  tabPanel(
    "Component 1",
    tagList(
      singleton(
        tags$head(
          tags$link(href="https://cdnjs.cloudflare.com/ajax/libs/c3/0.7.20/c3.css", rel="stylesheet", type="text/css"),
          tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/d3/5.16.0/d3.min.js"),
          tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/c3/0.7.20/c3.min.js"),
          tags$link(href="css/shiny_TRUMP.css", rel="stylesheet", type="text/css"),
          tags$script(src="js/shiny_TRUMP.js"),
          tags$link(href="css/disclamer.css", rel="stylesheet", type="text/css"),
          tags$script(src="js/disclamer.js")
        )
      )
    ),
    div(class="app-name", APP_NAME),
    sidebarPanel(
      class = "input-panel",
      actionButton("go","Calculation"),
      sliderInput(".Age",
                  h5(strong("Age")),
                  min = 16,
                  max = 70,
                  value = 30),
      selectInput(".Analysed_Disease",
                  h5(strong("Disease")),
                  choices = list("ALL","AML","MDS","ML","MPN","ATL")),
      selectInput(".Disease.satatus.all",
                  h5(strong("Disease Status")),
                  choices = list("CR","nonCR")),
      sliderInput(".PS24",
                  h5(strong("PS")),
                  min = 0,
                  max = 4,
                  value = 0),
      selectInput(".HCT.CI",
                  h5(strong("HCT-CI")),
                  choices = list("0" = 0, "1" = 1, "2" = 2, "over 3" = 3)),
      selectInput("R_CMVAntiB",
                  h5(strong("CMV antibody")),
                  choices = list("positive" = "yes","negative" = "no"))
    ),
    mainPanel(
      fluidRow(
        column(12,
          class="chart-panel",
          div(id="distPlot")
        )
      ),
      fluidRow(
        column(6,
          class="l-box",
          div(
            class = "ranking-panel",
            div(id="chart-ranking")
          )
        ),
        column(6,
          class="r-box",
          div(
            class="legend-panel",
            div(id="chart-legend")
          ),
          div(
            class="chart-desc",
            lapply(CHART_DESC, function(e) p(e))
          )
        )
      )
    )
  )
)


# Define server logic required to draw a histogram
server = function(input, output, session) {
  patient = observeEvent(input$go,{

    patient = data.frame(.Age=input$.Age,
                       .HCT.CI=input$.HCT.CI,
                       .PS24=input$.PS24,
                       .Analysed_Disease=input$.Analysed_Disease,
                       .Disease.satatus.all=input$.Disease.satatus.all,
                       R_CMVAntiB=input$R_CMVAntiB)

    calculate_10Txpattern_prediction <- function(patient){
      RFSRC_Pred_OS <- c()
      for( j in 1:10){
        patient1 <-transform(patient,Tx_pattern=c(j))
        print(patient1)
        str(patient1)
        RFSRC_Pred_OS[[(length(RFSRC_Pred_OS) + 1)]] <- predict(modelRFSRC, patient1)
        session$sendCustomMessage(
          type="calc",
          message=list(status="processing", total=10, rest=10-j))
      }
      return(RFSRC_Pred_OS)
    }
    
    extract_data <- function (preds) {
      x <- list(c('x', preds[[1]]$time.interest))
      data <- lapply(seq_along(preds), function (i, ds) {
        return (c(Tx_pattern.names[i], ds[[i]]$survival))
      }, ds=preds)
      return (c(x, data))
    }
    
    session$sendCustomMessage(
      type='calc',
      message=list(status="start", total=10, rest=10))
    
    RFSRC_Pred_OS <- calculate_10Txpattern_prediction(patient)
    data <- extract_data(RFSRC_Pred_OS)
    
    session$sendCustomMessage(
      type='calc',
      message=list(status="end", data=data))


    # When you update data.csv (= initial data for the OS prediction curve), comment out below and push calculation.
    # updateInitData(RFSRC_Pred_OS)

  })
}  

# Run the application 
shinyApp(ui = ui, server = server)
