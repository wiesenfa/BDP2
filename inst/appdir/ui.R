library(shinyBS)
library(BDP2)

fluidPage(
    # oben  
      titlePanel(paste0("BDP2 workflow (version ",packageVersion('BDP2'),")")),
      h3("Workflow to determine design parameters for a multi-stage single-arm phase II trial with binary endpoint. Declaration of efficacy and futility is based on Bayesian posterior distribution."),
      h4("Annette Kopp-Schneider, Manuel Wiesenfarth, Division of Biostatistics, German Cancer Research Center (DKFZ), Heidelberg, Germany"),
      br(),
      h4("1. Common settings: Select pF, pE and the prior distribution for the response rate"),#br(),
        fluidRow(
          column(3,
            h5("Uninteresting and target response rate"),
            numericInput("p0", "p0:", value=.12, min = 0, max = 1, step=.02),bsTooltip("p0", "Uninteresting (control) response rate", "right", options = list(container = "body")),
            numericInput("p1", "p1:", value=.3, min = 0, max = 1, step=.02),bsTooltip("p1", "Target response rate", "right", options = list(container = "body"))
          ),
          column(2,
            h5("pF and pE"),
            numericInput("pF", "pF:", value=.3, min = 0, max = 1, step=.02), bsTooltip("pF", "Response rate used for futility criterion P(p>pF|Data)", "right", options = list(container = "body")),
            numericInput("pE", "pE:", value=.12, min = 0, max = 1, step=.02),bsTooltip("pE", "Response rate used for efficacy criterion P(p>pE|Data)", "right", options = list(container = "body"))
          ),
          column(2,
            h5("Beta prior for futility criterion"),
            numericInput("shape1F", "shape1F:", value=.3, step=.02),bsTooltip("shape1F", "shape1 parameter for Beta prior for futility criterion", "right", options = list(container = "body")),
            numericInput("shape2F", "shape2F:", value=.7,  step=.02),bsTooltip("shape2F", "shape2 parameter for Beta prior for futility criterion", "right", options = list(container = "body"))
          ),
          column(2,
            h5("Beta prior for efficacy criterion"),
            numericInput("shape1E", "shape1E:", value=.12, step=.02),bsTooltip("shape1E", "shape1 parameter for Beta prior for efficacy criterion", "right", options = list(container = "body")),
            numericInput("shape2E", "shape2E:", value=.88,  step=.02),bsTooltip("shape2E", "shape2 parameter for Beta prior for efficacy criterion", "right", options = list(container = "body"))
          )
        ),
      br(),hr(),br(),#br(),

    sidebarLayout(
    # unten links
      sidebarPanel(
        h5("Note: For support for settings and for effect of settings on trial characteristics click on tabs to the right!"),hr(),
        h4( "2. Select timing of interims and cF"),
          numericInput("firstInterim.at", "First interim at:", value=10, min = 0, step=1),
          numericInput("cF", "cF (for futility criterion P(p>pF|Data) < cF):", value=.01, min = 0, max = 1, step=.01),
          br(),
          textInput("furtherInterims.at", "Second and potential further interims at (separate by blank):", value = "20"),
          hr(),
        h4("3. Select cE"),
          numericInput("cE", "cE (for efficacy criterion P(p>pE|Data) >= cE):", value=.9, min = 0, max = 1, step=.01),
          hr(),
        downloadButton("report", "Generate report"),
        radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
inline = TRUE)
      ),
    # unten rechts
      mainPanel(
        tabsetPanel(
          tabPanel(title="2. Select timing of interims and cF",
            h4("2a. Select timing of first interim analysis and futility stopping criterion, cF, such that the trial stops after a run of treatment failures:"),
              sliderInput("n.range_1", "Range of sample sizes at first interim:", min = 0, max = 200, value = c(4,15),width = 400),
              plotOutput("ProbSuccesses"),
              h5("The left figure shows true (for ptrue=p0, in green) and false (for ptrue=p1, in red) stopping probability. 
                 False stopping probability should be bounded to an acceptable limit, e.g., 5%. 
                 Any n with false stopping probability below this limit can be used as timing of first interim. 
                 The actual choice may be guided by logistical considerations."),
              h5("From the right figure read P(p>pF| 0 of n), in red, and P(p>pF| 1 of n), in green, for your selected timing of the first interim. 
                 The value of cF can be chosen anywhere between these two numbers and 
                  will ensure that the trial is stopped at first interim if no success has been observed so far.
                 Consider the operating characteristics (tab 4), power and expected patient number (tab 5) for final decision on cF."),
              h5("Note: If it is intended to stop the trial not only after a run of treatment failures but also in case 
                 a low number of successes is observed, the timing of first interim and cF can be chosen on the basis of the 
                 operating characteristics (tab 4), power and expected patient number (tab 5) without evaluating these plots.",style = "color:blue"),
              h4("2b. Timings of second and potential further interims are selected based on logistical considerations."),
              br()
             ),


          tabPanel(title="3. Select cE",
            h3("3. Selecting cE"),br(),
              fluidRow(
                column(2,numericInput("nfinal", "final analysis at :", value=20, min = 0, step=1)),
                column(2,sliderInput("cE.range_1", "Range of cE values:", min = 0, max = 1, value = c(0.75,0.99),width = 400))
              ),
              plotOutput("cE.vs.PEcall"),
              h5("The plot shows power (in red) and type I error (in green) as function of cE."),
              h5("Select threshold cE such that type I error probability can be tolerated and power is maximized.") ,
              textOutput("OCs.selected.cE"),
            br()
          ),
          
          tabPanel(title="4. Operating characteristics",
            h3("4. Investigating operating characteristics"),br(),
              sliderInput("n.range_2", "Range of final sample sizes n:", min = 0, max = 200, value = c(15,30),width = 400),
              plotOutput("n.vs.PFstopEcall"),
            h5("The left plot shows the probabilities for declaration of efficacy and cumulative stopping for futility at final for ptrue=p0."), 
            h5("The right plot shows the same quantities for ptrue=p1. "),
            br()
          ),
          
          tabPanel(title="5. Power and expected number of patients for multiple sample sizes",
            h3("5. Power function for multiple sample sizes at final analysis and corresponding expected number of patients in the trial as a function of ptrue"),br(),
              fluidRow(
                column(2,textInput("nfinal.vec", "final sample sizes to compare", value = "20 30")),
                #    column(2,numericInput("interim.atEvery", "interim at every ... patient :", value=10, min = 0, step=1)),
                column(2,sliderInput("ptrue.range_1", "Range of ptrue values:", min = 0, max = 1, value = c(0,0.5),width = 400))
              ),
              plotOutput("ptrue.vs.PEcall"),
              #       plotOutput("ptrue.vs.ExpectedNumber")
              br()
          ),
          
          tabPanel(title="6. Decision boundaries",
            h3("6. Decision boundaries for futility and efficacy (Design specification for clinician) and check for contradictory results"),
              sliderInput("n.range_4", "Range of final sample sizes n:", min = 0, max = 200, value = c(0,40),width = 400),
              plotOutput("n.vs.bFbE"),
            h5("The plot shows the decision boundaries for the selected trial design in terms of number of successes per number of enrolled patients."),
            h5("The trial is stopped for futility if the observed number of successes among enrolled patients is in the red part of the plot, including the boundary."), 
            h5("Efficacy can be called if the observed number of successes among patients at final analysis is in the green part of the plot, including the boundary."),
            h5("The range of final sample size should be increased to the maximal feasible number to check for contradictory results, i.e., for overlap between red and green areas."),
            br()
          )
        )

      )
    ),br(),br(),
h4("Reference"),
  h5("Kopp-Schneider, A., Wiesenfarth, M., Witt, R., Edelmann, D., Witt, O. and Abel, U. (2018). Monitoring futility and efficacy in phase II trials with Bayesian
posterior distributions - a calibration approach. 
Biometrical Journal, to appear.")
  )
