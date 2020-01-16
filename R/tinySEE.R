#' tinySEE
#'
#' A very simple shiny app for exploring a SummarizedExperiment, allowing the
#' browsing and plotting of genes (essentially an interface to the SEtools
#' heatmaps). For more functionality, see the `iSEE` package.
#'
#' @param se An optional object of class `SummarizedExperiment`. If omitted,
#' users will be able to upload one.
#' @param ... Any parameters passed to `shinyApp`
#'
#' @return a shiny app
#' @export
#'
#' @examples
#' data(SE)
#' app <- tinySEE(SE)
#'
#' @import shiny
tinySEE <- function(se=NULL, ...){
    shiny::shinyApp(ui=.tinySEE_ui(), server=.tinySEE_server(se), ...)
}

#' @import shiny shinydashboard shinycssloaders
#' @importFrom DT datatable renderDT DTOutput
.tinySEE_ui <- function(){
    shinyUI( dashboardPage(

      dashboardHeader(title="tinySEE"),
      dashboardSidebar(
        sidebarMenu(
          menuItem("File input", tabName="tab_fileinput"),
          menuItem("Samples", tabName="tab_samples"),
          menuItem("Features", tabName="tab_features"),
          menuItem("Plot gene", tabName="tab_gene"),
          menuItem("Heatmap", tabName="tab_heatmap")
        )
      ),
      dashboardBody(
        tabItems(
          tabItem("tab_fileinput",
            box(width=7,
                fileInput( "file", "Choose SE .rds file", multiple = FALSE,
                           accept=c(".rds")),
                withSpinner(textOutput("fileout")),
                withSpinner(verbatimTextOutput("fileout2")) )
          ),
          tabItem("tab_samples",
                  box(width=12, tags$div(style="overflow: auto;",
    	                        withSpinner(DTOutput("samples"))) )
                 ),
          tabItem("tab_features",
                  box(width=12, tags$div(style="overflow: auto;",
    	                        withSpinner(DTOutput("features")) ),
    	                        shiny::actionButton("transfer_geneSel",
    	                   "Transfer filtered genes (max 500) to the heatmap"))
                  ),
          tabItem("tab_gene", box( width=12,collapsible = TRUE,
               column(6, selectizeInput("gene_input", "Select Gene",
                                        choices=c(), multiple=FALSE),
                      selectInput("assay_input", "Assay", choices=c(),
                                  multiple=FALSE)
                      ),
               column(3, selectInput("plottype_input", "Type of Plot",
                                     choices=c("violin plot","box plot"),
                                     multiple=FALSE),
                      checkboxInput('select_plotpoints','Plot Points',
                                    value=TRUE),
                      checkboxInput('select_logaxis','Logarithmic Axis',
                                    value=FALSE),
                      numericInput("gp_height", "Plot height", 400, min=100,
                                   max=1000, step=50)
                      ),
               column(3, selectInput("select_groupvar", "Group by",
                                     choices=c(), multiple=FALSE),
                      checkboxInput('asfactor','As factor', value=TRUE),
                      selectInput("select_colorvar", "Color by", choices=c(),
                                  multiple=FALSE),
                      selectizeInput("select_gridvars", "Grid by", choices=c(),
                                     options=list(maxItems=2), multiple=TRUE),
                      checkboxInput('select_freeaxis','Free Axis', value=TRUE)
                      )
            ),
            box(width=12, withSpinner(plotOutput("gene_plot", height="auto")),
                collapsible = TRUE)
          ),
          tabItem("tab_heatmap",
               box( width=12, title="Genes", collapsible=TRUE,
                    textAreaInput('input_genes','Genes to plot', width="90%",
                                  rows=7, placeholder=
           "Enter genes symbols separated by commas, spaces, or line breaks...")
               ),
               box( width=12, title="Heatmap parameters", collapsible=TRUE,
                    column(4, selectInput("assay_input2", "Assay", choices=c(),
                                          multiple=FALSE),
                           checkboxInput('hm_scale', 'Scale rows'),
                           checkboxInput('hm_breaks', 'Symmetric scale')
                           ),
                    column(4, selectizeInput('hm_anno', "Column annotation",
                                             choices=c(), multiple=TRUE),
                           selectizeInput('hm_gaps', "Gaps at", choices=c(),
                                          multiple=T)
                    ),
                    column(4, selectizeInput('hm_order', "Column ordering",
                                             choices=c(), multiple=TRUE),
                           checkboxInput('hm_clusterCol','Cluster columns',
                                         value=F),
                           checkboxInput('hm_clusterRow','Sort rows',
                                         value=TRUE),
                           numericInput("hm_height", "Plot height", 400,
                                        min=100, max=2000, step=50)
                    )
               ),
               box(width=12, title="Heatmap", collapsible = TRUE,
                   withSpinner(plotOutput("heatmap", height="auto")) )
          )
      )
    )))
}


#' @import shiny shinydashboard shinycssloaders
#' @importFrom DT datatable renderDT DTOutput
.tinySEE_server <- function(se=NULL, maxSize=50*1024^2){
    library(shiny)
    library(DT)
    library(ggplot2)
    library(cowplot)
    theme_set(theme_cowplot())
    library(shinycssloaders)
    library(SummarizedExperiment)

    options(shiny.maxRequestSize=maxSize)

    grepGene <- function(x,g){
        if(!is.character(x)){
            g <- grepGene(row.names(x), g)
            return(x[g,drop=FALSE])
        }
        if(all(g %in% x)) return(g)
        g <- paste0("^",g,"\\.|^",g,"$|\\.",g,"$")
        g <- lapply(g,FUN=function(i) grep(i, x, value=TRUE, ignore.case=TRUE))
        return(unique(unlist(g)))
    }

    getDef <- function(se,var){
        if(length(var)>1){
            y <- unlist(lapply(var, FUN=function(x) getDef(se,x)))
            y <- y[!sapply(y,is.null)]
            if(length(y)==0) return(NULL)
            return(y)
        }
        if(is.null(se@metadata$default_view[[var]])) return(NULL)
        se@metadata$default_view[[var]]
    }

    shinyServer(function(input, output, session) {

        SEinit <- function(x){
            if(is.null(assayNames(x)))
                assayNames(x) <- paste0("assay",1:length(assays(x)))
            if(ncol(rowData(x))==0) rowData(x)$name <- row.names(x)
            updateSelectizeInput(session, "gene_input",
                                 choices=sort(unique(row.names(x))),
                                 server=TRUE)
            colvars <- colnames(colData(x))
            updateSelectInput(session, "assay_input",
                             choices=assayNames(x), selected=getDef(x, "assay"))
            updateSelectInput(session, "assay_input2",
                             choices=assayNames(x), selected=getDef(x, "assay"))
            updateSelectizeInput(session, "hm_order", choices=colvars)
            updateSelectizeInput(session, "hm_anno", choices=colvars,
                        selected=getDef(x, c("groupvar","colvar","gridvar")))
            updateSelectizeInput(session, "hm_gaps", choices=colvars,
                                 selected=getDef(x, "gridvar"))
            updateSelectInput(session, "select_groupvar", choices=colvars,
                              selected=getDef(x, "groupvar"))
            updateSelectInput(session, "select_colorvar", choices=colvars,
                              selected=getDef(x, "colvar"))
            updateSelectInput(session, "select_gridvars", choices=colvars,
                              selected=getDef(x, "gridvar"))
            if(!is.null(getDef(x, "assay")))
                updateCheckboxInput(session, "hm_scale",
                                    value=!grepl("FC$",getDef(x, "assay")))
            return(x)
        }

        observeEvent(input$transfer_geneSel, {
            g <- input$features_rows_all
            if(length(g)>0){
                g <- g[seq_len(min(length(g),500))]
                g <- row.names(SE())[g]
                g <- paste(g, collapse=", ")
                updateTextAreaInput(session, "input_genes", value=g)
            }
        })

        SE <- reactive({
            tryCatch({
                if(!is.null(input$file)){
                    x <- readRDS(input$file$datapath)
                    if(is(x,"SummarizedExperiment")){
                        return(SEinit(x))
                    }
                    stop("Object is not a SummarizedExperiment!")
                }else{
                    if(!is.null(se)){
                        return(SEinit(se))
                    }else if(file.exists("default.SE.rds")){
                        x <- readRDS("default.SE.rds")
                        return(SEinit(x))
                    }
                    return(NULL)
                }
            }, error=function(e){ stop(safeError(e)) }
            )
        })

        output$fileout <- renderText({
            if(is.null(SE())) return(NULL)
            tout <- "Loaded: "
            if(!is.null(SE()@metadata$title))
                tout <- paste(tout, SE()@metadata$title)
            tout
        })

        output$fileout2 <- renderPrint({
            if(is.null(SE())) return(NULL)
            print(SE())
        })

        output$features <- renderDT({
            if(is.null(SE())) return(NULL)
            datatable( as.data.frame(rowData(SE())), filter="top",
                       options=list( pageLength=30, dom = "fltBip" ),
                       extensions=c("ColReorder") )
        }, server = TRUE)

        output$samples <- renderDT({
            if(is.null(SE())) return(NULL)
            datatable( as.data.frame(colData(SE())), filter="top",
                       options=list( pageLength=30, dom = "fltBip" ),
                       extensions=c("ColReorder") )
        })

        ############
        ### BEGIN HEATMAP

        selGenes <- reactive({
            g <- gsub(",|\n|\r|;,"," ",input$input_genes)
            g <- strsplit(g," ",fixed=TRUE)[[1]]
            unique(g[which(g!="")])
        })

        output$heatmap <- renderPlot({
            if(is.null(SE())) return(NULL)
            g <- grepGene(row.names(SE()), selGenes())
            if(length(g)==0) return(NULL)
            if(length(g)>2 && input$hm_clusterRow){
                srow <- seq_length(ncol(SE()))
            }else{
                srow <- NULL
            }
            se <- SE()
            o <- input$hm_order
            if(is.null(o)) o <- c()
            for(f in rev(o)) se <- se[,order(colData(se)[[f]])]

            sehm(se, g, input$hm_scale, assayName=input$assay_input2,
                 sortRowsOn=srow, anno_columns=input$hm_anno,
                 gaps_at=input$hm_gaps, cluster_cols=input$hm_clusterCol,
                 cluster_rows=FALSE, breaks=input$hm_breaks)

        }, height=reactive(input$hm_height))

        ### END HEATMAP
        ############




        ############
        ### START GENE TAB
        output$gene_plot <- renderPlot({
            d <- tryCatch(meltSE(SE(), input$gene_input),
                          error=function(x) NULL)
            if(is.null(d)) return(NULL)
            gr <- input$select_groupvar
            if(input$asfactor) d[[gr]] <- factor(d[[gr]])
            if(input$select_plotpoints){
                p <- ggplot(d, aes_string(input$select_groupvar,
                                          input$assay_input,
                                          colour=input$select_colorvar))
            }else{
                p <- ggplot(d, aes_string(input$select_groupvar,
                                          input$assay_input,
                                          fill=input$select_colorvar))
            }
            if(input$plottype_input=="violin plot"){
                p <- p + geom_violin()
            }else{
                p <- p + geom_boxplot(outlier.shape = NA)
            }
            if(input$select_plotpoints)
                p <- p + geom_point(position = position_jitterdodge())
            p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle(input$gene_input)

            if(!is.null(input$select_gridvars)){
                form <- as.formula(paste0("~", paste(input$select_gridvars,
                                                     collapse="+")))
                if(input$select_freeaxis){
                    p <- p + facet_wrap(form, scales="free_y")
                }else{
                    p <- p + facet_wrap(form)
                }
            }
            p
        }, height=reactive(input$gp_height))

        ### END GENE TAB
        ############

    })


}
