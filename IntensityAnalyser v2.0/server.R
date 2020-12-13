library('shiny')
library('shinyjs')
library('readxl')
library('ggplot2')
library('rhandsontable')
library('rmarkdown')
library('plyr')
library('dplyr')
library('forcats')
library('splines')
library('reshape2')
library('emmeans')
library('multcompView')
library('multcomp')
library('PMCMRplus')
library('showtext')
library('showtextdb')
library('sysfonts')
library('curl')

#Function for finding peak maximum
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

shinyServer(function(input, output, session){
  shinyjs::disable(selector = '.navbar-nav a')
  #Data input
  filedata<- reactive({
    infile <- input$uploadfile
    if (is.null(infile)){
      return(NULL)      
    }
    data<- read_xlsx(infile$datapath, sheet =  1)
    data$pixel<- as.numeric(data$pixel)
    return(data)
  })
  
  #Upload plot font
  sysfonts::font_add_google('Roboto', 'Roboto')
  showtext_auto()

  #Run program button for moving to marker and enabling other tabs
  output$ui.run_program <- renderUI({
    if (is.null(filedata())) return()
    actionButton("run_program", "Run Program", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  })
  
  observeEvent(input$run_program, {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "marker")
    shinyjs::enable(selector = '.navbar-nav a')
    shinyjs::disable(selector = ".navbar-nav a[data-value=upload_file]")
    shinyjs::disable(selector = ".navbar-nav a[data-value=graphic_result]")
    shinyjs::disable(selector = ".navbar-nav a[data-value=pxmw_ratio]")
    shinyjs::disable(selector = ".navbar-nav a[data-value=length_calc]")
    shinyjs::disable(selector = ".navbar-nav a[data-value=marker_correction]")
  }) 
  
  #Dataframe to fill with bp and peakno in markertable (made as 1 kb ladder as default)
  marker<- reactiveValues()
  kbladder<- c(10000,8000,6000,5000,4000,3500,3000,2500,2000,1500,1000,750,500,250)
  peak0<- as.numeric(rep(NA,14))
  dftofill<- setNames(data.frame(kbladder, peak0), c("band size [bp]", "peak no."))
  marker$dataframe_marker<- reactiveVal({dftofill})
  
  #Observe event for the change of marker$dataframe if changed to different marker
  observeEvent(input$select_marker, {
    if (input$select_marker == 1) {
      kbladder<- c(10000,8000,6000,5000,4000,3500,3000,2500,2000,1500,1000,750,500,250)
      peak0<- as.numeric(rep(NA,14))
      dftofill<- setNames(data.frame(kbladder, peak0), c("band size [bp]", "peak no."))
      marker$dataframe_marker(dftofill)
    }
    if (input$select_marker == 2) {
      kbladder<- c(20000,10000,7000,5000,4000,3000,2000,1500,1000,700,500,400,300,200,75)
      peak0<- as.numeric(rep(NA,15))
      dftofill<- setNames(data.frame(kbladder, peak0), c("band size [bp]", "peak no."))
      marker$dataframe_marker(dftofill)
    }
    if (input$select_marker == 3) {
      kbladder<- c(20000,10000,8000,7000,6000,5000,4000,3500,3000,2500,2000,1500,1000,750,700,500,400,300,250,200,75)
      peak0<- as.numeric(rep(NA,21))
      dftofill<- setNames(data.frame(kbladder, peak0), c("band size [bp]", "peak no."))
      marker$dataframe_marker(dftofill)
    }
    if (input$select_marker == 4) {
      kbladder<- c(48502,24508,20555,17000,15258,13825,12119,10171)
      peak0<- as.numeric(rep(NA,8))
      dftofill<- setNames(data.frame(kbladder, peak0), c("band size [bp]", "peak no."))
      marker$dataframe_marker(dftofill)
    }
    if (input$select_marker == 5) {
      kbladder<- c(48502,24508,20555,17000,15258,13825,12119,10171,10000,8000,6000,5000,4000,3500,3000,2500,2000,1500,1000,750,500,250)
      peak0<- as.numeric(rep(NA,22))
      dftofill<- setNames(data.frame(kbladder, peak0), c("band size [bp]", "peak no."))
      marker$dataframe_marker(dftofill)
    }
    if (input$select_marker == 6) {
      kbladder<- c(48502,24508,20555,20000,17000,15258,13825,12119,10171,10000,7000,5000,4000,3000,2000,1500,1000,700,500,400,300,200,75)
      peak0<- as.numeric(rep(NA,23))
      dftofill<- setNames(data.frame(kbladder, peak0), c("band size [bp]", "peak no."))
      marker$dataframe_marker(dftofill)
    }
    if (input$select_marker == 7) {
      kbladder<- c(727500,679000,630500,582000,533500,485000,436500,388000,339500,291000,242500,194000,145500,97000,48500)
      peak0<- as.numeric(rep(NA,15))
      dftofill<- setNames(data.frame(kbladder, peak0), c("band size [bp]", "peak no."))
      marker$dataframe_marker(dftofill)
    }
    if (input$select_marker == 8) {
      kbladder<- c(23130,9416,6557,4361,2322,2027,564,125)
      peak0<- as.numeric(rep(NA,8))
      dftofill<- setNames(data.frame(kbladder, peak0), c("band size [bp]", "peak no."))
      marker$dataframe_marker(dftofill)
    }
  })
  
  #Editable dataframe by rhandsontable for marker bp and peak no.
  output$markertable<- renderRHandsontable({
    rhandsontable(marker$dataframe_marker(), rowHeaderWidth = 0)%>%
      hot_cols(colWidths = 100)%>%
      hot_cols(format = "0")
  })
  
  #Function that has saved dataframe later used for creating graph
  maxima_dataframe<- reactive({
    #Taking pixel and marker data out of dataframe
    pixel<- filedata()$pixel[!is.na(filedata()$marker) & !is.na(filedata()$pixel)]
    pixel<- pixel[-1]
    marker<- filedata()$marker[!is.na(filedata()$marker) & !is.na(filedata()$pixel)]
    marker<- marker[-1]
    
    #Smoothing
    marker_smooth <- loess(marker ~ pixel, span=0.02)$fitted
    
    #Creating false values for highlighting maxima in graph
    marker_TF <- marker_smooth
    marker_TF[find_peaks(marker_smooth, m = input$sens)] <- 0
    marker_TF[marker_TF < 0] <- 1
    marker_TF<- as.logical(marker_TF)
    
    #Creating dataframe
    marker_dataframe<- data.frame(pixel[1:length(marker_smooth)], marker_smooth, marker_TF)
    
    #Naming maxima in graph by numbers
    marker_TF_length<- (1:length(marker_dataframe$marker_TF[marker_dataframe$marker_TF == FALSE]))
    marker_dataframe$name <- 0
    marker_dataframe$name[marker_dataframe$marker_TF == FALSE] <- marker_TF_length
    
    #Renaming dataframe
    colnames(marker_dataframe)<- c("pixel", "marker_smooth", "marker_TF", "name")
    return(marker_dataframe)
  })
  
  #ggplot2 function for showing maxima of graph (used later as output)
  maxima_graph<- reactive({
    ggplot2::ggplot(maxima_dataframe(), aes(x = maxima_dataframe()$pixel ,y= maxima_dataframe()$marker_smooth)) +
      geom_line(colour= "grey") +
      geom_point(alpha= 0.6, colour= "grey", size= 2) +
      geom_point(alpha= 1, aes(colour= maxima_dataframe()$marker_TF == FALSE, size = maxima_dataframe()$marker_TF == FALSE)) +
      geom_text(label=ifelse(maxima_dataframe()$marker_TF == FALSE,as.character(maxima_dataframe()$name),''), vjust= -0.6, hjust= -0.15, family="Roboto") +
      scale_size_manual(values =c(-10, 3)) +
      scale_colour_manual(values = c("grey", "#4493d3")) +
      scale_y_continuous(expand = c(0,0),
                         limits = c(0,(max(maxima_dataframe()$marker_smooth) + max(maxima_dataframe()$marker_smooth)*0.07))) +
      ylab("marker intensity") +
      xlab("pixel") +
      theme_light() +
      theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), axis.line = element_line(colour = "grey"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major=element_blank(), panel.border = element_blank())
  })
  
  output$maxima_plot_for_marker_tab <- output$maxima_plot <- renderPlot(maxima_graph())
  
  #Saving markertable for the report
  dataframe_marker<- reactive({
    input_markertable<- as.data.frame(hot_to_r(input$markertable))
    input_markertable<- input_markertable[!is.na(input_markertable[[2]]) & !is.na(input_markertable[[1]]),]
    return(input_markertable)
  })
  
  #Proceed without marker correction -> proceed button with first marker to the pixel/MW ratio
  output$ui.proceed_without_correction <- renderUI({
    actionButton("proceed_not_correct", "Proceed without marker correction", style="color: #337ab7")
  })
  
  observeEvent(input$proceed_not_correct, {
    if (anyNA(dataframe_polynom()) == TRUE) {
      shinyalert("Non-existant peak maxima chosen", "Check peak no. column. It includes value of peak maxima not existing in the graph.", type = "error", confirmButtonCol =  "#337ab7")
    } else if (dim(dataframe_polynom())[1] == 0) {
      shinyalert("Table is not filled", "Fill out the table before continuing.", type = "error", confirmButtonCol =  "#337ab7")
    } else {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "pxmw_ratio")
    shinyjs::disable(selector = ".navbar-nav a[data-value=pxmw_ratio]")
    shinyjs::disable(selector = ".navbar-nav a[data-value=marker_correction]")
    values$marker_correction<- FALSE
    }
  })
  
  #Choose polynomial button to take user to the Marker correction tab
  output$ui.proceed_with_correction <- renderUI({
    actionButton("proceed_to_correct", "Proceed with marker correction", style="color: #ff4d4d")
  })
  
  observeEvent(input$proceed_to_correct, {
    if (anyNA(dataframe_polynom()) == TRUE) {
      shinyalert("Non-existant peak maxima chosen", "Check peak no. column. It includes value of peak maxima not existing in the graph.", type = "error", confirmButtonCol =  "#337ab7")
    } else if (dim(dataframe_polynom())[1] == 0) {
      shinyalert("Table is not filled", "Fill out the table before continuing.", type = "error", confirmButtonCol =  "#337ab7")
    } else {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "marker_correction")
    shinyjs::disable(selector = ".navbar-nav a[data-value=marker_correction]")
    shinyjs::disable(selector = ".navbar-nav a[data-value=pxmw_ratio]")
    }
  })
  
  
  
  ##
  ## START MARKER CORRECTION PART - parts related to the Marker correction tab (selecting second marker from samples, plot of the sample, sensitivity correction##
  ##
  
  observeEvent(input$info_marker_correction_tab, {
    wan_et_al <- 
      shinyalert(html = TRUE,paste(
        h6(tags$div("Here you can correct for the marker shift if you applied two markers on the gel:", tags$br(),
                    "1) Select the sample corresponding to the second marker.", tags$br(),
                    "2) Check the graph. Displayed points on the highest peaks should correspond to the bands of the ladder.", tags$br(),
                    "3) Write the peak number corresponding to the respective band in the table. In the row, peaks from both markers needs to correspond to the same band."), align = "justify"),
        h6("Note: You need to fill the whole table unlike in the first marker table on previous the Marker tab. In case that you do not have the same number of bands visible in the second marker,
           remove those bands that are not present in the second marker from the table in the Marker tab.", align = "justify")
      ), type = "info", confirmButtonCol =  "#337ab7")
  })
  
  #slider for selecting sample as second marker
  output$slider_second_marker_option <- renderUI({
    sliderInput("sample_slider_for_second_marker", "Which sample is the second marker?", min=1, max=length(colnames(filedata()))-2, value=length(colnames(filedata()))-2, step = 1, ticks = TRUE, width = 296)
  })
  #> button for sliding on slider for number of samples forward
  observeEvent(input$forward_button_slider_for_second_marker, {
    if(input$sample_slider_for_second_marker < length(colnames(filedata()))-2) {
      updateSliderInput(session, "sample_slider_for_second_marker", value=input$sample_slider_for_second_marker+1)
    }
  })
  
  #< button for sliding on slider for number of samples backward
  observeEvent(input$backward_button_slider_for_second_marker, {
    if(input$sample_slider_for_second_marker > 1) {
      updateSliderInput(session, "sample_slider_for_second_marker", value=input$sample_slider_for_second_marker-1)
    }
  })
  
  #Function that has saved dataframe later used for creating graph
  maxima_dataframe_marker_correction<- reactive({
    #Taking pixel and marker data out of dataframe
    pixel<- filedata()$pixel
    pixel<- pixel[-1]
    marker<- filedata()[[input$sample_slider_for_second_marker+2]]
    marker<- marker[-1]
    #Smoothing
    marker_smooth <- loess(marker ~ pixel, span=0.02)$fitted
    
    #Creating false values for highlighting maxima in graph
    marker_TF <- marker_smooth
    marker_TF[find_peaks(marker_smooth, m = input$sens_second_marker)] <- 0
    marker_TF[marker_TF < 0] <- 1
    marker_TF<- as.logical(marker_TF)
    
    #Creating dataframe
    marker_dataframe<- data.frame(pixel[1:length(marker_smooth)], marker_smooth, marker_TF)
    
    #Naming maxima in graph by numbers
    marker_TF_length<- (1:length(marker_dataframe$marker_TF[marker_dataframe$marker_TF == FALSE]))
    marker_dataframe$name <- 0
    marker_dataframe$name[marker_dataframe$marker_TF == FALSE] <- marker_TF_length
    
    #Renaming dataframe
    colnames(marker_dataframe)<- c("pixel", "marker_smooth", "marker_TF", "name")
    return(marker_dataframe)
  })
  
  #ggplot2 function for showing maxima of sample graph for selecting points of the marker correction (used later as output)
  maxima_graph_second_marker<- reactive({
    ggplot2::ggplot(maxima_dataframe_marker_correction(), aes(x = maxima_dataframe_marker_correction()$pixel ,y= maxima_dataframe_marker_correction()$marker_smooth)) +
      geom_line(colour= "grey") +
      geom_point(alpha= 0.6, colour= "grey", size= 2) +
      geom_point(alpha= 1, aes(colour= maxima_dataframe_marker_correction()$marker_TF == FALSE, size = maxima_dataframe_marker_correction()$marker_TF == FALSE)) +
      geom_text(label=ifelse(maxima_dataframe_marker_correction()$marker_TF == FALSE,as.character(maxima_dataframe_marker_correction()$name),''), vjust= -0.6, hjust= -0.15, family="Roboto") +
      scale_size_manual(values =c(-10, 3)) +
      scale_colour_manual(values = c("grey", "#4493d3")) +
      scale_y_continuous(expand = c(0,0),
                         limits = c(0,(max(maxima_dataframe_marker_correction()$marker_smooth) + max(maxima_dataframe_marker_correction()$marker_smooth)*0.07))) +
      ylab("marker intensity") +
      xlab("pixel") +
      theme_light() +
      theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), axis.line = element_line(colour = "grey"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major=element_blank(), panel.border = element_blank())
  })
  
  output$maxima_plot_second_marker <- renderPlot(maxima_graph_second_marker())
  
  #Dataframe to fill with second and group name in markertable
  both_markers_dataframe<- reactive({
    input_markertable<- as.data.frame(hot_to_r(input$markertable))
    input_markertable<- input_markertable[!is.na(input_markertable[[2]]) & !is.na(input_markertable[[1]]),]
    bm_dataframe<- data.frame(input_markertable[[1]], input_markertable[[2]], as.numeric(rep(NA, length(input_markertable[[1]]))), stringsAsFactors = FALSE)
    bm_dataframe<- setNames(bm_dataframe, c("band size [bp]", "1st marker peak no.", "peak no."))
    return(bm_dataframe)
  })
  
  #rhandsontable to fill the 2nd marker maxima corresponding to the bands
  output$second_marker_table<- renderRHandsontable({
    rhandsontable(both_markers_dataframe(), rowHeaderWidth = 0, contextMenu = FALSE)%>%
      hot_cols(colWidths = 89)%>%
      hot_cols(format = "0")%>%
      hot_col("band size [bp]", readOnly = TRUE)%>%
      hot_col("1st marker peak no.", readOnly = TRUE)
  })
  
  #Taking rhandsontable data and creating dataframe for later calculation of the x,y for each sample
  dataframe_pixel_pos_both_markers<- reactive({
    input_markertable<- as.data.frame(hot_to_r(input$markertable))
    input_markertable<- input_markertable[!is.na(input_markertable[[2]]) & !is.na(input_markertable[[1]]),]
    input_second_marker_table<- as.data.frame(hot_to_r(input$second_marker_table))
    input_second_marker_table<- input_second_marker_table[!is.na(input_second_marker_table[[2]]) & !is.na(input_second_marker_table[[1]]),]
    df_pixel_pos_both_markers<- data.frame(maxima_dataframe()$pixel[maxima_dataframe()$marker_TF == FALSE][input_markertable[[2]]],
                                           maxima_dataframe_marker_correction()$pixel[maxima_dataframe_marker_correction()$marker_TF == FALSE][input_second_marker_table[[3]]],
                                           input_markertable[[1]])
    colnames(df_pixel_pos_both_markers)<- c("pixel_pos_marker_one", "pixel_pos_marker_two", "band_size")
    return(df_pixel_pos_both_markers)
  })
  
  #Taking the pixel positions of every sample and marker to vector
  row_with_pixel_gel_positions<- reactive({filedata()[1,-1]}) 
  
  #Dataframe with values calculated using linear model to predict the y-pixel position based on the x-position on the TRF scan
  predicted_markers_datatable<- reactive({
    predicted_markers_dataframe<- dataframe_pixel_pos_both_markers()
    for (i in 1:(ncol(row_with_pixel_gel_positions())-1)) {
      pixel_band_value<- NA
      for (j in 1:nrow(predicted_markers_dataframe)) {
        x<- as.numeric(c(row_with_pixel_gel_positions()[[1]], row_with_pixel_gel_positions()[colnames(row_with_pixel_gel_positions()) == as.character(input$sample_slider_for_second_marker)][1,1]))
        y<- dataframe_pixel_pos_both_markers()[j,]
        y<- as.vector(y)
        y<- as.numeric(y[-length(y)])
        linear_model<- lm(y ~ x)
        pixel_band_value[j]<- round(linear_model$coefficients[[1]]+linear_model$coefficients[[2]]*row_with_pixel_gel_positions()[colnames(row_with_pixel_gel_positions()) == as.character(i)][1,1], digits = 0)
      }
      pixel_band_value<- as.data.frame(as.numeric(pixel_band_value))
      colnames(pixel_band_value)<- as.character(i)
      predicted_markers_dataframe<- cbind(predicted_markers_dataframe, pixel_band_value)
    }
    predicted_markers_dataframe<- predicted_markers_dataframe[,-(input$sample_slider_for_second_marker+3)]
    colnames(predicted_markers_dataframe)<- c("pixel_pos_marker_one", "pixel_pos_marker_two", "band_size", 1:(ncol(predicted_markers_dataframe)-3))
    return(predicted_markers_dataframe)
  })
  
  observeEvent(input$proceed_after_second_marker, {
    if (anyNA(dataframe_pixel_pos_both_markers()[[2]]) == TRUE) {
      shinyalert("Table not filled or non-existant peak maxima chosen", "The table either not fully filled or check the peak no. column. It might include value(s) of peak maxima not existing in the graph.", type = "error", confirmButtonCol =  "#337ab7")
    } else {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "pxmw_ratio")
    values$marker_correction<- TRUE
    }
  })
  
  ##
  ## END MARKER CORRECTION PART ##
  ##
  
  ##Info for report.rmd about senstivity value(s) and about marker table (either with or without marker correction)
  #sensitivity values
  sensitivity_vals<- reactive({
    if (values$marker_correction == FALSE) {
      input$sens
    } else {
      c(input$sens, input$sens_second_marker)
    }
  })
  #information about filled points in the rhandsontable
  markers_filled_table<- reactive({
    if (values$marker_correction == FALSE) {
      dataframe_marker()
    } else {
      input_second_marker_table<- as.data.frame(hot_to_r(input$second_marker_table))
      input_second_marker_table<- input_second_marker_table[!is.na(input_second_marker_table[[2]]) & !is.na(input_second_marker_table[[1]]),]
      colnames(input_second_marker_table)<- c("band size [bp]", "1st marker peak no.", "2nd marker peak no.")
      return(input_second_marker_table)
    }
  })
  
  ##
  ## Sample slider in pixel/MW ratio tab
  ##
  
  #slider to show calculated polynomial regression and ANOVA statistics for each table
  output$slider_pxmw_tab <- renderUI({
    if (values$marker_correction == TRUE) {
      sliderInput("slider_sample_pxmw_tab", "Which sample you wish to see?", min=1, max=ncol(predicted_markers_datatable())-3, value=1, step = 1, ticks = TRUE, width = 296)
    }
  })
  
  #show > button
  output$forward_button_slider_pxmw_tab <- renderUI({
    if (values$marker_correction == TRUE) {
      actionButton("forward_button_slider_sample_pxmw_tab", ">", style='padding:3px')
    }
  })
  
  #show < button
  output$backward_button_slider_pxmw_tab <- renderUI({
    if (values$marker_correction == TRUE) {
      actionButton("backward_button_slider_sample_pxmw_tab", "<", style='padding:3px')
    }
  })
  
  #> button for sliding on slider for number of samples forward
  observeEvent(input$forward_button_slider_sample_pxmw_tab, {
    if(input$slider_sample_pxmw_tab < ncol(predicted_markers_datatable())-3) {
      updateSliderInput(session, "slider_sample_pxmw_tab", value=input$slider_sample_pxmw_tab+1)
    }
  })
  
  #< button for sliding on slider for number of samples backward
  observeEvent(input$backward_button_slider_sample_pxmw_tab, {
    if(input$slider_sample_pxmw_tab > 1) {
      updateSliderInput(session, "slider_sample_pxmw_tab", value=input$slider_sample_pxmw_tab-1)
    }
  })
  
  ##
  ##
  ##
  
  
  #Taking rhandsontable data and creating dataframe for polynomial regression (but first removing possibility of error message due to different length of arguments)
  dataframe_polynom<- reactive({
    input_markertable<- as.data.frame(hot_to_r(input$markertable))
    input_markertable<- input_markertable[!is.na(input_markertable[[2]]) & !is.na(input_markertable[[1]]),]
    df_polynom<- data.frame(maxima_dataframe()$pixel[maxima_dataframe()$marker_TF == FALSE][input_markertable[[2]]],input_markertable[[1]])
    colnames(df_polynom)<- c("x", "y")
    return(df_polynom)
  })
  

  #Function for polynomial regression with correction
  polynom_corrected<- function(data, sample_number) {
    lm(as.numeric(predicted_markers_datatable()[[3]]) ~ poly(as.numeric(predicted_markers_datatable()[[sample_number]]), data, raw=TRUE))
  }
  
  #Function for polynomial regression
  polynom<- function(data) {
    lm(dataframe_polynom()$y ~ poly(dataframe_polynom()$x, data, raw=TRUE))
  }
  
  #Polynom functions
  polynom1<- reactive({polynom(1)})
  polynom2<- reactive({polynom(2)})
  polynom3<- reactive({polynom(3)})
  polynom4<- reactive({polynom(4)})  
  polynom5<- reactive({polynom(5)})
  polynom6<- reactive({polynom(6)})
  polynom7<- reactive({polynom(7)})  
  polynom8<- reactive({polynom(8)})
  polynom9<- reactive({polynom(9)}) 
  
  #mode function for selecting the most frequent polynomial number
  mode_calc <- function(data) {
    unique_data <- unique(data)
    unique_data[which.max(tabulate(match(data, unique_data)))]
  }
  
  #Finding the best polynom by anova method (polynom with lower number is compared to polynom with higher until the best polynom is found)
  best_polynom_number<- reactive({
    if (values$marker_correction == FALSE) {
      if(min(which(c(anova(polynom2(),polynom1())$`Pr(>F)`[2], anova(polynom3(),polynom2())$`Pr(>F)`[2], anova(polynom4(),polynom3())$`Pr(>F)`[2], anova(polynom5(),polynom4())$`Pr(>F)`[2], anova(polynom6(),polynom5())$`Pr(>F)`[2], anova(polynom7(),polynom6())$`Pr(>F)`[2], anova(polynom8(),polynom7())$`Pr(>F)`[2], anova(polynom9(),polynom8())$`Pr(>F)`[2]) > 0.05)) == Inf) {
        max(which(c(anova(polynom2(),polynom1())$`Pr(>F)`[2], anova(polynom3(),polynom2())$`Pr(>F)`[2], anova(polynom4(),polynom3())$`Pr(>F)`[2], anova(polynom5(),polynom4())$`Pr(>F)`[2], anova(polynom6(),polynom5())$`Pr(>F)`[2], anova(polynom7(),polynom6())$`Pr(>F)`[2], anova(polynom8(),polynom7())$`Pr(>F)`[2], anova(polynom9(),polynom8())$`Pr(>F)`[2]) <= 0.05)) + 1
      } else {
        min(which(c(anova(polynom2(),polynom1())$`Pr(>F)`[2], anova(polynom3(),polynom2())$`Pr(>F)`[2], anova(polynom4(),polynom3())$`Pr(>F)`[2], anova(polynom5(),polynom4())$`Pr(>F)`[2], anova(polynom6(),polynom5())$`Pr(>F)`[2], anova(polynom7(),polynom6())$`Pr(>F)`[2], anova(polynom8(),polynom7())$`Pr(>F)`[2], anova(polynom9(),polynom8())$`Pr(>F)`[2]) > 0.05))
      }
    } else {
      #Mark corrected best_polynom_number
      vector_best_polynoms_each<- NA
      for (i in 4:ncol(predicted_markers_datatable())) {
        vector_best_polynoms_each[i-3]<- if(min(which(c(anova(polynom_corrected(2,i),polynom_corrected(1,i))$`Pr(>F)`[2], anova(polynom_corrected(3,i),polynom_corrected(2,i))$`Pr(>F)`[2], anova(polynom_corrected(4,i),polynom_corrected(3,i))$`Pr(>F)`[2], anova(polynom_corrected(5,i),polynom_corrected(4,i))$`Pr(>F)`[2], 
                                                        anova(polynom_corrected(6,i),polynom_corrected(5,i))$`Pr(>F)`[2], anova(polynom_corrected(7,i),polynom_corrected(6,i))$`Pr(>F)`[2], anova(polynom_corrected(8,i),polynom_corrected(7,i))$`Pr(>F)`[2], anova(polynom_corrected(9,i),polynom_corrected(8,i))$`Pr(>F)`[2]) > 0.05)) == Inf) {
          max(which(c(anova(polynom_corrected(2,i),polynom_corrected(1,i))$`Pr(>F)`[2], anova(polynom_corrected(3,i),polynom_corrected(2,i))$`Pr(>F)`[2], anova(polynom_corrected(4,i),polynom_corrected(3,i))$`Pr(>F)`[2], anova(polynom_corrected(5,i),polynom_corrected(4,i))$`Pr(>F)`[2], anova(polynom_corrected(6,i),polynom_corrected(5,i))$`Pr(>F)`[2], 
                      anova(polynom_corrected(7,i),polynom_corrected(6,i))$`Pr(>F)`[2], anova(polynom_corrected(8,i),polynom_corrected(7,i))$`Pr(>F)`[2], anova(polynom_corrected(9,i),polynom_corrected(8,i))$`Pr(>F)`[2]) <= 0.05)) + 1
        } else {
          min(which(c(anova(polynom_corrected(2,i),polynom_corrected(1,i))$`Pr(>F)`[2], anova(polynom_corrected(3,i),polynom_corrected(2,i))$`Pr(>F)`[2], anova(polynom_corrected(4,i),polynom_corrected(3,i))$`Pr(>F)`[2], anova(polynom_corrected(5,i),polynom_corrected(4,i))$`Pr(>F)`[2], anova(polynom_corrected(6,i),polynom_corrected(5,i))$`Pr(>F)`[2], 
                      anova(polynom_corrected(7,i),polynom_corrected(6,i))$`Pr(>F)`[2], anova(polynom_corrected(8,i),polynom_corrected(7,i))$`Pr(>F)`[2], anova(polynom_corrected(9,i),polynom_corrected(8,i))$`Pr(>F)`[2]) > 0.05))
        }
      }
      vector_best_polynoms_each<- mode_calc(vector_best_polynoms_each)
      return(vector_best_polynoms_each)
    }
    
  })
  
  #Recommending best polynom + anova value with the previous
  output$polynom_recommendation<- renderText({
    paste("For the best fit, you should choose number ","<b>", best_polynom_number(),"</b>", ".", sep="")
  })
  
  dataframe_anova_polynom_test<- reactive({
    if (values$marker_correction == FALSE) {
      dataframe_anova_poly_test<- data.frame(c("pol. 2 vs pol. 1","pol. 3 vs pol. 2","pol. 4 vs pol. 3","pol. 5 vs pol. 4","pol. 6 vs pol. 5","pol. 7 vs pol. 6","pol. 8 vs pol. 7","pol. 9 vs pol. 8"),
                                             c(anova(polynom1(),polynom2())$`Pr(>F)`[2], anova(polynom2(),polynom3())$`Pr(>F)`[2], anova(polynom3(),polynom4())$`Pr(>F)`[2], anova(polynom4(),polynom5())$`Pr(>F)`[2], 
                                               anova(polynom5(),polynom6())$`Pr(>F)`[2], anova(polynom6(),polynom7())$`Pr(>F)`[2], anova(polynom7(),polynom8())$`Pr(>F)`[2], anova(polynom8(),polynom9())$`Pr(>F)`[2]))
      colnames(dataframe_anova_poly_test)<- c("model vs model", "p-value")
      return(dataframe_anova_poly_test)
    } else {
      #Mark corrected anova table for the best polynom based on the slider
      dataframe_anova_poly_test<- data.frame(c("pol. 2 vs pol. 1","pol. 3 vs pol. 2","pol. 4 vs pol. 3","pol. 5 vs pol. 4","pol. 6 vs pol. 5","pol. 7 vs pol. 6","pol. 8 vs pol. 7","pol. 9 vs pol. 8"),
                                             c(anova(polynom_corrected(1,input$slider_sample_pxmw_tab+3),polynom_corrected(2,input$slider_sample_pxmw_tab+3))$`Pr(>F)`[2], anova(polynom_corrected(2,input$slider_sample_pxmw_tab+3),polynom_corrected(3,input$slider_sample_pxmw_tab+3))$`Pr(>F)`[2], 
                                               anova(polynom_corrected(3,input$slider_sample_pxmw_tab+3),polynom_corrected(4,input$slider_sample_pxmw_tab+3))$`Pr(>F)`[2], anova(polynom_corrected(4,input$slider_sample_pxmw_tab+3),polynom_corrected(5,input$slider_sample_pxmw_tab+3))$`Pr(>F)`[2], 
                                               anova(polynom_corrected(5,input$slider_sample_pxmw_tab+3),polynom_corrected(6,input$slider_sample_pxmw_tab+3))$`Pr(>F)`[2], anova(polynom_corrected(6,input$slider_sample_pxmw_tab+3),polynom_corrected(7,input$slider_sample_pxmw_tab+3))$`Pr(>F)`[2], 
                                               anova(polynom_corrected(7,input$slider_sample_pxmw_tab+3),polynom_corrected(8,input$slider_sample_pxmw_tab+3))$`Pr(>F)`[2], anova(polynom_corrected(8,input$slider_sample_pxmw_tab+3),polynom_corrected(9,input$slider_sample_pxmw_tab+3))$`Pr(>F)`[2]))
      colnames(dataframe_anova_poly_test)<- c("model vs model", "p-value")
      return(dataframe_anova_poly_test)
    }
  })
  
  output$dataframe_anova_test<- renderTable(dataframe_anova_polynom_test(), align = "c")
  
  #Choices for polynom selection
  output$polynom_number_ui<- renderUI({
    selectInput("polynom_number", label = "Select polynomial degree", choices = 1:9, selected = best_polynom_number(), width = 187)
  })
  
  #Calculating polynom with input selected number
  best_polynom<- reactive({polynom(input$polynom_number)})
  
  
  #Making vector with coefficient values of the best polynom
  coefficient_values<- reactive({
    coefficient_val<-as.numeric(best_polynom()$coefficients)
    length(best_polynom()$coefficients)
    coefficient_val<- c(coefficient_val, rep(0,10-length(best_polynom()$coefficients)))
    return(coefficient_val)
  })
  
  #Function for recalculating pixel to bp
  polynom_equation<- function(x) {
    coefficient_values()[10]*x^9 + coefficient_values()[9]*x^8 + coefficient_values()[8]*x^7 + coefficient_values()[7]*x^6 + coefficient_values()[6]*x^5 + coefficient_values()[5]*x^4 + coefficient_values()[4]*x^3 + coefficient_values()[3]*x^2 + coefficient_values()[2]*x + coefficient_values()[1]
  }
  
  #Choose polynomial button to take user to the Length calc. tab
  output$ui.choose <- renderUI({
    if (is.na(best_polynom())) return()
    actionButton("choose", "Choose", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  })
  
  observeEvent(input$choose, {
    
    #Two lines to accept if analysis is aborted and markertable is changed
    restart_intensity<- maxima_dataframe_final()
    dataframe$dataframe(restart_intensity)
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "length_calc")
    shinyjs::enable(selector = ".navbar-nav a[data-value=length_calc]")
    shinyjs::disable(selector = ".navbar-nav a[data-value=marker]")
    shinyjs::disable(selector = ".navbar-nav a[data-value=pxmw_ratio]")
    shinyjs::disable(selector = ".navbar-nav a[data-value=marker_correction]")
  })
  
  #Trimming of the basic dataframe with all values to the last selected band of the marker (basic dataframe was cut of the row with pixel positions)
  maxima_dataframe_trimmed<- reactive({
    if (values$marker_correction == FALSE) {
    filedata_trimmed<- filedata()[-1,]
    input_markertable<- as.data.frame(hot_to_r(input$markertable))
    input_markertable<- input_markertable[!is.na(input_markertable[[2]]) & !is.na(input_markertable[[1]]),]
    last_peak_maximum<- maxima_dataframe()[1][maxima_dataframe()[4] == max(input_markertable[[2]])]
    maxima_dataframe_trimm<- filedata_trimmed[filedata_trimmed$pixel <= last_peak_maximum,]
    return(maxima_dataframe_trimm)
    } else {
      #!!! After marker correction!!! trimming of the basic dataframe with all values to the last selected band of the marker (the lowest of the two used) (basic dataframe was cut of the row with pixel positions)
      filedata_trimmed<- filedata()[-1,]
      #take the last row from the dataframe with x-values corresponding to bands
      last_peak_maximum<- as.numeric(predicted_markers_datatable()[nrow(predicted_markers_datatable()),])
      last_peak_maximum<- min(last_peak_maximum[-3])
      maxima_dataframe_trimm<- filedata_trimmed[filedata_trimmed$pixel <= last_peak_maximum,]
      maxima_dataframe_trimm<- maxima_dataframe_trimm[,-(input$sample_slider_for_second_marker+2)]
      colnames(maxima_dataframe_trimm)<- c("pixel", "marker", 1:(ncol(maxima_dataframe_trimm)-2))
      return(maxima_dataframe_trimm)
    }
  })
  
  #Dataframe with recalculated pixel to bp
  dataframe_bp_recalculated<- reactive({data.frame(maxima_dataframe_trimmed()[1],polynom_equation(maxima_dataframe_trimmed()[1]))
  })
  
  #!!!After mark correction!!! Dataframe with recalculated pixel to bp 
  dataframe_bp_recalculated_after_mark_correction<- reactive({
    df_bp_ramc<- data.frame(maxima_dataframe_trimmed()[1])
    for (i in 4:ncol(predicted_markers_datatable())) {
      #Calculating polynom with input selected number
      best_polynom<- polynom_corrected(input$polynom_number,i)
      #Making vector with coefficient values of the best polynom
      coefficient_val<-as.numeric(best_polynom$coefficients)
      coefficient_val<- c(coefficient_val, rep(0,10-length(best_polynom$coefficients)))
      coefficient_vals<- coefficient_val
      px_to_bp_sample<- coefficient_vals[10]*(maxima_dataframe_trimmed()[1])^9 + coefficient_vals[9]*(maxima_dataframe_trimmed()[1])^8 + coefficient_vals[8]*(maxima_dataframe_trimmed()[1])^7 + coefficient_vals[7]*(maxima_dataframe_trimmed()[1])^6 + coefficient_vals[6]*(maxima_dataframe_trimmed()[1])^5 + coefficient_vals[5]*(maxima_dataframe_trimmed()[1])^4 + 
        coefficient_vals[4]*(maxima_dataframe_trimmed()[1])^3 + coefficient_vals[3]*(maxima_dataframe_trimmed()[1])^2 + coefficient_vals[2]*(maxima_dataframe_trimmed()[1]) + coefficient_vals[1]
      colnames(px_to_bp_sample)<- as.character(i-3)
      df_bp_ramc<- cbind(df_bp_ramc, px_to_bp_sample)
    }
    return(df_bp_ramc)
  })
  
  #ggplot2 function for showing polynomial regression (used later as output)
  polynom_graph<- reactive({
    if (values$marker_correction == FALSE) {
      ggplot2::ggplot(data = dataframe_polynom(), aes(x = as.numeric(dataframe_polynom()[[1]]), y = as.numeric(dataframe_polynom()[[2]]))) +
        geom_point(data = dataframe_polynom(), colour= "#4493d3", size= 3) +
        geom_line(data = dataframe_bp_recalculated(), aes(x = as.numeric(dataframe_bp_recalculated()[[1]]), y = as.numeric(dataframe_bp_recalculated()[[2]])), colour = "grey") +
        ylab("bp") +
        xlab("pixel") +
        theme_light() +
        theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), axis.line = element_line(colour = "grey"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major=element_blank(), panel.border = element_blank())
    } else {
      #plot for polynomial regressions after marker correction
      ggplot2::ggplot(data = predicted_markers_datatable(), aes(x = as.numeric(predicted_markers_datatable()[[input$slider_sample_pxmw_tab+3]]), y = as.numeric(predicted_markers_datatable()[[3]]))) +
        geom_point(data = predicted_markers_datatable(), colour= "#4493d3", size= 3) +
        geom_line(data = dataframe_bp_recalculated_after_mark_correction(), aes(x = as.numeric(dataframe_bp_recalculated_after_mark_correction()[[1]]), y = as.numeric(dataframe_bp_recalculated_after_mark_correction()[[input$slider_sample_pxmw_tab+1]])), colour = "grey") +
        ylab("bp") +
        xlab("pixel") +
        theme_light() +
        theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), axis.line = element_line(colour = "grey"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major=element_blank(), panel.border = element_blank())
    }
  })
  
  output$polynom_plot<- renderPlot(polynom_graph(), height = 462)
  
  #Final dataframe with basic data used for taking x-coordinates from a graph
  maxima_dataframe_final<- reactive({
    if (values$marker_correction == FALSE) {
      maxima_dataframe_fin<- data.frame(maxima_dataframe_trimmed()[3:length(colnames(maxima_dataframe_trimmed()))], dataframe_bp_recalculated()[2], dataframe_bp_recalculated()[1])
      colnames(maxima_dataframe_fin)<- c(1:(length(colnames(maxima_dataframe_trimmed()))-2),"bp","pixel")
      maxima_dataframe_fin[] <- lapply(maxima_dataframe_fin, function(x) as.numeric(as.character(x)))
      return(maxima_dataframe_fin)
    } else {
      maxima_dataframe_fin<- data.frame(maxima_dataframe_trimmed()[3:length(colnames(maxima_dataframe_trimmed()))], dataframe_bp_recalculated_after_mark_correction()[2], dataframe_bp_recalculated_after_mark_correction()[1])
      colnames(maxima_dataframe_fin)<- c(1:(length(colnames(maxima_dataframe_trimmed()))-2),"bp","pixel")
      maxima_dataframe_fin[] <- lapply(maxima_dataframe_fin, function(x) as.numeric(as.character(x)))
      return(maxima_dataframe_fin)
    }
  })
  
  #Slider for number of samples
  output$slider <- renderUI({
    sliderInput("sample_slider", "Which sample you want to plot?", min=1, max=length(colnames(maxima_dataframe_trimmed()))-2, value=1, step = 1, ticks = TRUE, width = 217)
  })
  
  #> button for sliding on slider for number of samples forward
  observeEvent(input$forward_button_slider, {
    if(input$sample_slider < length(colnames(maxima_dataframe_trimmed()))-2) {
      updateSliderInput(session, "sample_slider", value=input$sample_slider+1)
    }
  })
  
  #< button for sliding on slider for number of samples backward
  observeEvent(input$backward_button_slider, {
    if(input$sample_slider > 1) {
      updateSliderInput(session, "sample_slider", value=input$sample_slider-1)
    }
  })
  
  #Ticks for plot sample_coordinates_graph
  ticks_sample_coordinates_graph<- reactive({
    seq(from=min(maxima_dataframe_final()$pixel), to=max(maxima_dataframe_final()$pixel), by=max(maxima_dataframe_final()$pixel)*0.1)
  })
  
  #Ticks for plot sample_coordinates_graph_bp_scale
  ticks_sample_coordinates_graph_bp_scale<- reactive({
    seq(from=min(maxima_dataframe_final()$pixel), to=max(maxima_dataframe_final()$pixel), by=max(maxima_dataframe_final()$pixel)*0.1)
  })
  
  #Background correction calculation based on previously selected areas by later introduced first and last peak of telomeres

    #Left side interval until automatically selected left border
  interval_bg_left_side<- reactive({
    interval_bg_left<- as.numeric(maxima_dataframe_final()[[input$sample_slider]])
    interval_bg_left_fin<- interval_bg_left[0:point_value_left()$x_value]
    return(interval_bg_left_fin)
  })
    #Right side interval from automatically selected right border until the end of the TRF scan
  interval_bg_right_side<- reactive({
    interval_bg_right<- as.numeric(maxima_dataframe_final()[[input$sample_slider]])
    interval_bg_right_fin<- interval_bg_right[
      point_value_right()$x_value
      :max(maxima_dataframe_final()$pixel)]
    return(interval_bg_right_fin)
  })
  
    #Rewriting the dataframe with values with the BG corrected one by clicking on the BG correction button
  dataframe<- reactiveValues()
  dataframe$dataframe<- reactiveVal({0})

  observeEvent(input$background_correction, {
      #Create dataframe from maxima_pixel_final dataframe to reactive value
    test_if_null<- dataframe$dataframe()
    maxima_df_final<- if(!is.data.frame(test_if_null)) {
      maxima_dataframe_final()
      } else {
      test_if_null
    }
      #Calculation of the left background minimum by finding local minimum, taking out every local minium that 
      #is higher than median and then calculating point in the middle of the linear regression that 
      #goes through the rest of the points (needs at least 3 points, otherwise local minima used)
    left_side_minima <- find_peaks(-interval_bg_left_side(), m = 2)
    if (length(left_side_minima) >= 3) {
      dataframe_left<- data.frame(left_side_minima,interval_bg_left_side()[left_side_minima])
      colnames(dataframe_left)<- c("x","y")
      points_left<- dataframe_left$x[dataframe_left$y <= median(dataframe_left$y)]
      regression_left<- lm(interval_bg_left_side()[points_left] ~ points_left)
      x_value_left<- (((max(points_left)-min(points_left))/2)+min(points_left))
      y_value_left<- regression_left$coefficients[[1]]+regression_left$coefficients[[2]]*(((max(points_left)-min(points_left))/2)+min(points_left))
    } else if (length(left_side_minima) >= 1) {
      dataframe_left<- data.frame(which.min(interval_bg_left_side()),interval_bg_left_side()[which.min(interval_bg_left_side())])
      colnames(dataframe_left)<- c("x","y")
      x_value_left<- dataframe_left$x
      y_value_left<- dataframe_left$y
    } else {
      dataframe_left<- data.frame(point_value_left()$x_value,point_value_left()$y_value)
      colnames(dataframe_left)<- c("x","y")
      x_value_left<- dataframe_left$x
      y_value_left<- dataframe_left$y
    }
    
  
      #Calculation of the right background minimum by finding local minimum, taking out every that 
      #is higher than median and then calculating point in the middle of the linear regression that 
      #goes through the rest of the points (needs at least 3 points, otherwise local minima used)
    right_side_minima <- find_peaks(-interval_bg_right_side(), m = 2)
    if (length(right_side_minima) >= 3) {
      dataframe_right<- data.frame(right_side_minima,interval_bg_right_side()[right_side_minima])
      colnames(dataframe_right)<- c("x","y")
      points_right<- dataframe_right$x[dataframe_right$y <= median(dataframe_right$y)]
      regression_right<- lm(interval_bg_right_side()[points_right] ~ points_right)
      x_value_right<- ((((max(points_right)-min(points_right))/2)+min(points_right))+ point_value_right()$x_value)
      y_value_right<- regression_right$coefficients[[1]]+regression_right$coefficients[[2]]*(((max(points_right)-min(points_right))/2)+min(points_right))
    } else if (length(right_side_minima) >= 1) {
      dataframe_right<- data.frame(which.min(interval_bg_right_side()),interval_bg_right_side()[which.min(interval_bg_right_side())])
      colnames(dataframe_right)<- c("x","y")
      x_value_right<- dataframe_right$x
      y_value_right<- dataframe_right$y
    } else {
      dataframe_right<- data.frame(point_value_right()$x_value,point_value_right()$y_value)
      colnames(dataframe_right)<- c("x","y")
      x_value_right<- dataframe_right$x
      y_value_right<- dataframe_right$y
    }
    
      
      #x and y corrdinates + afterward calculation of the background based on the linear model through those points
    x_values<- c(x_value_left,x_value_right)
    y_values<- c(y_value_left,y_value_right)
    regression_total<- lm(y_values ~ x_values)
    vector_for_recalculation<- 0
    for (x in 1:(max(maxima_dataframe_final()$pixel)+1)) {
      vector_for_recalculation[x]<- regression_total$coefficients[[1]]+regression_total$coefficients[[2]]*x
    }
    vector_for_recalculation[vector_for_recalculation < 0]<- 0
    vector_for_recalculation_final<- as.numeric(maxima_dataframe_final()[[input$sample_slider]]) - vector_for_recalculation
    vector_for_recalculation_final[vector_for_recalculation_final < 0]<- 0
    maxima_df_final[[input$sample_slider]]<- vector_for_recalculation_final
    
    dataframe$dataframe(maxima_df_final)
  })
  
  
  #Restart button
  observeEvent(input$restart, {
    restart_intensity<- maxima_dataframe_final()[[input$sample_slider]]
    corrected_dataframe<- maxima_dataframe_corrected()
    corrected_dataframe[[input$sample_slider]]<- restart_intensity
    dataframe$dataframe(corrected_dataframe)
  })
 
  #Dataframe after BG correction
  maxima_dataframe_corrected<- reactive({
    test_if_null<- dataframe$dataframe()
    maxima_df_final<- if(!is.data.frame(test_if_null)) {
      maxima_dataframe_final()
    } else {
      test_if_null
    }
    return(maxima_df_final)
  })
  
  #Button for BG correction
  output$background_correction_button <- renderUI({
    if (sum(maxima_dataframe_final()[[input$sample_slider]]) == sum(maxima_dataframe_corrected()[[input$sample_slider]])) {
      actionButton("background_correction", "BG Correction", width = 115)
    } else {
      actionButton("nothing", "BG Correction", width = 115, style="background-color: #E5E5E5")
    }
    
  })
  
  #Slider for threshold
  output$threshold <- renderUI({
    sliderInput("threshold_slider", "Threshold value (%) (default 28%):", min=1, max=100, value=28, step = 1, ticks = TRUE, width = 281)
  })
  
  #Making input$sample_slider a value
  observeEvent(input$threshold_slider, {
    values$threshold_slider<- reactive({input$threshold_slider})
  })   
  
  #point calculation for left side based on minimum baseline value based on slider
  point_value_left<- reactive({
    smooth_curve<- loess(as.numeric(maxima_dataframe_corrected()[[input$sample_slider]])~maxima_dataframe_corrected()$pixel, span = 0.05)
    threshold_val_max<- predict(smooth_curve)
    threshold_val_max<- threshold_val_max[
      if(!is.na(limits_left_peak()[1])){
        limits_left_peak()[1]:limits_left_peak()[3]
      } else {
        0:max(maxima_dataframe_corrected()$pixel)
      }
      ]
    threshold_val_all<- predict(smooth_curve)
    threshold_val<- threshold_val_all[0:
                                        (if(!is.na(limits_left_peak()[1])){
                                          limits_left_peak()[1] + which.max(threshold_val_max)
                                        } else {
                                          which.max(threshold_val_max)
                                        })]
    threshold_val_fin<- (max(threshold_val_max) - min(threshold_val))*(as.numeric(values$threshold_slider())/100) + min(threshold_val)
    x_values<- rev(threshold_val)-threshold_val_fin
    x_pre_value<- which(x_values <= 0)[1]
    x_value<- length(threshold_val) - x_pre_value
    y_value<- threshold_val[if (x_value <= 0) {
      1
    }else {x_value}]
    point_left_dataframe<- data.frame(x_value,y_value)
    return(point_left_dataframe)
  })
  
  #point calculation for right side based on minimum baseline value based on slider
  point_value_right<- reactive({
    smooth_curve<- loess(as.numeric(maxima_dataframe_corrected()[[input$sample_slider]])~maxima_dataframe_corrected()$pixel, span = 0.05)
    threshold_val_max<- predict(smooth_curve)
    threshold_val_max<- threshold_val_max[
      if(!is.na(limits_right_peak()[1])){
        limits_right_peak()[1]:limits_right_peak()[3]
      } else {
        0:max(maxima_dataframe_corrected()$pixel)
      }
      ]
    threshold_val_all<- predict(smooth_curve)
    threshold_val<- threshold_val_all[
      ((if(!is.na(limits_right_peak()[1])){
        limits_right_peak()[1]
      } else {
        0
      }) + which.max(threshold_val_max))
      :max(maxima_dataframe_corrected()$pixel)]
    threshold_val_fin<- (max(threshold_val_max) -  min(threshold_val))*(as.numeric(values$threshold_slider())/100) + min(threshold_val)
    x_values<- threshold_val-threshold_val_fin
    x_value<- if (sum(x_values <= 0) == 0) {
      (if(!is.na(limits_right_peak()[1])){
        limits_right_peak()[1] + which.max(threshold_val_max)
      } else {
        0 + which.max(threshold_val_max)
      })
    } else {which(x_values <= 0)[1] + 
      (if(!is.na(limits_right_peak()[1])){
      limits_right_peak()[1] + which.max(threshold_val_max)
    } else {
      0 + which.max(threshold_val_max)
    })}
    y_value<- threshold_val_all[x_value]
    point_right_dataframe<- data.frame(x_value,y_value)
    return(point_right_dataframe)
  })
  
  #function for showing a graph of a sample based on a slider input sample number
  sample_coordinates_graph<- reactive({
    if (values$marker_correction == FALSE) {
      smooth_curve<- loess(as.numeric(maxima_dataframe_corrected()[[input$sample_slider]])~maxima_dataframe_corrected()$pixel, span = 0.05)
      par(family = "Roboto")
      plot(x= as.numeric(maxima_dataframe_corrected()$pixel), y= as.numeric(maxima_dataframe_corrected()[[input$sample_slider]]), xaxt= "n", type= "l", xaxs= "i", ylab= "", xlab="pixel", col="grey", fg="grey", mgp=c(1.35, 0.2, 0), col.axis="gray30", cex.axis=0.75, cex.lab=1.1, las=1, tck="-0.005")
      rect(xleft = point_value_left()$x_value,ybottom = -1000, xright = point_value_right()$x_value, ytop = max(maxima_dataframe_corrected()[[input$sample_slider]])*2,
           col= rgb(1,0,0,alpha=0.02), border="#fa3939", lty=NULL, lwd=par("lwd"), xpd=FALSE)
      lines(predict(smooth_curve), col="#4493d3",lwd=2)
      title(ylab = "intensity", mgp = c(2.35, 0.2, 0), cex.lab=1.1)
      axis(3, at = ticks_sample_coordinates_graph_bp_scale(), labels = round(maxima_dataframe_corrected()$bp[(ticks_sample_coordinates_graph_bp_scale()+1)],digits = 2), col="grey", tck="-0.005", col.axis="gray30", cex.axis=0.75, mgp = c(1.45, 0.3, 0))
      axis(1, at = ticks_sample_coordinates_graph(), labels = ticks_sample_coordinates_graph(), col="grey", tck="-0.005", col.axis="gray30", cex.axis=0.75, mgp = c(1.45, 0.2, 0))
      mtext(side=3,text="bp", cex=1.1, line=1.45)
      points(x = point_value_left()$x_value, y = point_value_left()$y_value, col="#ff4d4d", pch=15, cex = 1.65)
      points(x = point_value_left()$x_value, y = point_value_left()$y_value, col="white", pch=-8250, cex = 1.25)
      points(x = point_value_right()$x_value, y = point_value_right()$y_value, col="#ff4d4d", pch=15, cex = 1.65)
      points(x = point_value_right()$x_value, y = point_value_right()$y_value, col="white", pch=-8249, cex = 1.25)
    } else {
      smooth_curve<- loess(as.numeric(maxima_dataframe_corrected()[[input$sample_slider]])~maxima_dataframe_corrected()$pixel, span = 0.05)
      par(family = "Roboto")
      plot(x= as.numeric(maxima_dataframe_corrected()$pixel), y= as.numeric(maxima_dataframe_corrected()[[input$sample_slider]]), xaxt= "n", type= "l", xaxs= "i", ylab= "", xlab="pixel", col="grey", fg="grey", mgp=c(1.35, 0.2, 0), col.axis="gray30", cex.axis=0.75, cex.lab=1.1, las=1, tck="-0.005")
      rect(xleft = point_value_left()$x_value,ybottom = -1000, xright = point_value_right()$x_value, ytop = max(maxima_dataframe_corrected()[[input$sample_slider]])*2,
           col= rgb(1,0,0,alpha=0.02), border="#fa3939", lty=NULL, lwd=par("lwd"), xpd=FALSE)
      lines(predict(smooth_curve), col="#4493d3",lwd=2)
      title(ylab = "intensity", mgp = c(2.35, 0.2, 0), cex.lab=1.1)
      axis(3, at = ticks_sample_coordinates_graph_bp_scale(), labels = round(dataframe_bp_recalculated_after_mark_correction()[[input$sample_slider+1]][(ticks_sample_coordinates_graph_bp_scale()+1)],digits = 2), col="grey", tck="-0.005", col.axis="gray30", cex.axis=0.75, mgp = c(1.45, 0.3, 0))
      axis(1, at = ticks_sample_coordinates_graph(), labels = ticks_sample_coordinates_graph(), col="grey", tck="-0.005", col.axis="gray30", cex.axis=0.75, mgp = c(1.45, 0.2, 0))
      mtext(side=3,text="bp", cex=1.1, line=1.45)
      points(x = point_value_left()$x_value, y = point_value_left()$y_value, col="#ff4d4d", pch=15, cex = 1.65)
      points(x = point_value_left()$x_value, y = point_value_left()$y_value, col="white", pch=-8250, cex = 1.25)
      points(x = point_value_right()$x_value, y = point_value_right()$y_value, col="#ff4d4d", pch=15, cex = 1.65)
      points(x = point_value_right()$x_value, y = point_value_right()$y_value, col="white", pch=-8249, cex = 1.25)
    }
  })
  
  output$sample_coordinates_plot<- renderPlot(sample_coordinates_graph(), height = 474)
  
  #Making input$brush_sample_coordinates_plot a value
  values<-reactiveValues()
  observeEvent(input$brush_sample_coordinates_plot,{
    values$brush_sample_coordinates<- reactive({input$brush_sample_coordinates_plot})
  })
  
  #Making input$sample_slider a value
  observeEvent(input$sample_slider, {
    values$sample_slider<- reactive({input$sample_slider})
  })   
  
  #Saving data for the left peak after clicking action button
  values$dataframe_coordinates_left_side <- reactiveVal({data.frame(xmin=NA, xmax=NA, sample_nu=NA, counter=NA)})
  
  observeEvent(input$save_sample_coordinates_left_side, {
    if(!is.null(input$brush_sample_coordinates_plot)) {
      old_values <- values$dataframe_coordinates_left_side()
      new_values <- 
        data.frame(xmin=values$brush_sample_coordinates()$xmin, xmax=values$brush_sample_coordinates()$xmax, sample_nu=values$sample_slider(), counter=input$save_sample_coordinates_left_side)
      new_dataframe <- rbind(old_values, new_values)
      values$dataframe_coordinates_left_side(new_dataframe)
    }
  })
  
  #Select last measurements for each sample number from "memory" created above for the left peak
  coordinates_from_memory_left_side<- reactive({
    coordinates_memory<- values$dataframe_coordinates_left_side()
    require(dplyr)
    coordinates_memory %>% group_by(sample_nu) %>% filter(counter %in% max(counter)) -> coordinates_memory
    coordinates_memory_ordered<- coordinates_memory[order(coordinates_memory$sample_nu),]
    return(coordinates_memory_ordered)
  })
  
  #Creating reactive value with last measured coordinates of the left peak for selected sample
  limits_left_peak<- reactive({
  limits_left<- as.numeric(c(coordinates_from_memory_left_side()$xmin[coordinates_from_memory_left_side()$sample_nu == as.numeric(input$sample_slider)],coordinates_from_memory_left_side()$xmax[coordinates_from_memory_left_side()$sample_nu == as.numeric(input$sample_slider)]))
  return(limits_left)
  })
  
  #Values for the selection of the maxima for the left threshold calculation
  output$coordinates_left_peak_left_side<- renderText({
    paste("Left [px]: ","<b>",     
      if(!is.na(limits_left_peak()[1])){
      round(limits_left_peak()[1], digits = 0)
    } else {
      0
    }
    , sep="")
  })
  output$coordinates_left_peak_right_side<- renderText({
    paste("Right [px]: ","<b>",
    if(!is.na(limits_left_peak()[1])){
      round(limits_left_peak()[3], digits = 0)
    } else {
      max(maxima_dataframe_corrected()$pixel)
    }
    , sep="")
  })
  
  #Saving data for the right peak after clicking action button
  values$dataframe_coordinates_right_side <- reactiveVal({data.frame(xmin=NA, xmax=NA, sample_nu=NA, counter=NA)})
  
  observeEvent(input$save_sample_coordinates_right_side, {
    if(!is.null(input$brush_sample_coordinates_plot)) {
      old_values <- values$dataframe_coordinates_right_side()
      new_values <- 
        data.frame(xmin=values$brush_sample_coordinates()$xmin, xmax=values$brush_sample_coordinates()$xmax, sample_nu=values$sample_slider(), counter=input$save_sample_coordinates_right_side)
      new_dataframe <- rbind(old_values, new_values)
      values$dataframe_coordinates_right_side(new_dataframe)
    }
  })
  
  #Select last measurements for each sample number from "memory" created above for the right peak
  coordinates_from_memory_right_side<- reactive({
    coordinates_memory<- values$dataframe_coordinates_right_side()
    require(dplyr)
    coordinates_memory %>% group_by(sample_nu) %>% filter(counter %in% max(counter)) -> coordinates_memory
    coordinates_memory_ordered<- coordinates_memory[order(coordinates_memory$sample_nu),]
    return(coordinates_memory_ordered)
  })
  
  #Creating reactive value with last measured coordinates of the left peak for selected sample
  limits_right_peak<- reactive({
    limits_right<- as.numeric(c(coordinates_from_memory_right_side()$xmin[coordinates_from_memory_right_side()$sample_nu == as.numeric(input$sample_slider)],coordinates_from_memory_right_side()$xmax[coordinates_from_memory_right_side()$sample_nu == as.numeric(input$sample_slider)]))
    return(limits_right)
  })
  
  #Values for the selection of the maxima for the right threshold calculation
  output$coordinates_right_peak_left_side<- renderText({
    paste("Left [px]: ","<b>",     
          if(!is.na(limits_right_peak()[1])){
            round(limits_right_peak()[1], digits = 0)
          } else {
            0
          }
          , sep="")
    })
  output$coordinates_right_peak_right_side<- renderText({
    paste("Right [px]: ","<b>",
          if(!is.na(limits_right_peak()[1])){
            round(limits_right_peak()[3], digits = 0)
          } else {
            max(maxima_dataframe_corrected()$pixel)
          }
          , sep="")
  })
  
  #Saving data after clicking Manual action button
  values$dataframe_coordinates <- reactiveVal({data.frame(xmin=NA, xmax=NA, sample_nu=NA, counter=NA, auto_save=NA, background_correction=NA, threshold_value=NA)})
  observeEvent(input$save_sample_coordinates, {
    if(!is.null(input$brush_sample_coordinates_plot)) {
      remove_sample_dataframe_all<- values$dataframe_coordinates()
      remove_sample_dataframe <- remove_sample_dataframe_all[remove_sample_dataframe_all$sample_nu != input$sample_slider,]
      values$dataframe_coordinates(remove_sample_dataframe)
      shinyjs::enable(selector = ".navbar-nav a[data-value=graphic_result]")
      old_values <- values$dataframe_coordinates()
      new_values <- 
        data.frame(xmin=as.numeric(values$brush_sample_coordinates()$xmin), xmax=as.numeric(values$brush_sample_coordinates()$xmax), sample_nu=values$sample_slider(), counter=input$save_sample_coordinates, auto_save=(0 == 1), background_correction=(sum(maxima_dataframe_final()[[input$sample_slider]]) != sum(maxima_dataframe_corrected()[[input$sample_slider]])), threshold_value=input$threshold_slider)
      new_dataframe <- rbind(old_values, new_values)
      values$dataframe_coordinates(new_dataframe)
      #Later code is for updating name of samples in Graphic result while clicking on the button
      col_numbers<- colnames(cumsum_maxima_datatable())
      col_numbers<- as.numeric(col_numbers)
      grouptable$column_numbers(col_numbers)
      NAs<- as.character(rep("",length(col_numbers)))
      grouptable$NAs(NAs)
      shinyjs::enable("control_group_name")
      values$statistical_value<- FALSE
      values$statistical_counter<- 0
    }
  })
  
  #Saving data after clicking Automatic action button
  observeEvent(input$save_sample_coordinates_auto, {
    shinyjs::enable(selector = ".navbar-nav a[data-value=graphic_result]")
    remove_sample_dataframe_all<- values$dataframe_coordinates()
    remove_sample_dataframe <- remove_sample_dataframe_all[remove_sample_dataframe_all$sample_nu != input$sample_slider,]
    values$dataframe_coordinates(remove_sample_dataframe)
      old_values <- values$dataframe_coordinates()
      new_values <- 
        data.frame(xmin=as.numeric(point_value_left()$x_value), xmax=as.numeric(point_value_right()$x_value), sample_nu=values$sample_slider(), counter=input$save_sample_coordinates, auto_save=(1 == 1), background_correction=(sum(maxima_dataframe_final()[[input$sample_slider]]) != sum(maxima_dataframe_corrected()[[input$sample_slider]])), threshold_value=input$threshold_slider)
      new_dataframe <- rbind(old_values, new_values)
      values$dataframe_coordinates(new_dataframe)
      #Later code is for updating name of samples in Graphic result while clicking on the button
      col_numbers<- colnames(cumsum_maxima_datatable())
      col_numbers<- as.numeric(col_numbers)
      grouptable$column_numbers(col_numbers)
      NAs<- as.character(rep("",length(col_numbers)))
      grouptable$NAs(NAs)
      shinyjs::enable("control_group_name")
      values$statistical_value<- FALSE
      values$statistical_counter<- 0
  })
  
  #Information if the data was already saved and if automatically or manually
  output$info_sample<- renderText({
    if (as.character(sample_datatable()[[4]][sample_datatable()[[1]] == input$sample_slider]) == "TRUE") {
      paste("Sample was saved using ",'<span style=\"color:#ff4d4d\">',"<b>", "AUTO SAVE","</b>",'</span>',".", sep="")
    }
    else if (as.character(sample_datatable()[[4]][sample_datatable()[[1]] == input$sample_slider]) == "FALSE") {
      paste("Sample was saved using ",'<span style=\"color:#337ab7\">',"<b>", "MANUAL SAVE","</b>",'</span>',".", sep="")
    }
  })
  
  #Restart current sample
  observeEvent(input$restart_sample_coordinates, {
    remove_sample_dataframe_all<- values$dataframe_coordinates()
    remove_sample_dataframe <- remove_sample_dataframe_all[remove_sample_dataframe_all$sample_nu != input$sample_slider,]
    values$dataframe_coordinates(remove_sample_dataframe)
    
    dataframe_coordinates_left_side_all<- values$dataframe_coordinates_left_side()
    dataframe_coordinates_left_side<- dataframe_coordinates_left_side_all[dataframe_coordinates_left_side_all$sample_nu != input$sample_slider,]
    values$dataframe_coordinates_left_side(dataframe_coordinates_left_side)
    
    dataframe_coordinates_right_side_all<- values$dataframe_coordinates_right_side()
    dataframe_coordinates_right_side<- dataframe_coordinates_right_side_all[dataframe_coordinates_right_side_all$sample_nu != input$sample_slider,]
    values$dataframe_coordinates_right_side(dataframe_coordinates_right_side)
    shinyjs::enable("control_group_name")
    values$statistical_value<- FALSE
    values$statistical_counter<- 0
    
    if (nrow(values$dataframe_coordinates()) == 1) {
      shinyjs::disable(selector = ".navbar-nav a[data-value=graphic_result]")
    }
    
    if (nrow(values$dataframe_coordinates()) > 1) {
    #Later code is for updating name of samples in Graphic result while clicking on the button
    col_numbers<- colnames(cumsum_maxima_datatable())
    col_numbers<- as.numeric(col_numbers)
    grouptable$column_numbers(col_numbers)
    NAs<- as.character(rep("",length(col_numbers)))
    grouptable$NAs(NAs)
    }
  })
  
  #Abort analysis
  observeEvent(input$restart_all_coordinates, {
    shinyalert_abort_confirm <- function(value) {
      if (value == TRUE) {
      remove_sample_dataframe_all<- values$dataframe_coordinates()
      remove_sample_dataframe <- data.frame(xmin=NA, xmax=NA, sample_nu=NA, counter=NA, auto_save=NA, background_correction=NA, threshold_value=NA)
      values$dataframe_coordinates(remove_sample_dataframe)
      
      dataframe_coordinates_left_side<- values$dataframe_coordinates_left_side()
      dataframe_coordinates_left_side<- data.frame(xmin=NA, xmax=NA, sample_nu=NA, counter=NA)
      values$dataframe_coordinates_left_side(dataframe_coordinates_left_side)
      
      dataframe_coordinates_right_side<- values$dataframe_coordinates_right_side()
      dataframe_coordinates_right_side<- data.frame(xmin=NA, xmax=NA, sample_nu=NA, counter=NA)
      values$dataframe_coordinates_right_side(dataframe_coordinates_right_side)
      
      updateTabsetPanel(session = session, inputId = "tabs", selected = "marker")
      updateSliderInput(session = session, "sample_slider", value = 1)
      shinyjs::disable(selector = ".navbar-nav a[data-value=length_calc]")
      shinyjs::enable(selector = ".navbar-nav a[data-value=marker]")
      shinyjs::disable(selector = ".navbar-nav a[data-value=pxmw_ratio]")
      shinyjs::disable(selector = ".navbar-nav a[data-value=graphic_result]")
      shinyjs::disable(selector = ".navbar-nav a[data-value=marker_correction]")
      values$marker_correction<- FALSE
      }}
    shinyalert(callbackR = shinyalert_abort_confirm, "Abort analysis", "You are about to abort the analysis, do you want to proceed?", type = "warning", closeOnEsc = TRUE, confirmButtonCol =  "#337ab7", closeOnClickOutside = TRUE, showCancelButton = TRUE, showConfirmButton = TRUE)
    
  })
  
  #Select last measurements for each sample number from "memory" created above
  coordinates_from_memory<- reactive({
    coordinates_memory<- values$dataframe_coordinates()
    require(dplyr)
    coordinates_memory %>% group_by(sample_nu) %>% filter(counter %in% max(counter)) -> coordinates_memory
    coordinates_memory_ordered<- coordinates_memory[order(coordinates_memory$sample_nu),]
    return(coordinates_memory_ordered)
  })
  
  #Ticks for plot sample_coordinates_graph
  ticks_sample_coordinates_graph_selection<- reactive({
    limits_x<- as.numeric(c(coordinates_from_memory()$xmin[coordinates_from_memory()$sample_nu == as.numeric(input$sample_slider)],coordinates_from_memory()$xmax[coordinates_from_memory()$sample_nu == as.numeric(input$sample_slider)]))
    if(!is.na(limits_x[1])){
      ticks_sample_selection<- seq(from=limits_x[1], to=limits_x[3], by=(limits_x[3]-limits_x[1])/4)
    } else {
      ticks_sample_selection<- seq(from=0, to=max(maxima_dataframe_corrected()$pixel), by=(max(maxima_dataframe_corrected()$pixel)-0)/4)
    }
    return(ticks_sample_selection)
  })
  
  #function for showing a graph of a sample based on a slider input with limits based on brush
  sample_coordinates_graph_selection<- reactive({
    if (values$marker_correction == FALSE) {
      smooth_curve<- loess(as.numeric(maxima_dataframe_corrected()[[input$sample_slider]])~maxima_dataframe_corrected()$pixel, span = 0.05)
      limits_x<- as.numeric(c(coordinates_from_memory()$xmin[coordinates_from_memory()$sample_nu == as.numeric(input$sample_slider)],coordinates_from_memory()$xmax[coordinates_from_memory()$sample_nu == as.numeric(input$sample_slider)]))
      par(family = "Roboto")
      plot(x= maxima_dataframe_corrected()$pixel, y= as.numeric(maxima_dataframe_corrected()[[input$sample_slider]]), type= "l", xaxt= "n", xaxs= "i", ylab= "", xlab="pixel", col="white", fg="grey", mgp=c(1.35, 0.2, 0), col.axis="gray30", cex.axis=0.75, cex.lab=1.1, las=1, tck="-0.005", xlim= (if(!is.na(limits_x[1])){
        limits_x[c(1,3)]
      } else {
        as.numeric(c(0,max(maxima_dataframe_corrected()$pixel)))
      }))
      lines(predict(smooth_curve), col="#4493d3",lwd=2)
      title(ylab = "intensity", mgp = c(2.35, 0.2, 0), cex.lab=1.1)
      axis(1, at = ticks_sample_coordinates_graph_selection(), labels = round(ticks_sample_coordinates_graph_selection(),digits = 0), col="grey", tck="-0.005", col.axis="gray30", cex.axis=0.75, mgp = c(1.45, 0.2, 0))
      axis(3, at = ticks_sample_coordinates_graph_selection(), labels = round(maxima_dataframe_corrected()$bp[(ticks_sample_coordinates_graph_selection()+1)],digits = 2), col="grey", tck="-0.005", col.axis="gray30", cex.axis=0.75, mgp = c(1.45, 0.3, 0))
      mtext(side=3,text="bp", cex=1.1, line=1.45)
    } else {
      smooth_curve<- loess(as.numeric(maxima_dataframe_corrected()[[input$sample_slider]])~maxima_dataframe_corrected()$pixel, span = 0.05)
      limits_x<- as.numeric(c(coordinates_from_memory()$xmin[coordinates_from_memory()$sample_nu == as.numeric(input$sample_slider)],coordinates_from_memory()$xmax[coordinates_from_memory()$sample_nu == as.numeric(input$sample_slider)]))
      par(family = "Roboto")
      plot(x= maxima_dataframe_corrected()$pixel, y= as.numeric(maxima_dataframe_corrected()[[input$sample_slider]]), type= "l", xaxt= "n", xaxs= "i", ylab= "", xlab="pixel", col="white", fg="grey", mgp=c(1.35, 0.2, 0), col.axis="gray30", cex.axis=0.75, cex.lab=1.1, las=1, tck="-0.005", xlim= (if(!is.na(limits_x[1])){
        limits_x[c(1,3)]
      } else {
        as.numeric(c(0,max(maxima_dataframe_corrected()$pixel)))
      }))
      lines(predict(smooth_curve), col="#4493d3",lwd=2)
      title(ylab = "intensity", mgp = c(2.35, 0.2, 0), cex.lab=1.1)
      axis(1, at = ticks_sample_coordinates_graph_selection(), labels = round(ticks_sample_coordinates_graph_selection(),digits = 0), col="grey", tck="-0.005", col.axis="gray30", cex.axis=0.75, mgp = c(1.45, 0.2, 0))
      axis(3, at = ticks_sample_coordinates_graph_selection(), labels = round(dataframe_bp_recalculated_after_mark_correction()[[input$sample_slider+1]][(ticks_sample_coordinates_graph_selection()+1)],digits = 2), col="grey", tck="-0.005", col.axis="gray30", cex.axis=0.75, mgp = c(1.45, 0.3, 0))
      mtext(side=3,text="bp", cex=1.1, line=1.45)
    }
  })
  
  output$sample_coordinates_plot_selection<- renderPlot(sample_coordinates_graph_selection(), height = 307)
  
  #Making datatable with maxima_dataframe_corrected where values outside of selection are changed to null values -> for the cumsum calculation later
  maxima_datatable_final_0<- reactive({
    #smoothing
    maxima_dataframe_fin_0<- maxima_dataframe_corrected()
    for (i in c(1:(length(colnames(maxima_dataframe_fin_0))-2))) {
      smooth_curve<- loess(as.numeric(maxima_dataframe_fin_0[[i]])~maxima_dataframe_fin_0$pixel, span = 0.05)
      maxima_dataframe_fin_0[[i]]<- predict(smooth_curve)
    }
    
    for (i in 1:(nrow(coordinates_from_memory())-1)) {
      maxima_dataframe_fin_0[0:round(coordinates_from_memory()$xmin[i],digits = 0),coordinates_from_memory()$sample_nu[i]]<- 0
      maxima_dataframe_fin_0[round(coordinates_from_memory()$xmax[i],digits = 0) < maxima_dataframe_fin_0$pixel,coordinates_from_memory()$sample_nu[i]]<- 0
    } 
    return(maxima_dataframe_fin_0)
  })
  
  #Recalculation to cumsums for the datatable above with only selected samples
  cumsum_maxima_datatable<- reactive({
    if (values$marker_correction == FALSE) {
      cumsum_maxima_dataframe<- maxima_datatable_final_0()[,as.numeric(coordinates_from_memory()$sample_nu[!is.na(coordinates_from_memory()$sample_nu)])]
      cumsum_maxima_dataframe <- cumsum(cumsum_maxima_dataframe/maxima_dataframe_corrected()$bp)
      return(cumsum_maxima_dataframe)
    } else {
      dataframe_bp_for_division<- dataframe_bp_recalculated_after_mark_correction()[,-1]
      cumsum_maxima_dataframe<- maxima_datatable_final_0()[,as.numeric(coordinates_from_memory()$sample_nu[!is.na(coordinates_from_memory()$sample_nu)])]
      cumsum_maxima_dataframe <- cumsum(cumsum_maxima_dataframe/dataframe_bp_for_division[,as.numeric(coordinates_from_memory()$sample_nu[!is.na(coordinates_from_memory()$sample_nu)])])
      return(cumsum_maxima_dataframe)
    }
  })
  
  #Half of the maximum values of each selected sample
  max_values_cumsum_datatable<- reactive({
    max_val_cumsum_dataframe<- if(is.data.frame(cumsum_maxima_datatable())) {
      cumsum_maxima_datatable()[nrow(cumsum_maxima_datatable()),]/2
    } else {
      max(cumsum_maxima_datatable())/2
    }
    max_val_cumsum_dataframe<- as.data.frame(max_val_cumsum_dataframe)
    return(max_val_cumsum_dataframe)
  })
 
  #Final table with values min, 0.25, median, 0.75 and max for later boxplot
  final_boxplot_datatable<- reactive({
    if (values$marker_correction == FALSE) {
      final_boxplot_dataframe<- max_values_cumsum_datatable()
      final_boxplot_dataframe[1,]<- maxima_dataframe_corrected()$bp[round(coordinates_from_memory()$xmax[1:(nrow(coordinates_from_memory())-1)], digits = 0)+1]
      final_boxplot_dataframe[5,]<- maxima_dataframe_corrected()$bp[round(coordinates_from_memory()$xmin[1:(nrow(coordinates_from_memory())-1)], digits = 0)+1]
      
      #Median calculation
      for (i in 1:ncol(max_values_cumsum_datatable())) {
        
        position_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable()[[i]] <= max_values_cumsum_datatable()[[i]]),i])
        } else {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable() <= max_values_cumsum_datatable()[[i]])])
        }
        cumsum_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          as.vector(cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1),i])
        } else {
          cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1)]
        }
        bp_as_vector<- as.vector(maxima_dataframe_corrected()$bp)
        bp_by_position<- bp_as_vector[c(position_in_bp, position_in_bp+1)]
        
        final_boxplot_dataframe[3,i]<- lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[2]]*max_values_cumsum_datatable()[[i]]+lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[1]]
        
      }
      
      #0.75 value calculation
      for (i in 1:ncol(max_values_cumsum_datatable())) {
        
        position_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable()[[i]] <= ((max_values_cumsum_datatable()[[i]]*2)*0.75)),i])
        } else {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable() <= ((max_values_cumsum_datatable()[[i]]*2)*0.75))])
        }
        cumsum_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          as.vector(cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1),i])
        } else {
          cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1)]
        }
        bp_as_vector<- as.vector(maxima_dataframe_corrected()$bp)
        bp_by_position<- bp_as_vector[c(position_in_bp, position_in_bp+1)]
        
        final_boxplot_dataframe[2,i]<- lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[2]]*((max_values_cumsum_datatable()[[i]]*2)*0.75)+lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[1]]
        
      }
      
      #0.25 value calculation
      for (i in 1:ncol(max_values_cumsum_datatable())) {
        
        position_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable()[[i]] <= ((max_values_cumsum_datatable()[[i]]*2)*0.25)),i])
        } else {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable() <= ((max_values_cumsum_datatable()[[i]]*2)*0.25))])
        }
        cumsum_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          as.vector(cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1),i])
        } else {
          cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1)]
        }
        bp_as_vector<- as.vector(maxima_dataframe_corrected()$bp)
        bp_by_position<- bp_as_vector[c(position_in_bp, position_in_bp+1)]
        
        final_boxplot_dataframe[4,i]<- lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[2]]*((max_values_cumsum_datatable()[[i]]*2)*0.25)+lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[1]]
        
      }
      rownames(final_boxplot_dataframe)<- c("min","1st quartile","median", "3rd quartile", "max")
      if (length(colnames(final_boxplot_dataframe)) == 1){
        colnames(final_boxplot_dataframe) <- c("")
      }
      return(final_boxplot_dataframe)
    } else {
      
      ###After filling Mark correction tab ###
      
      final_boxplot_dataframe<- max_values_cumsum_datatable()
      final_boxplot_dataframe[1,]<- maxima_dataframe_corrected()$bp[round(coordinates_from_memory()$xmax[1:(nrow(coordinates_from_memory())-1)], digits = 0)+1]
      final_boxplot_dataframe[5,]<- maxima_dataframe_corrected()$bp[round(coordinates_from_memory()$xmin[1:(nrow(coordinates_from_memory())-1)], digits = 0)+1]
      
      dataframe_corresponding_bp<- dataframe_bp_recalculated_after_mark_correction()[,-1]
      dataframe_corresponding_bp<- dataframe_corresponding_bp[,as.numeric(coordinates_from_memory()$sample_nu[!is.na(coordinates_from_memory()$sample_nu)])]
      for (i in 1:ncol(max_values_cumsum_datatable())) {
        if (is.data.frame(dataframe_corresponding_bp) == TRUE) {
          vector_corresponding_bp<- dataframe_corresponding_bp[[i]]
        } else {
          vector_corresponding_bp<- as.numeric(dataframe_corresponding_bp)
        }
        final_boxplot_dataframe[1,i]<- vector_corresponding_bp[round(coordinates_from_memory()$xmax[i], digits = 0)+1]
        final_boxplot_dataframe[5,i]<- vector_corresponding_bp[round(coordinates_from_memory()$xmin[i], digits = 0)+1]
      }
      
      #Median calculation
      for (i in 1:ncol(max_values_cumsum_datatable())) {
        
        position_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable()[[i]] <= max_values_cumsum_datatable()[[i]]),i])
        } else {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable() <= max_values_cumsum_datatable()[[i]])])
        }
        cumsum_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          as.vector(cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1),i])
        } else {
          cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1)]
        }
        
        if (is.data.frame(dataframe_corresponding_bp) == TRUE) {
          vector_corresponding_bp<- dataframe_corresponding_bp[[i]]
        } else {
          vector_corresponding_bp<- as.numeric(dataframe_corresponding_bp)
        }
        bp_as_vector<- vector_corresponding_bp
        bp_by_position<- bp_as_vector[c(position_in_bp, position_in_bp+1)]
        
        final_boxplot_dataframe[3,i]<- lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[2]]*max_values_cumsum_datatable()[[i]]+lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[1]]
        
      }
      
      #0.75 value calculation
      for (i in 1:ncol(max_values_cumsum_datatable())) {
        
        position_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable()[[i]] <= ((max_values_cumsum_datatable()[[i]]*2)*0.75)),i])
        } else {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable() <= ((max_values_cumsum_datatable()[[i]]*2)*0.75))])
        }
        cumsum_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          as.vector(cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1),i])
        } else {
          cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1)]
        }
        if (is.data.frame(dataframe_corresponding_bp) == TRUE) {
          vector_corresponding_bp<- dataframe_corresponding_bp[[i]]
        } else {
          vector_corresponding_bp<- as.numeric(dataframe_corresponding_bp)
        }
        bp_as_vector<- vector_corresponding_bp
        bp_by_position<- bp_as_vector[c(position_in_bp, position_in_bp+1)]
        
        final_boxplot_dataframe[2,i]<- lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[2]]*((max_values_cumsum_datatable()[[i]]*2)*0.75)+lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[1]]
        
      }
      
      #0.25 value calculation
      for (i in 1:ncol(max_values_cumsum_datatable())) {
        
        position_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable()[[i]] <= ((max_values_cumsum_datatable()[[i]]*2)*0.25)),i])
        } else {
          length(cumsum_maxima_datatable()[(cumsum_maxima_datatable() <= ((max_values_cumsum_datatable()[[i]]*2)*0.25))])
        }
        cumsum_in_bp<- if(is.data.frame(cumsum_maxima_datatable())) {
          as.vector(cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1),i])
        } else {
          cumsum_maxima_datatable()[position_in_bp:(position_in_bp+1)]
        }
        if (is.data.frame(dataframe_corresponding_bp) == TRUE) {
          vector_corresponding_bp<- dataframe_corresponding_bp[[i]]
        } else {
          vector_corresponding_bp<- as.numeric(dataframe_corresponding_bp)
        }
        bp_as_vector<- vector_corresponding_bp
        bp_by_position<- bp_as_vector[c(position_in_bp, position_in_bp+1)]
        
        final_boxplot_dataframe[4,i]<- lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[2]]*((max_values_cumsum_datatable()[[i]]*2)*0.25)+lm(bp_by_position ~ cumsum_in_bp)$`coefficients`[[1]]
        
      }
      rownames(final_boxplot_dataframe)<- c("min","1st quartile","median", "3rd quartile", "max")
      if (length(colnames(final_boxplot_dataframe)) == 1){
        colnames(final_boxplot_dataframe) <- c("")
      }
      
      return(final_boxplot_dataframe)
    }
  })
  #coordinates_from_memory() to datatable for rmarkdown report
  sample_datatable<- reactive({
    sample_dataframe<- as.data.frame(coordinates_from_memory()[,c(3,1,2,5,6,7)])
    sample_dataframe<- sample_dataframe[-nrow(sample_dataframe),]
    colnames(sample_dataframe)[1] <- "sample_number"
    return(sample_dataframe)
  })
  
  #Download rmarkdown report
  output$report <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_report_IntensityAnalyser_WALTER", ".html", sep = "")
    },
    content = function(file) {

      tempReport <- file.path(tempdir(), "report.rmd")
      file.copy("report.rmd", tempReport, overwrite = TRUE)
      
      params <- list(marker_table = markers_filled_table(), sensitivity_value = sensitivity_vals(), polynomial_value = best_polynom_number(), sample_table = sample_datatable(), boxplot_table = final_boxplot_datatable(), recalculated_table = mean_sd_n_datatable(), groups_table = datatable_groups(), mean_groups_table = datatable_mean_groups(), used_test = test_used(), significance_table = datatable_significance())
      
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  #Transpose final_boxplot_datatable for ggplot2
  final_ggplot2_boxplot_datatable_transposed<- reactive ({
    final_ggplot2_boxplot_dataframe<- data.frame(t(final_boxplot_datatable()))
    final_ggplot2_boxplot_dataframe<- as.matrix(final_ggplot2_boxplot_dataframe)
    class(final_ggplot2_boxplot_dataframe) <- "numeric"
    final_ggplot2_boxplot_dataframe<- as.data.frame(final_ggplot2_boxplot_dataframe)
    final_ggplot2_boxplot_dataframe<- cbind(colnames(final_boxplot_datatable()), final_ggplot2_boxplot_dataframe)
    colnames(final_ggplot2_boxplot_dataframe)<- c("X1","0%","25%","50%","75%","100%")
    final_ggplot2_boxplot_dataframe$X1 <- factor(final_ggplot2_boxplot_dataframe$X1, levels = final_ggplot2_boxplot_dataframe$X1[order(final_ggplot2_boxplot_dataframe$X1)])
    return(final_ggplot2_boxplot_dataframe)
  })
  
  #Recalculation of fin_0 to only selected samples and then using it to create dataframe that can be used for density plotting
  densityplot_datatable<- reactive({
    if (values$marker_correction == FALSE) {
      selectedsamples_dataframe<- as.data.frame(maxima_datatable_final_0()[,as.numeric(coordinates_from_memory()$sample_nu[!is.na(coordinates_from_memory()$sample_nu)])])
      slopesamples_dataframe<- data.frame(NA,NA)
      colnames(slopesamples_dataframe)<- c("bp","sample")
      for (i in c(1:ncol(selectedsamples_dataframe))) {
        selectedsamples_dataframe[[i]]
        df_for_repeat_i_sample<- as.data.frame(cbind(maxima_dataframe_corrected()$bp, selectedsamples_dataframe[[i]]))
        df_for_repeat_i_sample<- df_for_repeat_i_sample[sample_datatable()[i,2]:sample_datatable()[i,3],]#select just the xmin,xmax area with the use of sample_datatable information
        if (max(df_for_repeat_i_sample[[2]]) < 100) {
          df_for_repeat_i_sample[[2]]<- df_for_repeat_i_sample[[2]]*20
        } else if (max(df_for_repeat_i_sample[[2]]) > 100 & max(df_for_repeat_i_sample[[2]]) < 1000) {
          df_for_repeat_i_sample[[2]]<- df_for_repeat_i_sample[[2]]
        } else {
          df_for_repeat_i_sample[[2]]<- df_for_repeat_i_sample[[2]]/(10^(floor(log10(max(df_for_repeat_i_sample[[2]])))-2))
        }
        df_for_repeat_i_sample[df_for_repeat_i_sample < 1]<- 1
        repeated_nu<- rep(df_for_repeat_i_sample[[1]],df_for_repeat_i_sample[[2]])
        slope_with_sample<- data.frame(repeated_nu,colnames(selectedsamples_dataframe)[i])
        colnames(slope_with_sample)<- c("bp","sample")
        slopesamples_dataframe<- rbind(slopesamples_dataframe, slope_with_sample)
      }
      slopesamples_dataframe<- slopesamples_dataframe[-1,]
      if (length(colnames(final_boxplot_datatable())) == 1){
        slopesamples_dataframe[2]<- c("")
      }
      return(slopesamples_dataframe)
    } else {
      selectedsamples_dataframe<- as.data.frame(maxima_datatable_final_0()[,as.numeric(coordinates_from_memory()$sample_nu[!is.na(coordinates_from_memory()$sample_nu)])])
      slopesamples_dataframe<- data.frame(NA,NA)
      colnames(slopesamples_dataframe)<- c("bp","sample")
      dataframe_corresponding_bp<- dataframe_bp_recalculated_after_mark_correction()[,-1]
      dataframe_corresponding_bp<- dataframe_corresponding_bp[,as.numeric(coordinates_from_memory()$sample_nu[!is.na(coordinates_from_memory()$sample_nu)])]
      for (i in c(1:ncol(selectedsamples_dataframe))) {
        selectedsamples_dataframe[[i]]
        if (is.data.frame(dataframe_corresponding_bp) == TRUE) {
          vector_corresponding_bp<- dataframe_corresponding_bp[[i]]
        } else {
          vector_corresponding_bp<- as.numeric(dataframe_corresponding_bp)
        }
        bp_as_vector<- vector_corresponding_bp
        df_for_repeat_i_sample<- as.data.frame(cbind(bp_as_vector, selectedsamples_dataframe[[i]]))
        df_for_repeat_i_sample<- df_for_repeat_i_sample[sample_datatable()[i,2]:sample_datatable()[i,3],]#select just the xmin,xmax area with the use of sample_datatable information
        if (max(df_for_repeat_i_sample[[2]]) < 100) {
          df_for_repeat_i_sample[[2]]<- df_for_repeat_i_sample[[2]]*20
        } else if (max(df_for_repeat_i_sample[[2]]) > 100 & max(df_for_repeat_i_sample[[2]]) < 1000) {
          df_for_repeat_i_sample[[2]]<- df_for_repeat_i_sample[[2]]
        } else {
          df_for_repeat_i_sample[[2]]<- df_for_repeat_i_sample[[2]]/(10^(floor(log10(max(df_for_repeat_i_sample[[2]])))-2))
        }
        df_for_repeat_i_sample[df_for_repeat_i_sample < 1]<- 1
        repeated_nu<- rep(df_for_repeat_i_sample[[1]],df_for_repeat_i_sample[[2]])
        slope_with_sample<- data.frame(repeated_nu,colnames(selectedsamples_dataframe)[i])
        colnames(slope_with_sample)<- c("bp","sample")
        slopesamples_dataframe<- rbind(slopesamples_dataframe, slope_with_sample)
      }
      slopesamples_dataframe<- slopesamples_dataframe[-1,]
      if (length(colnames(final_boxplot_datatable())) == 1){
        slopesamples_dataframe[2]<- c("")
      }
      return(slopesamples_dataframe)
    }
  })
  
  #Final_boxplot_datatable with only medians for ggplot2
  final_ggplot2_boxplot_datatable_median<- reactive ({
    final_ggplot2_boxplot_dataframe<- as.data.frame(melt(final_boxplot_datatable()[3,]))
    if (length(colnames(final_ggplot2_boxplot_dataframe)) == 1) {
      final_ggplot2_boxplot_dataframe<- cbind(final_ggplot2_boxplot_dataframe,c(""))
      colnames(final_ggplot2_boxplot_dataframe)<- c("value","variable")
    }
    return(final_ggplot2_boxplot_dataframe)
  })
  
  #Final_boxplot_datatable with only 1st quartile for ggplot2
  final_ggplot2_boxplot_datatable_stquartile<- reactive ({
    final_ggplot2_boxplot_dataframe<- as.data.frame(melt(final_boxplot_datatable()[2,]))
    if (length(colnames(final_ggplot2_boxplot_dataframe)) == 1) {
      final_ggplot2_boxplot_dataframe<- cbind(final_ggplot2_boxplot_dataframe,c(""))
      colnames(final_ggplot2_boxplot_dataframe)<- c("value","variable")
    }
    return(final_ggplot2_boxplot_dataframe)
  })
  #Final_boxplot_datatable with only 3rd quartile for ggplot2
  final_ggplot2_boxplot_datatable_rdquartile<- reactive ({
    final_ggplot2_boxplot_dataframe<- as.data.frame(melt(final_boxplot_datatable()[4,]))
    if (length(colnames(final_ggplot2_boxplot_dataframe)) == 1) {
      final_ggplot2_boxplot_dataframe<- cbind(final_ggplot2_boxplot_dataframe,c(""))
      colnames(final_ggplot2_boxplot_dataframe)<- c("value","variable")
    }
    return(final_ggplot2_boxplot_dataframe)
  })
  
  #Making own breaks for the plot by the user -> at first auto
  max_y_graph_break<- reactive({
    max_y_break<- round_any(max(unlist(final_boxplot_datatable())), 1000, f = ceiling)
    return(max_y_break)
  })
  
  #Input max_y_scale of the result plot
  output$max_input<- renderUI({numericInput("max", label = "Max (bp)", max_y_graph_break(), 
                                              width = 77)})
  
  #Making own breaks for the plot by the user -> at first auto
  min_y_graph_break<- reactive({
    min_y_break<- round_any(min(unlist(final_boxplot_datatable())), 1000, f = floor)
    return(min_y_break)
  })
  
  #Input max_y_scale of the result plot
  output$min_input<- renderUI({numericInput("min", label = "Min (bp)", min_y_graph_break(), 
                                            width = 77)})
  
  
  output$result_plot<- renderPlot({
    all_boxplot_function <- function (x){
      ggplot2::ggplot(final_ggplot2_boxplot_datatable_transposed()) +
        geom_boxplot(aes(x = fct_inorder(X1), ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`), stat = "identity", colour = "black", fill = "white") +
        ylab("length [bp]") +
        scale_y_continuous(breaks=seq(0, input$max, x), limits = c(input$min, input$max), expand = c(0, 0)) +
        xlab("") +
        theme_bw() +
        theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
    }
    all_violinplot_function <- function (x){
      ggplot2::ggplot(densityplot_datatable(), aes(x = fct_inorder(sample), y = bp)) +
        geom_violin(colour = "black", fill = "white") +
        geom_point(data = final_ggplot2_boxplot_datatable_median(), aes(x = fct_inorder(variable), y = value), shape = 15, colour = "black", cex=2) +
        geom_point(data = final_ggplot2_boxplot_datatable_stquartile(), aes(x = fct_inorder(variable), y = value), shape = 6, colour = "black") +
        geom_point(data = final_ggplot2_boxplot_datatable_rdquartile(), aes(x = fct_inorder(variable), y = value), shape = 2, colour = "black") +
        ylab("length [bp]") +
        scale_y_continuous(breaks=seq(0, input$max, x), limits = c(input$min, input$max), expand = c(0, 0)) +
        xlab("") +
        theme_bw() +
        theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
    }
    group_boxplot_function <- function (x){
      ggplot2::ggplot(final_ggplot2_boxplot_datatable_transposed()) +
        geom_boxplot(aes(x = fct_inorder(X1), ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`), stat = "identity", colour = "black", fill = "white") +
        ylab("length [bp]") +
        scale_y_continuous(breaks=seq(0, input$max, x), limits = c(input$min, input$max), expand = c(0, 0)) +
        xlab("") +
        theme_bw() +
        facet_wrap(~ datatable_groups_renamed()$group_name, nrow = 1, scales = "free_x") +
        theme(axis.text=element_text(family="Roboto"), strip.text.x = element_text(family="Roboto", size = 11), axis.title=element_text(family="Roboto", size=14), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
    }
    group_violinplot_function <- function (x){
      ggplot2::ggplot(densityplot_groups_datatable(), aes(x = fct_inorder(sample), y = bp)) +
        geom_violin(colour = "black", fill = "white") +
        geom_point(aes(x = fct_inorder(sample), y = median), shape = 15, colour = "black", cex=2) +
        geom_point(aes(x = fct_inorder(sample), y = stquartile), shape = 6, colour = "black") +
        geom_point(aes(x = fct_inorder(sample), y = rdquartile), shape = 2, colour = "black") +
        facet_wrap(~ datatable_groups_density_renamed()$group_name, nrow = 1, scales = "free_x") +
        ylab("length [bp]") +
        scale_y_continuous(breaks=seq(0, input$max, x), limits = c(input$min, input$max), expand = c(0, 0)) +
        xlab("") +
        theme_bw() +
        theme(axis.text=element_text(family="Roboto"), strip.text.x = element_text(family="Roboto", size = 11), axis.title=element_text(family="Roboto", size=14), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
    }
    if (statisticalanalysis_on() == FALSE) {
      if (input$box_or_densityplots == "Boxplot") {
        if (input$yaxisrange == "Auto") {
          ggplot2::ggplot(final_ggplot2_boxplot_datatable_transposed()) +
            geom_boxplot(aes(x = fct_inorder(X1), ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`), stat = "identity", colour = "black", fill = "white") +
            ylab("length [bp]") +
            scale_y_continuous(limits = c(input$min, input$max), expand = c(0, 0)) +
            xlab("") +
            theme_bw() +
            theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
        } else if (input$yaxisrange == "250") {
          all_boxplot_function(250)
        } else if (input$yaxisrange == "500") {
          all_boxplot_function(500)
        } else if (input$yaxisrange == "750") {
          all_boxplot_function(750)
        } else if (input$yaxisrange == "1000") {
          all_boxplot_function(1000)
        } else if (input$yaxisrange == "2.5K") {
          all_boxplot_function(2500)
        }
      } else {
        if (input$yaxisrange == "Auto") {
          ggplot2::ggplot(densityplot_datatable(), aes(x = fct_inorder(sample), y = bp)) +
            geom_violin(colour = "black", fill = "white") +
            geom_point(data = final_ggplot2_boxplot_datatable_median(), aes(x = fct_inorder(variable), y = value), shape = 15, colour = "black", cex=2) +
            geom_point(data = final_ggplot2_boxplot_datatable_stquartile(), aes(x = fct_inorder(variable), y = value), shape = 6, colour = "black") +
            geom_point(data = final_ggplot2_boxplot_datatable_rdquartile(), aes(x = fct_inorder(variable), y = value), shape = 2, colour = "black") +
            ylab("length [bp]") +
            scale_y_continuous(limits = c(input$min, input$max), expand = c(0, 0)) +
            xlab("") +
            theme_bw() +
            theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
        } else if (input$yaxisrange == "250") {
          all_violinplot_function(250)
        } else if (input$yaxisrange == "500") {
          all_violinplot_function(500)
        } else if (input$yaxisrange == "750") {
          all_violinplot_function(750)
        } else if (input$yaxisrange == "1000") {
          all_violinplot_function(1000)
        } else if (input$yaxisrange == "2.5K") {
          all_violinplot_function(2500)
        }
      }
    } else {
      if (input$box_or_densityplots == "Boxplot") {
        if (input$yaxisrange == "Auto") {
          ggplot2::ggplot(final_ggplot2_boxplot_datatable_transposed()) +
            geom_boxplot(aes(x = fct_inorder(X1), ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`), stat = "identity", colour = "black", fill = "white") +
            ylab("length [bp]") +
            scale_y_continuous(limits = c(input$min, input$max), expand = c(0, 0)) +
            xlab("") +
            theme_bw() +
            facet_wrap(~ datatable_groups_renamed()$group_name, nrow = 1, scales = "free_x") +
            theme(axis.text=element_text(family="Roboto"), strip.text.x = element_text(family="Roboto", size = 11), axis.title=element_text(family="Roboto", size=14), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
        } else if (input$yaxisrange == "250") {
          group_boxplot_function(250)
        } else if (input$yaxisrange == "500") {
          group_boxplot_function(500)
        } else if (input$yaxisrange == "750") {
          group_boxplot_function(750)
        } else if (input$yaxisrange == "1000") {
          group_boxplot_function(1000)
        } else if (input$yaxisrange == "2.5K") {
          group_boxplot_function(2500)
        }
      } else {
        if (input$yaxisrange == "Auto") {
          ggplot2::ggplot(densityplot_groups_datatable(), aes(x = fct_inorder(sample), y = bp)) +
            geom_violin(colour = "black", fill = "white") +
            geom_point(aes(x = fct_inorder(sample), y = median), shape = 15, colour = "black", cex=2) +
            geom_point(aes(x = fct_inorder(sample), y = stquartile), shape = 6, colour = "black") +
            geom_point(aes(x = fct_inorder(sample), y = rdquartile), shape = 2, colour = "black") +
            facet_wrap(~ datatable_groups_density_renamed()$group_name, nrow = 1, scales = "free_x") +
            ylab("length [bp]") +
            scale_y_continuous(limits = c(input$min, input$max), expand = c(0, 0)) +
            xlab("") +
            theme_bw() +
            theme(axis.text=element_text(family="Roboto"), strip.text.x = element_text(family="Roboto", size = 11), axis.title=element_text(family="Roboto", size=14), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
        } else if (input$yaxisrange == "250") {
          group_violinplot_function(250)
        } else if (input$yaxisrange == "500") {
          group_violinplot_function(500)
        } else if (input$yaxisrange == "750") {
          group_violinplot_function(750)
        } else if (input$yaxisrange == "1000") {
          group_violinplot_function(1000)
        } else if (input$yaxisrange == "2.5K") {
          group_violinplot_function(2500)
        }
      }
    }
  })
  
  #300DPI plot to render in the new tab
  output$result_plot_300DPI<- renderPlot({
    all_boxplot_function <- function (x){
      ggplot2::ggplot(final_ggplot2_boxplot_datatable_transposed()) +
        geom_boxplot(aes(x = fct_inorder(X1), ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`), stat = "identity", colour = "black", fill = "white") +
        ylab("length [bp]") +
        scale_y_continuous(breaks=seq(0, input$max, x), limits = c(input$min, input$max), expand = c(0, 0)) +
        xlab("") +
        theme_bw() +
        theme(axis.text=element_text(family="Roboto", size=34.375), axis.title=element_text(family="Roboto", size=43.75), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
    }
    all_violinplot_function <- function (x){
      ggplot2::ggplot(densityplot_datatable(), aes(x = fct_inorder(sample), y = bp)) +
        geom_violin(colour = "black", fill = "white") +
        geom_point(data = final_ggplot2_boxplot_datatable_median(), aes(x = fct_inorder(variable), y = value), shape = 15, colour = "black", cex=2) +
        geom_point(data = final_ggplot2_boxplot_datatable_stquartile(), aes(x = fct_inorder(variable), y = value), shape = 6, colour = "black") +
        geom_point(data = final_ggplot2_boxplot_datatable_rdquartile(), aes(x = fct_inorder(variable), y = value), shape = 2, colour = "black") +
        ylab("length [bp]") +
        scale_y_continuous(breaks=seq(0, input$max, x), limits = c(input$min, input$max), expand = c(0, 0)) +
        xlab("") +
        theme_bw() +
        theme(axis.text=element_text(family="Roboto", size=34.375), axis.title=element_text(family="Roboto", size=43.75), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
    }
    group_boxplot_function <- function (x){
      ggplot2::ggplot(final_ggplot2_boxplot_datatable_transposed()) +
        geom_boxplot(aes(x = fct_inorder(X1), ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`), stat = "identity", colour = "black", fill = "white") +
        ylab("length [bp]") +
        scale_y_continuous(breaks=seq(0, input$max, x), limits = c(input$min, input$max), expand = c(0, 0)) +
        xlab("") +
        theme_bw() +
        facet_wrap(~ datatable_groups_renamed()$group_name, nrow = 1, scales = "free_x") +
        theme(axis.text=element_text(family="Roboto", size=34.375), strip.text.x = element_text(family="Roboto", size=34.375), axis.title=element_text(family="Roboto", size=43.75), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
    }
    group_violinplot_function <- function (x){
      ggplot2::ggplot(densityplot_groups_datatable(), aes(x = fct_inorder(sample), y = bp)) +
        geom_violin(colour = "black", fill = "white") +
        geom_point(aes(x = fct_inorder(sample), y = median), shape = 15, colour = "black", cex=2) +
        geom_point(aes(x = fct_inorder(sample), y = stquartile), shape = 6, colour = "black") +
        geom_point(aes(x = fct_inorder(sample), y = rdquartile), shape = 2, colour = "black") +
        facet_wrap(~ datatable_groups_density_renamed()$group_name, nrow = 1, scales = "free_x") +
        ylab("length [bp]") +
        scale_y_continuous(breaks=seq(0, input$max, x), limits = c(input$min, input$max), expand = c(0, 0)) +
        xlab("") +
        theme_bw() +
        theme(axis.text=element_text(family="Roboto", size=34.375), strip.text.x = element_text(family="Roboto", size=34.375), axis.title=element_text(family="Roboto", size=43.75), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
    }
    if (statisticalanalysis_on() == FALSE) {
      if (input$box_or_densityplots == "Boxplot") {
        if (input$yaxisrange == "Auto") {
          ggplot2::ggplot(final_ggplot2_boxplot_datatable_transposed()) +
            geom_boxplot(aes(x = fct_inorder(X1), ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`), stat = "identity", colour = "black", fill = "white") +
            ylab("length [bp]") +
            scale_y_continuous(limits = c(input$min, input$max), expand = c(0, 0)) +
            xlab("") +
            theme_bw() +
            theme(axis.text=element_text(family="Roboto", size=34.375), axis.title=element_text(family="Roboto", size=43.75), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
        } else if (input$yaxisrange == "250") {
          all_boxplot_function(250)
        } else if (input$yaxisrange == "500") {
          all_boxplot_function(500)
        } else if (input$yaxisrange == "750") {
          all_boxplot_function(750)
        } else if (input$yaxisrange == "1000") {
          all_boxplot_function(1000)
        } else if (input$yaxisrange == "2.5K") {
          all_boxplot_function(2500)
        }
      } else {
        if (input$yaxisrange == "Auto") {
          ggplot2::ggplot(densityplot_datatable(), aes(x = fct_inorder(sample), y = bp)) +
            geom_violin(colour = "black", fill = "white") +
            geom_point(data = final_ggplot2_boxplot_datatable_median(), aes(x = fct_inorder(variable), y = value), shape = 15, colour = "black", cex=2) +
            geom_point(data = final_ggplot2_boxplot_datatable_stquartile(), aes(x = fct_inorder(variable), y = value), shape = 6, colour = "black") +
            geom_point(data = final_ggplot2_boxplot_datatable_rdquartile(), aes(x = fct_inorder(variable), y = value), shape = 2, colour = "black") +
            ylab("length [bp]") +
            scale_y_continuous(limits = c(input$min, input$max), expand = c(0, 0)) +
            xlab("") +
            theme_bw() +
            theme(axis.text=element_text(family="Roboto", size=34.375), axis.title=element_text(family="Roboto", size=43.75), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
        } else if (input$yaxisrange == "250") {
          all_violinplot_function(250)
        } else if (input$yaxisrange == "500") {
          all_violinplot_function(500)
        } else if (input$yaxisrange == "750") {
          all_violinplot_function(750)
        } else if (input$yaxisrange == "1000") {
          all_violinplot_function(1000)
        } else if (input$yaxisrange == "2.5K") {
          all_violinplot_function(2500)
        }
      }
    } else {
      if (input$box_or_densityplots == "Boxplot") {
        if (input$yaxisrange == "Auto") {
          ggplot2::ggplot(final_ggplot2_boxplot_datatable_transposed()) +
            geom_boxplot(aes(x = fct_inorder(X1), ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`), stat = "identity", colour = "black", fill = "white") +
            ylab("length [bp]") +
            scale_y_continuous(limits = c(input$min, input$max), expand = c(0, 0)) +
            xlab("") +
            theme_bw() +
            facet_wrap(~ datatable_groups_renamed()$group_name, nrow = 1, scales = "free_x") +
            theme(axis.text=element_text(family="Roboto", size=34.375), strip.text.x = element_text(family="Roboto", size=34.375), axis.title=element_text(family="Roboto", size=43.75), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
        } else if (input$yaxisrange == "250") {
          group_boxplot_function(250)
        } else if (input$yaxisrange == "500") {
          group_boxplot_function(500)
        } else if (input$yaxisrange == "750") {
          group_boxplot_function(750)
        } else if (input$yaxisrange == "1000") {
          group_boxplot_function(1000)
        } else if (input$yaxisrange == "2.5K") {
          group_boxplot_function(2500)
        }
      } else {
        if (input$yaxisrange == "Auto") {
          ggplot2::ggplot(densityplot_groups_datatable(), aes(x = fct_inorder(sample), y = bp)) +
            geom_violin(colour = "black", fill = "white") +
            geom_point(aes(x = fct_inorder(sample), y = median), shape = 15, colour = "black", cex=2) +
            geom_point(aes(x = fct_inorder(sample), y = stquartile), shape = 6, colour = "black") +
            geom_point(aes(x = fct_inorder(sample), y = rdquartile), shape = 2, colour = "black") +
            facet_wrap(~ datatable_groups_density_renamed()$group_name, nrow = 1, scales = "free_x") +
            ylab("length [bp]") +
            scale_y_continuous(limits = c(input$min, input$max), expand = c(0, 0)) +
            xlab("") +
            theme_bw() +
            theme(axis.text=element_text(family="Roboto", size=34.375), strip.text.x = element_text(family="Roboto", size=34.375), axis.title=element_text(family="Roboto", size=43.75), axis.ticks = element_line(colour = "black"), axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major.x= element_blank())
        } else if (input$yaxisrange == "250") {
          group_violinplot_function(250)
        } else if (input$yaxisrange == "500") {
          group_violinplot_function(500)
        } else if (input$yaxisrange == "750") {
          group_violinplot_function(750)
        } else if (input$yaxisrange == "1000") {
          group_violinplot_function(1000)
        } else if (input$yaxisrange == "2.5K") {
          group_violinplot_function(2500)
        }
      }
    }
  }, res = 300)
  
  #Input width of the result plot
  output$width_input<- renderUI({numericInput("width", label = "Width (px)", 
                                              value = if(length(colnames(final_boxplot_datatable())) == 1){
                                                92
                                              } else {
                                                60+40*(length(colnames(final_boxplot_datatable())))
                                              }, 
                                              width = 77)})
  
  #Output for the IntensityAnalyser
  output$result<- renderUI({plotOutput("result_plot", width = input$width, height = input$height)})
  
  #Output for the new plot of hi-res
  output$result_300DPI<- renderUI({plotOutput("result_plot_300DPI", width = input$width*3.525, height = input$height*3.125)})
  
  #Text in Graphic result based on selected plot type
  output$explanation_graphicresult_ui<- renderUI({
    if (input$box_or_densityplots == "Boxplot") {
      img(src="help_graphic_result_boxplot.png")
    } else {
      img(src="help_graphic_result_violinplot.png")
    }
  })
  
  #Recalculation of median, Q1 and Q3(taken as (medain-Q1+median)) to mean SD based on Wan et al. 2014 doi: 10.1186/1471-2288-14-135 for each line, n samples taken as a range for each sample
  mean_sd_n_datatable<- reactive({
    if (is.data.frame(cumsum_maxima_datatable())) {
      mean_sd_n_dataframe<- final_boxplot_datatable()[1:3,]
      for (i in 1:ncol(final_boxplot_datatable())) {
        mean<- (final_boxplot_datatable()[2,i]+final_boxplot_datatable()[3,i]+final_boxplot_datatable()[4,i])/3
        mean_sd_n_dataframe[1,i]<- mean
      }
      for (i in 1:ncol(final_boxplot_datatable())) {
        n<- (sample_datatable()[i,3]-sample_datatable()[i,2])
        mean_sd_n_dataframe[3,i]<- n
      }
      for (i in 1:ncol(final_boxplot_datatable())) {
        sd<- (final_boxplot_datatable()[4,i]-final_boxplot_datatable()[2,i])/(2*qnorm((0.75*mean_sd_n_dataframe[3,i]-0.125)/(mean_sd_n_dataframe[3,i]+0.25)))
        mean_sd_n_dataframe[2,i]<- sd
      }
      rownames(mean_sd_n_dataframe)<- c("mean","SD","n")
    } else {
      n<- (sample_datatable()[1,3]-sample_datatable()[1,2])
      mean_sd_n_dataframe<- data.frame(c((final_boxplot_datatable()[2,1]+final_boxplot_datatable()[3,1]+final_boxplot_datatable()[4,1])/3,
                                         (final_boxplot_datatable()[4,1]-final_boxplot_datatable()[2,1])/(2*qnorm((0.75*n-0.125)/(n+0.25))),
                                       n
                                       ))
      rownames(mean_sd_n_dataframe)<- c("mean","SD","n")
      colnames(mean_sd_n_dataframe)<- ""
    }
    return(mean_sd_n_dataframe)
  })
  
  #Button for getting information about statistical analysis
  observeEvent(input$info_statistics, {
    wan_et_al <- 
      shinyalert(html = TRUE,paste(
        h6(tags$div("For statistical analysis you need to:", tags$br(),
        "- have at least 3 samples in each group (however, it is STRONGLY recommended to have at least 5 samples)", tags$br(),
        "- fill the group name for all samples", tags$br(),
        "- have at least 2 groups", tags$br(),
        "- fill the Control group name box with a name used in the table"), align = "justify"),
        h6("Groups can be compared by 4 different types of tests (2 groups are always compared by Welch's t-test). There is a possibility to compare data by multiple Welch's t-test either vs control group (CTR) or 
           between each group. However, those tests do NOT adjust p-value for multiple testing. More stringent testing is provided by Tamhane-Dunnett test
           that compares groups to the control group and by Games-Howell test that compares groups between each other. Due to the possibility of type II error, it is recommended to use multiple 
           Welch's t-test. For the same reason, significance is set as: * < 0.1, ** < 0.05, *** < 0.01 p-value.", align = "justify"),
        tags$hr(),
        h6("Before the statistical evaluation, samples are recalculated according to", a("Wan et al. 2014", href="https://bmcmedresmethodol.biomedcentral.com/track/pdf/10.1186/1471-2288-14-135/", target="_blank"),
           "using information about their interquartile range and median (n = number of pixels in the whole range where the boxplot is calculated).
           Calculated mean, SD and n of each sample is then used to combine samples according to their group name according to", a("Cochrane Handbook", href="https://handbook-5-1.cochrane.org/chapter_7/table_7_7_a_formulae_for_combining_groups.htm", target="_blank"),
           "Resulting mean and SD (n = number of samples in a group) for each group are used to create mock data with these parameters for the statistical calculation mentioned above.", align = "justify")
      ), type = "info", confirmButtonCol =  "#337ab7")
    })
    
  #Button for Evaluation of rhandsontable for statistical analysis
  output$evaluate_button <- renderUI({
    if (values$statistical_value == FALSE) {
      actionButton("evaluate_table", "Evaluate", width = 115)
    } else {
      actionButton("nothing", "Evaluate", width = 115, style="background-color: #E5E5E5")
    }
    
  })
  
  #Evaluate current table button with analysing, if you can even do it
  values$statistical_value<- FALSE
  values$statistical_counter<- 0
  values$grouptable_datatable_after_click<- 0
  observeEvent(input$evaluate_table, {
    if(nrow(datatable_groups()) < 6) {
      shinyalert("Amount of samples too low", "Table doesn't include enough samples to carry out statistical evaluation. There must be at least three samples with the group name same as in the box Control group name and three samples belonging to other group.", type = "error", confirmButtonCol =  "#337ab7")
    } else if (!is.element(input$control_group_name,datatable_groups()$group_name)) {
      shinyalert("No control group included", "Table doesn't include control group (samples with the same name as in a box Control group name) to which other groups are compared to. There must be at least three samples with the group name same as in the box Control group name and three samples belonging to other group.", type = "error", confirmButtonCol =  "#337ab7")
    } else if (any(datatable_groups()$group_name == "")) {
      shinyalert("Sample(s) not belonging to any group", "Table has samples that are not part of any group. There must be at least three samples with the group name same as in the box Control group name and three samples belonging to other group.", type = "error", confirmButtonCol =  "#337ab7")
    } else if (all(as.data.frame(table(datatable_groups()$group_name))[[2]] >= 3, na.rm = TRUE) == FALSE) {
      shinyalert("Groups with just one/two sample(s)", "Table contains group names with just one/two sample(s). There must be at least three samples with the group name same as in the box Control group name and three samples belonging to other group.", type = "error", confirmButtonCol =  "#337ab7")
    } else if (length(unique(datatable_groups()$group_name)) < 2) {
      shinyalert("Only one group", "Table contains just one group and therefore there is nothing to compare. There must be at least three samples with the group name same as in the box Control group name and three samples belonging to other group.", type = "error", confirmButtonCol =  "#337ab7")
    } else {
      shinyjs::disable("control_group_name")
      values$statistical_value<- TRUE
      satistical_counting<- values$statistical_counter
      values$statistical_counter<- satistical_counting + 1
      grouptable_datatable_click<- as.data.frame(hot_to_r(input$grouptable))
      values$grouptable_datatable_after_click<- grouptable_datatable_click
    }
  })
  
  #Evaluate current table button with analysing, if you can even do it
  observeEvent(input$restart_statistical_evaluation, {
    shinyjs::enable("control_group_name")
    values$statistical_value<- FALSE
    satistical_counting<- values$statistical_counter
    values$statistical_counter<- satistical_counting + 1
    grouptable_datatable_click<- as.data.frame(hot_to_r(input$grouptable))
    values$grouptable_datatable_after_click<- grouptable_datatable_click
  })
  
  #If to allow statistical_evaluation
  statisticalanalysis_on<- reactive({
    if (values$statistical_value == TRUE) {
      statisticalanalysis_value<- TRUE
    } else if (values$statistical_value == FALSE) {
      statisticalanalysis_value<- FALSE
    }
    return(statisticalanalysis_value)
  })
  
  #Dataframe to fill with sample number and group name
  grouptable<- reactiveValues()
  grouptable$column_numbers <- reactiveVal({0})
  grouptable$NAs <- reactiveVal({""})
  group_dataframe<- reactive({
    setNames(data.frame(grouptable$column_numbers(), grouptable$NAs(), stringsAsFactors = FALSE), c("sample no.", "group name"))
    })
  
  #Editable dataframe by rhandsontable for group names
  output$grouptable<- renderRHandsontable({
    if (values$statistical_counter < 1) {
      if (length(grouptable$column_numbers()) >= 12) {
        rhandsontable(group_dataframe(), readOnly = statisticalanalysis_on(), rowHeaderWidth = 0, height = 282, contextMenu = FALSE)%>%
          hot_cols(colWidths = 90)%>%
          hot_cols(format = "0")%>%
          hot_col("sample no.", halign = "htCenter", readOnly = TRUE)
      } else {
        rhandsontable(group_dataframe(), readOnly = statisticalanalysis_on(), rowHeaderWidth = 0, height = 282, contextMenu = FALSE)%>%
          hot_cols(colWidths = 100)%>%
          hot_cols(format = "0")%>%
          hot_col("sample no.", halign = "htCenter", readOnly = TRUE)
      }
    } else {
      if (length(grouptable$column_numbers()) >= 12) {
        rhandsontable(values$grouptable_datatable_after_click, readOnly = statisticalanalysis_on(), rowHeaderWidth = 0, height = 282, contextMenu = FALSE)%>%
          hot_cols(colWidths = 90)%>%
          hot_cols(format = "0")%>%
          hot_col("sample no.", halign = "htCenter", readOnly = TRUE)
      } else {
        rhandsontable(values$grouptable_datatable_after_click, readOnly = statisticalanalysis_on(), rowHeaderWidth = 0, height = 282, contextMenu = FALSE)%>%
          hot_cols(colWidths = 100)%>%
          hot_cols(format = "0")%>%
          hot_col("sample no.", halign = "htCenter", readOnly = TRUE)
      }
    }
  })
  
  #Taking rhandsontable data and creating dataframe for statistical evaluation (but first removing possibility of error message due to different length of arguments)
  datatable_groups<- reactive({
    input_grouptable<- as.data.frame(hot_to_r(input$grouptable))
    input_grouptable<- input_grouptable[1:length(grouptable$column_numbers()),]
    colnames(input_grouptable)<- c("sample_number", "group_name")
    return(input_grouptable)
  })
  
  #Combining mean and SD and n in a group to one mean and SD (n still used those from mean_sd_n_datatable) for report
  datatable_mean_groups<- reactive({
    if(nrow(datatable_groups()) < 6) {
      dataframe_significance <- "Can't be calculated."
      return(dataframe_significance)
    } else if (!is.element(input$control_group_name,datatable_groups()$group_name)) {
      dataframe_significance <- "Can't be calculated."
      return(dataframe_significance)
    } else if (all(as.data.frame(table(datatable_groups()$group_name))[[2]] >= 3, na.rm = TRUE) == FALSE) {
      dataframe_significance <- "Can't be calculated."
      return(dataframe_significance)
    } else if (length(unique(datatable_groups()$group_name)) < 2) {
      dataframe_significance <- "Can't be calculated."
      return(dataframe_significance)
    } else if (any(datatable_groups()$group_name == "")) {
      dataframe_significance <- "Can't be calculated."
    } else {
      #Grouping of samples mean SD to one mean SD based on Cochrane's handbook (7.7.a)
      #Functions for calculating mean and SD and generation of random data for ANOVA and post-hoc analysis
      mean_calc<- function (mean_one, mean_two, n_one, n_two){
        (n_one*mean_one+n_two*mean_two)/(n_one+n_two)
      }
      sd_calc<- function (mean_one, mean_two, n_one, n_two, sd_one, sd_two){
        ((((n_one-1)*(sd_one^2))+((n_two-1)*(sd_two^2))+((n_one*n_two)/(n_one+n_two))*((mean_one^2)+(mean_two^2)-(2*mean_one*mean_two)))/(n_one+n_two-1))^0.5
      } 
      gen_data <- function(means, sds, samplesizes, unique_group_names){
        group_name <- factor(rep(unique_group_names, samplesizes), levels=unique_group_names)
        dat <- lapply(1:length(unique_group_names), function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
        values <- do.call(rbind, dat)
        out <- data.frame(group_name = group_name, values = values)
        out
      }
      
      group_name_vector<- datatable_groups()[[2]]
      unique_groups_all<- unique(group_name_vector)
      mean_sd_n_group_dataframe<- as.data.frame(matrix(nrow = 3, ncol = length(unique_groups_all)))
      colnames(mean_sd_n_group_dataframe)<- unique_groups_all
      rownames(mean_sd_n_group_dataframe)<- c("mean","SD","n")
      for (a in 1:length(unique_groups_all)) {
        sample_nu_group<- datatable_groups()$sample_number[datatable_groups()$group_name == unique_groups_all[a]]
        #Save values of mean, sd and n of the first from sample_nu_group so they will be recalculated as 1+2, (1+2)+3...  
        m<- mean_sd_n_datatable()[1, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[1])]
        s<- mean_sd_n_datatable()[2, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[1])]
        n<- mean_sd_n_datatable()[3, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[1])]
        n_true<- length(group_name_vector[group_name_vector == unique_groups_all[a]])
        for (b in 2:length(group_name_vector[group_name_vector == unique_groups_all[a]])) {
          m_next<- mean_sd_n_datatable()[1, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[b])] 
          s_next<- mean_sd_n_datatable()[2, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[b])]
          n_next<- mean_sd_n_datatable()[3, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[b])]
          s<- sd_calc(m,m_next,n,n_next,s,s_next)
          m<- mean_calc(m,m_next,n,n_next)
          n<- n+n_next
          
        }
        mean_sd_n_group_dataframe[1,a]<- m
        mean_sd_n_group_dataframe[2,a]<- s
        mean_sd_n_group_dataframe[3,a]<- n_true
      }
      return(mean_sd_n_group_dataframe)
    }
  })
  
  #Combining mean and SD and n in a group to one mean and SD (n still used those from mean_sd_n_datatable) +
  #statistical analysis either by multiple Welch's t-test (vs control, vs all) or by ANOVA followed by Tamhane-Dunnett or Games-Howell 
  #(can be chosen) -> here n considered number of samples in each group
  datatable_significance<- reactive({
    if(nrow(datatable_groups()) < 6) {
      dataframe_significance <- "Can't be calculated."
      return(dataframe_significance)
    } else if (!is.element(input$control_group_name,datatable_groups()$group_name)) {
      dataframe_significance <- "Can't be calculated."
      return(dataframe_significance)
    } else if (all(as.data.frame(table(datatable_groups()$group_name))[[2]] >= 3, na.rm = TRUE) == FALSE) {
      dataframe_significance <- "Can't be calculated."
      return(dataframe_significance)
    } else if (length(unique(datatable_groups()$group_name)) < 2) {
      dataframe_significance <- "Can't be calculated."
      return(dataframe_significance)
    } else if (any(datatable_groups()$group_name == "")) {
      dataframe_significance <- "Can't be calculated."
    } else {
      #Grouping of samples mean SD to one mean SD based on Cochrane's handbook (7.7.a)
      #Functions for calculating mean and SD and generation of random data for ANOVA and post-hoc analysis
      mean_calc<- function (mean_one, mean_two, n_one, n_two){
        (n_one*mean_one+n_two*mean_two)/(n_one+n_two)
      }
      sd_calc<- function (mean_one, mean_two, n_one, n_two, sd_one, sd_two){
        ((((n_one-1)*(sd_one^2))+((n_two-1)*(sd_two^2))+((n_one*n_two)/(n_one+n_two))*((mean_one^2)+(mean_two^2)-(2*mean_one*mean_two)))/(n_one+n_two-1))^0.5
      } 
      gen_data <- function(means, sds, samplesizes, unique_group_names){
        group_name <- factor(rep(unique_group_names, samplesizes), levels=unique_group_names)
        dat <- lapply(1:length(unique_group_names), function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
        values <- do.call(rbind, dat)
        out <- data.frame(group_name = group_name, values = values)
        out
      }
      
      group_name_vector<- datatable_groups()[[2]]
      unique_groups_all<- unique(group_name_vector)
      mean_sd_n_group_dataframe<- as.data.frame(matrix(nrow = 3, ncol = length(unique_groups_all)))
      colnames(mean_sd_n_group_dataframe)<- unique_groups_all
      rownames(mean_sd_n_group_dataframe)<- c("mean","SD","n")
      for (a in 1:length(unique_groups_all)) {
        sample_nu_group<- datatable_groups()$sample_number[datatable_groups()$group_name == unique_groups_all[a]]
        #Save values of mean, sd and n of the first from sample_nu_group so they will be recalculated as 1+2, (1+2)+3...  
        m<- mean_sd_n_datatable()[1, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[1])]
        s<- mean_sd_n_datatable()[2, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[1])]
        n<- mean_sd_n_datatable()[3, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[1])]
        n_true<- length(group_name_vector[group_name_vector == unique_groups_all[a]])
        for (b in 2:length(group_name_vector[group_name_vector == unique_groups_all[a]])) {
          m_next<- mean_sd_n_datatable()[1, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[b])] 
          s_next<- mean_sd_n_datatable()[2, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[b])]
          n_next<- mean_sd_n_datatable()[3, colnames(mean_sd_n_datatable()) == as.character(sample_nu_group[b])]
          s<- sd_calc(m,m_next,n,n_next,s,s_next)
          m<- mean_calc(m,m_next,n,n_next)
          n<- n+n_next
          
        }
        mean_sd_n_group_dataframe[1,a]<- m
        mean_sd_n_group_dataframe[2,a]<- s
        mean_sd_n_group_dataframe[3,a]<- n_true
      }
      mean_sd_n_group_dataframe <- mean_sd_n_group_dataframe[,c(which(colnames(mean_sd_n_group_dataframe)==input$control_group_name),which(colnames(mean_sd_n_group_dataframe)!=input$control_group_name))]
      simulated_data <- gen_data(unname(unlist(mean_sd_n_group_dataframe[1,])), 
                                 unname(unlist(mean_sd_n_group_dataframe[2,])),
                                 unname(unlist(mean_sd_n_group_dataframe[3,])),
                                 colnames(mean_sd_n_group_dataframe))
      #Analysis of simulated_data by Welch's test (number of groups = 2), or ANOVA with posthoc Tukey or Dunnett based on what user has selected 
      if (length(colnames(mean_sd_n_group_dataframe)) < 3) {
        welch_result_df<- data.frame(colnames(mean_sd_n_group_dataframe), rep("",(length(colnames(mean_sd_n_group_dataframe)))), stringsAsFactors=FALSE)
        colnames(welch_result_df) <- c("group_name", "significance")
        welch_result<- t.test(simulated_data$values[simulated_data$group_name == colnames(mean_sd_n_group_dataframe)[1]],
                            simulated_data$values[simulated_data$group_name == colnames(mean_sd_n_group_dataframe)[2]]
                            )
        welch_result<- welch_result$p.value
        welch_result[as.numeric(welch_result) <= 0.01]<- "***"
        welch_result[as.numeric(welch_result) <= 0.05]<- "**"
        welch_result[as.numeric(welch_result) <= 0.1]<- "*"
        welch_result[as.numeric(welch_result) > 0.1]<- "ns"
        welch_result<- c("",welch_result)
        welch_result_df[,2]<- welch_result
        return(welch_result_df)
      } else {
        if (input$variant_of_statistics == "Tamhane-Dunnett (vs control)") {
          dunnett_result_df<- data.frame(colnames(mean_sd_n_group_dataframe), rep("",(length(colnames(mean_sd_n_group_dataframe)))), stringsAsFactors=FALSE)
          colnames(dunnett_result_df) <- c("group_name", "significance")
          set.seed(20140123)
          dunnett_result<- tamhaneDunnettTest(values ~ group_name, data = simulated_data)
          dunnett_result<- as.vector(dunnett_result[[3]])
          dunnett_result[as.numeric(dunnett_result) <= 0.01]<- "***"
          dunnett_result[as.numeric(dunnett_result) <= 0.05]<- "**"
          dunnett_result[as.numeric(dunnett_result) <= 0.1]<- "*"
          dunnett_result[as.numeric(dunnett_result) > 0.1]<- "ns"
          dunnett_result<- c("",dunnett_result)
          dunnett_result_df[,2]<- dunnett_result
          return(dunnett_result_df)
        } else if (input$variant_of_statistics == "Games-Howell (vs all)") {
          tukey <- function(	data,					
                             group,					
                             method=c("Tukey", "Games-Howell"))	
          {
            OK <- complete.cases(data, group)			
            data <- data[OK]
            group <- factor(group[OK])
            n <- tapply(data, group, length)			
            a <- length(n)						
            phi.e <- sum(n)-a					
            Mean <- tapply(data, group, mean)			
            Variance <- tapply(data, group, var)			
            result1 <- cbind(n, Mean, Variance)			
            rownames(result1) <- paste("Group", 1:a, sep="")
            method <- match.arg(method)
            if (method == "Tukey") {				
              v.e <- sum((n-1)*Variance)/phi.e		
              t <- combn(a, 2, function(ij)			
                abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
              p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)
              Tukey <- cbind(t, p)					
              rownames(Tukey) <- combn(a, 2, paste, collapse="-")
              return(list(result1=result1, Tukey=Tukey, phi=phi.e, v=v.e))
            }
            else {							
              t.df <- combn(a, 2, function(ij) {		
                t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
                df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
                return(c(t, df))} )
              t <- t.df[1,]
              df <- t.df[2,]
              p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)	
              Games.Howell <- cbind(t, df, p)			
              rownames(Games.Howell) <- combn(a, 2, paste, collapse="-")
              return(list(result1=result1, Games.Howell=Games.Howell))
            }
          }
          
          games_howell_result_df<- data.frame(colnames(mean_sd_n_group_dataframe), rep("",(length(colnames(mean_sd_n_group_dataframe)))), stringsAsFactors=FALSE)
          colnames(games_howell_result_df) <- c("group_name", "significance")
          games_howell_result<- tukey(data = simulated_data$values, group = simulated_data$group_name, method = "Games-Howell")
          games_howell_result<- games_howell_result$`Games.Howell`
          p_values_from_games_howell<- games_howell_result[,3]
          games_howell_letters<- multcompLetters(p_values_from_games_howell, threshold = 0.1)
          games_howell_letters<- games_howell_letters$`Letters`
          games_howell_result_df[,2]<- games_howell_letters
          return(games_howell_result_df) 
        } else if (input$variant_of_statistics == "multiple Welch's t-test (vs control)") {
          multiple_welch_c_df<- data.frame(colnames(mean_sd_n_group_dataframe), rep("",(length(colnames(mean_sd_n_group_dataframe)))), stringsAsFactors=FALSE)
          colnames(multiple_welch_c_df) <- c("group_name", "significance")
          multiple_welch_c_result<- NA
          for (i in 2:length(colnames(mean_sd_n_group_dataframe))) {
            t_test<- t.test(simulated_data[[2]][simulated_data[[1]] == colnames(mean_sd_n_group_dataframe)[1]],simulated_data[[2]][simulated_data[[1]] == colnames(mean_sd_n_group_dataframe)[i]])
            multiple_welch_c_result[i-1]<- t_test$p.value
          }
          multiple_welch_c_result[as.numeric(multiple_welch_c_result) <= 0.01]<- "***"
          multiple_welch_c_result[as.numeric(multiple_welch_c_result) <= 0.05]<- "**"
          multiple_welch_c_result[as.numeric(multiple_welch_c_result) <= 0.1]<- "*"
          multiple_welch_c_result[as.numeric(multiple_welch_c_result) > 0.1]<- "ns"
          multiple_welch_c_result<- c("",multiple_welch_c_result)
          multiple_welch_c_df[,2]<- multiple_welch_c_result
          return(multiple_welch_c_df) 
        } else if (input$variant_of_statistics == "multiple Welch's t-test (vs all)") {
          multiple_welch_all_df<- data.frame(colnames(mean_sd_n_group_dataframe), rep("",(length(colnames(mean_sd_n_group_dataframe)))), stringsAsFactors=FALSE)
          colnames(multiple_welch_all_df) <- c("group_name", "significance")
          combinations<- combn(unique(simulated_data[[1]]), 2)
          combinations_vector<- NA
          for (i in 1:ncol(combinations)) {
            combinations_vector[i]<- paste(combinations[1,i], combinations[2,i], sep="-")
          }
          multiple_welch_all_result<- NA
          for (i in 1:ncol(combinations)) {
            t_test<- t.test(simulated_data[[2]][simulated_data[[1]] == as.character(combinations[1,i])],simulated_data[[2]][simulated_data[[1]] == as.character(combinations[2,i])])
            multiple_welch_all_result[i]<- t_test$p.value
          }
          names(multiple_welch_all_result)<- combinations_vector
          multiple_welch_all_letters<- multcompLetters(multiple_welch_all_result, threshold = 0.1)
          multiple_welch_all_letters<- multiple_welch_all_letters$`Letters`
          multiple_welch_all_df[,2]<- multiple_welch_all_letters
          return(multiple_welch_all_df) 
        }
      }
    }
  })
  
  #Information which test was used for the calculation
  test_used<- reactive({
    if(nrow(datatable_groups()) < 6) {
      test_in_use <- "Statistical analysis wasn't done."
      return(test_in_use)
    } else if (!is.element(input$control_group_name,datatable_groups()$group_name)) {
      test_in_use <- "Statistical analysis wasn't done."
      return(test_in_use)
    } else if (all(as.data.frame(table(datatable_groups()$group_name))[[2]] >= 3, na.rm = TRUE) == FALSE) {
      test_in_use <- "Statistical analysis wasn't done."
      return(test_in_use)
    } else if (length(unique(datatable_groups()$group_name)) < 2) {
      test_in_use <- "Statistical analysis wasn't done."
      return(test_in_use)
    } else if (any(datatable_groups()$group_name == "")) {
      test_in_use <- "Statistical analysis wasn't done."
      return(test_in_use)
    } else {
      if (length(unique(datatable_groups()$group_name)) < 3) {
        test_in_use<- "Welch's t-test for 2 samples"
        return(test_in_use)
      } else {
        test_in_use<- input$variant_of_statistics
        return(test_in_use)
      }
    }
  })
  
  #Rewriting datatable_groups and datatable_significance for the purpose of the boxplot (include \n at the ends of the names of groups)
  significance_groups_renamed<- reactive({
    groups_renamed_vector<- NA
    for (i in 1:nrow(datatable_significance())) {
      groups_renamed_vector[i]<- paste(datatable_significance()$group_name[i], datatable_significance()$significance[i], sep = "\n")
    }
    return(groups_renamed_vector)
  })
  
  datatable_groups_renamed<- reactive({
    df_groups<- datatable_groups()
    significance_names_complete_vector<- datatable_significance()$group_name
    for (i in 1:length(significance_groups_renamed())) {
      df_groups$group_name[df_groups$group_name == significance_names_complete_vector[i]] <- significance_groups_renamed()[i]
    }
    df_groups$group_name<- factor(df_groups$group_name, levels = significance_groups_renamed())
    return(df_groups)
  })
  
  #Rewriting datatable_groups and densityplot_datatable for the purpose of the violinplot (problem with using geom_point with facet_wrap solved by making df with same rows by copying values as many)
  datatable_groups_density_renamed<- reactive({
    vector_density<- densityplot_datatable()$sample
    for (i in 1:length(datatable_groups_renamed()$group_name)) {
      vector_density[vector_density == datatable_groups_renamed()$sample_number[i]]<- as.vector(datatable_groups_renamed()$group_name[i])
    }
    vector_density<- as.factor(vector_density)
    df_groups_density_renamed<- data.frame(densityplot_datatable()$sample, vector_density)
    colnames(df_groups_density_renamed)<- c("sample_number", "group_name")
    df_groups_density_renamed$group_name<- factor(df_groups_density_renamed$group_name, levels = significance_groups_renamed())
    return(df_groups_density_renamed)
  })
  
  densityplot_groups_datatable<- reactive({
    vector_median<- rep(NA, length(densityplot_datatable()$sample))
    vector_stquartile<- rep(NA, length(densityplot_datatable()$sample))
    vector_rdquartile<- rep(NA, length(densityplot_datatable()$sample))
    for (i in 1:length(final_ggplot2_boxplot_datatable_median()$variable)) {
      vector_median[densityplot_datatable()$sample == datatable_groups_renamed()$sample_number[i]][1]<- final_ggplot2_boxplot_datatable_median()$value[i]
      vector_stquartile[densityplot_datatable()$sample == datatable_groups_renamed()$sample_number[i]][1]<- final_ggplot2_boxplot_datatable_stquartile()$value[i]
      vector_rdquartile[densityplot_datatable()$sample == datatable_groups_renamed()$sample_number[i]][1]<- final_ggplot2_boxplot_datatable_rdquartile()$value[i]
    }
    vector_median<- as.numeric(vector_median)
    vector_stquartile<- as.numeric(vector_stquartile)
    vector_rdquartile<- as.numeric(vector_rdquartile)
    df_densityplot_groups<- data.frame(vector_median, vector_stquartile, vector_rdquartile)
    colnames(df_densityplot_groups)<- c("median", "stquartile","rdquartile")
    df_densityplot_groups<- cbind(densityplot_datatable(),df_densityplot_groups)
    return(df_densityplot_groups)
  })
  
})
