if (!require('shiny')) install.packages('shiny'); library('shiny')
if (!require('shinyjs')) install.packages('shinyjs'); library('shinyjs')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('rJava')) install.packages('rJava'); library('rJava')
if (!require('xlsx')) install.packages('xlsx'); library('xlsx')
if (!require('imager')) install.packages('imager'); library('imager')
if (!require('grid')) install.packages('grid'); library('grid')
if(!require('showtext')) install.packages('showtext'); library('showtext')
if(!require('showtextdb')) install.packages('showtextdb'); library('showtextdb')
if(!require('sysfonts')) install.packages('sysfonts'); library('sysfonts')
if(!require('curl')) install.packages('curl'); library('curl')

options(shiny.maxRequestSize = 100*1024^2)
shinyServer(function(input, output, session){
  shinyjs::disable(selector = '.navbar-nav a')
  #Data input
  fileimage<- reactive({
    infile <- input$uploadpicture
    if (is.null(infile)){
      return(NULL)      
    }
    image<- load.image(infile$datapath)
    if(dim(image)[4] == 3){
     image<- grayscale(image, method = "Luma", drop = TRUE)
    }
    return(image)
  })
  
  #Upload plot font
  sysfonts::font_add_google('Roboto', 'Roboto')
  showtext_auto()
  
  #Run program button for moving to Sample transformation and enabling other tabs
  output$ui.run_program <- renderUI({
    if (is.null(fileimage())) return()
    actionButton("run_program", "Run program", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  })
  
  observeEvent(input$run_program, {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Sample transformation")
    shinyjs::enable(selector = '.navbar-nav a')
  })
  
  #Create dataframe out of picture
  intensity_datatable<- reactive({
    dataframe_picture_raw<- as.data.frame(fileimage())
  })

  #Horizontal value of the picture
  horizontal_pixel<- reactive({
    dataframe_picture_raw<- as.data.frame(fileimage())
    horizontal_pixel_max<- max(dataframe_picture_raw$x)
    return(horizontal_pixel_max)
  })
  
  #Vertical value of the picture
  vertical_pixel<- reactive({
    dataframe_picture_raw<- as.data.frame(fileimage())
    vertical_pixel_max<- max(dataframe_picture_raw$y)
    return(vertical_pixel_max)
  })
  
  #Invert color option later for the plot
  invert_y_n<- reactive({
    if (input$invert == TRUE) {
      invert<- 1-fileimage()
      return(invert)
    } else {
      invert <- fileimage()
      return(invert)
    }
  })
  
  #Plot the picture
  output$picture<- renderPlot({
    horizontal_breaks<- seq(horizontal_pixel()/10, horizontal_pixel(), by = horizontal_pixel()/10)
    vertical_breaks<- seq(vertical_pixel()/10, vertical_pixel(), by = vertical_pixel()/10)
    ggplot(data.frame()) +
      geom_point() +
      annotation_custom(rasterGrob(invert_y_n(),
                                   width = unit(1,"npc"),
                                   height = unit(1,"npc")),
                        -Inf, Inf, -Inf, Inf) +
      scale_x_continuous(expand = c(0,0), limits = c(0, horizontal_pixel()), breaks = c(0,horizontal_breaks), position = "top") +
      scale_y_reverse(expand = c(0,0), limits = c(vertical_pixel(), 0), breaks = c(0,vertical_breaks)) +
      xlab("pixel") +
      ylab("pixel") +
      theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14))
  })
  
  output$size_picture_slider <- renderUI({
    sliderInput("picture_size", "Plot size:", min = 0, max = 300, value = (630/vertical_pixel())*100, step = 1, post = "%", width = 221)
  })
  output$ui_picture_plot<- renderUI({
    plotOutput("picture", click = "plot_click", brush = brushOpts(id = "brush_sample_coordinates_plot", direction = "x", resetOnNew = TRUE),
    width= horizontal_pixel()*(input$picture_size/100), height = vertical_pixel()*(input$picture_size/100))
    })
  
  #Reactive values
  slider<- reactiveValues()
  slider$vector<- reactiveVal({0})
  slider$plot_history <- reactiveVal({0})
  
  
  #Slider for number of samples
  output$slider <- renderUI({
    sliderInput("sample_slider", "Selected sample:", min=0, max=last(slider$vector()), value=last(slider$vector()), step = 1, width = 296)
  })
  
  #Making input$brush_sample_coordinates_plot a value
  values<-reactiveValues()
  observeEvent(input$brush_sample_coordinates_plot,{
    values$brush_sample_coordinates<-reactive({input$brush_sample_coordinates_plot})
  })
  
  #Making input$sample_slider+1 a value
  observeEvent(input$sample_slider+1, {
    values$sample_slider<- reactive({input$sample_slider+1})
  })
  
  #Saving data after clicking action button
  values$dataframe_coordinates <- reactiveVal({data.frame(xmin=NA, xmax=NA, sample_nu=NA, counter=NA)})
  
  observeEvent(input$save_sample_coordinates, {
    if(!is.null(input$brush_sample_coordinates_plot)) {
      old_values <- values$dataframe_coordinates()
      new_values <- 
        data.frame(xmin=values$brush_sample_coordinates()$xmin, xmax=values$brush_sample_coordinates()$xmax, sample_nu=values$sample_slider(), counter=input$save_sample_coordinates)
      new_dataframe <- rbind(old_values, new_values)
      values$dataframe_coordinates(new_dataframe)
      counter<- slider$vector()
      plot_history <- input$sample_slider
      if(input$sample_slider == last(counter)) {
        counter<- c(counter, last(counter) + 1)
      }
      slider$vector(counter)
      slider$plot_history(plot_history)
    }
  })
  
  #Select last measurements for each sample number from "memory" created above
  coordinates_from_memory<- reactive({
    coordinates_memory<- values$dataframe_coordinates()
    require(dplyr)
    coordinates_memory %>% group_by(sample_nu) %>% filter(counter %in% max(counter)) -> coordinates_memory
    coordinates_memory_ordered<- coordinates_memory[order(coordinates_memory$sample_nu),]
    return(coordinates_memory_ordered)
  })
  
  #Filtering picture datatable according to the selection for each sample and making rowSums, alas intensity profiles for each sample
  #from min repeat is 1 more beacause it is already valued in the first for loop
  intensity_profile_datatable<- reactive({
    intensity_profile_dataframe<- as.data.frame(matrix(seq(from = 0,to = vertical_pixel()-1,by = 1), nrow = vertical_pixel(), ncol = 1))
    colnames(intensity_profile_dataframe)[1]<- "pixel"
    for (i in 1:(nrow(coordinates_from_memory())-1)) { 
      intensity_profile_dataframe_sample<- as.data.frame(matrix(rep(0, vertical_pixel()), nrow = vertical_pixel(), ncol = 1))
      min_repeat_number<- if (round(as.numeric(coordinates_from_memory()$xmin[i]),digits = 0) == 0) {
        round(as.numeric(coordinates_from_memory()$xmin[i]),digits = 0)+1
      } else {
        round(as.numeric(coordinates_from_memory()$xmin[i]),digits = 0)
      }
      max_repeat_number<- round(as.numeric(coordinates_from_memory()$xmax[i]),digits = 0)
      for (j in min_repeat_number:max_repeat_number) {
        intensity_pixel_row<- intensity_datatable()$value[intensity_datatable()$x == j]
        intensity_profile_dataframe_sample[,j-(min_repeat_number-1)]<- intensity_pixel_row
      }
      intensity_profile_dataframe_sample<- rowSums(intensity_profile_dataframe_sample)
      intensity_profile_dataframe[,i+1]<- intensity_profile_dataframe_sample
      if (coordinates_from_memory()$sample_nu[i] == 1) {
        colnames(intensity_profile_dataframe)[i+1]<- "marker"
      } else {
      colnames(intensity_profile_dataframe)[i+1]<- coordinates_from_memory()$sample_nu[i]-1
      }
    }
   
    return(intensity_profile_dataframe)
  })
  
  #Download the table
  output$downloadExcel <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_intensityprofiles", ".xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(intensity_profile_datatable(), file, row.names = FALSE)
    }
  )
  
  #Create intensity profile plot based on selected area of the picture for each sample for last saved profile
  sample_graph_intensity_history<- reactive({
    vertical_breaks<- seq(vertical_pixel()/10, vertical_pixel(), by = vertical_pixel()/10)
    ggplot2::qplot(x= intensity_profile_datatable()$pixel, y= as.numeric(intensity_profile_datatable()[[slider$plot_history()+1+1]]), geom = "line", color = I("#337ab7")) +
      ylab("intensity") +
      xlab("pixel") +
      theme_light() +
      theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major=element_blank())
  })
  
  output$intensity_sample_plot_history<- renderPlot(sample_graph_intensity_history(), height = 200, width = 300)
  
  #Text for history plot
  output$text_history_plot<- renderText({
    paste("<b>","Your last saved profile (sample no. ", slider$plot_history(), "):","</b>", sep="")
  })
  
  #Create intensity profile plot based on selected area of the picture for each sample for current sample
  sample_graph_intensity<- reactive({
    vertical_breaks<- seq(vertical_pixel()/10, vertical_pixel(), by = vertical_pixel()/10)
    ggplot2::qplot(x= intensity_profile_datatable()$pixel, y= as.numeric(intensity_profile_datatable()[[input$sample_slider+1+1]]), geom = "line", color = I("#337ab7")) +
      ylab("intensity") +
      xlab("pixel") +
      theme_light() +
      theme(axis.text=element_text(family="Roboto"), axis.title=element_text(family="Roboto", size=14), legend.position= "none", panel.grid.minor= element_blank(), panel.grid.major=element_blank())
  })
  
  output$intensity_sample_plot<- renderPlot(sample_graph_intensity(), height = 200, width = 300)
  
})