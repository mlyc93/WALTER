library('shiny')
library('shinyjs')
library('readxl')
library('ggplot2')
library('dplyr')
library('rJava')
library('xlsx')
library('imager')
library('grid')
library('showtext')
library('showtextdb')
library('sysfonts')
library('curl')
library('gtools')

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
    shinyjs::disable(selector = ".navbar-nav a[data-value=upload_file]")
  })
  
  #Create dataframe out of picture
  intensity_datatable<- reactive({
    dataframe_picture_raw<- as.data.frame(fileimage())
  })

  #Horizontal value of the picture
  horizontal_pixel<- reactive({
    dataframe_picture_raw<- intensity_datatable()
    horizontal_pixel_max<- max(dataframe_picture_raw$x)
    return(horizontal_pixel_max)
  })
  
  #Vertical value of the picture
  vertical_pixel<- reactive({
    dataframe_picture_raw<- intensity_datatable()
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
  
  #Reactive values for slider information
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
  
  #Creating value of intensity profile dataframe for which cbind function will bind columns with intensity values for selected samples when saved
  intensitydataframe<- reactiveValues()
  intensitydataframe$intensity_profile_datatable<- reactiveVal({0})
  
  #Creating value of x_position profile dataframe for which cbind function will bind columns with x position values for selected samples when saved
  x_position<- reactiveValues()
  x_position$x_position_datatable<- reactiveVal({0})
  
  #Saving data after clicking action button
  observeEvent(input$save_sample_coordinates, {
    if(!is.null(input$brush_sample_coordinates_plot)) {
      coordinates <- data.frame(xmin=values$brush_sample_coordinates()$xmin, xmax=values$brush_sample_coordinates()$xmax, sample_nu=values$sample_slider(), counter=input$save_sample_coordinates)
      counter<- slider$vector()
      plot_history <- input$sample_slider
      if(input$sample_slider == last(counter)) {
        counter<- c(counter, last(counter) + 1)
      }
      slider$vector(counter)
      slider$plot_history(plot_history)
      
      #saving selected area as intensity profile to the reactive value intensitydataframe$intensity_profile_datatable
      if(last(slider$vector())-1 == 0) {
        intensity_profile_dataframe<- as.data.frame(matrix(seq(from = 0,to = vertical_pixel()-1,by = 1), nrow = vertical_pixel(), ncol = 1))
        colnames(intensity_profile_dataframe)[1]<- "pixel"
      } else {
        intensity_profile_dataframe<- intensitydataframe$intensity_profile_datatable()
        }
      intensity_profile_dataframe_sample<- as.data.frame(matrix(rep(0, vertical_pixel()), nrow = vertical_pixel(), ncol = 1))
      min_repeat_number<- if (round(as.numeric(coordinates$xmin[1]),digits = 0) == 0) {
        round(as.numeric(coordinates$xmin[1]),digits = 0)+1
      } else {
        round(as.numeric(coordinates$xmin[1]),digits = 0)
      }
      max_repeat_number<- round(as.numeric(coordinates$xmax[1]),digits = 0)
      for (j in min_repeat_number:max_repeat_number) {
        intensity_pixel_row<- intensity_datatable()$value[intensity_datatable()$x == j]
        intensity_profile_dataframe_sample[,j-(min_repeat_number-1)]<- intensity_pixel_row
      }
      intensity_profile_dataframe_sample<- as.data.frame(rowSums(intensity_profile_dataframe_sample))
      if (input$sample_slider == 0) {
        colnames(intensity_profile_dataframe_sample)[1]<- "marker"
      } else {
        colnames(intensity_profile_dataframe_sample)[1]<- as.character(input$sample_slider)
      }
      intensity_profile_dataframe<- as.data.frame(cbind(intensity_profile_dataframe,intensity_profile_dataframe_sample))
      intensitydataframe$intensity_profile_datatable(intensity_profile_dataframe)
      
      #saving selected area as x_position to the reactive value x_position$x_position_datatable
      if(last(slider$vector())-1 == 0) {
        x_position_profile_dataframe<- as.data.frame(0)
        colnames(x_position_profile_dataframe)[1]<- "pixel"
      } else {
        x_position_profile_dataframe<- x_position$x_position_datatable()
      }
      x_position_profile_dataframe_sample<- as.data.frame((round(as.numeric((max_repeat_number-min_repeat_number)/2+min_repeat_number), digits = 0)))
      if (input$sample_slider == 0) {
        colnames(x_position_profile_dataframe_sample)[1]<- "marker"
      } else {
        colnames(x_position_profile_dataframe_sample)[1]<- as.character(input$sample_slider)
      }
      intensity_profile_dataframe<- as.data.frame(cbind(x_position_profile_dataframe,x_position_profile_dataframe_sample))
      x_position$x_position_datatable(intensity_profile_dataframe)
    }
  })
  
  
  #Making datatable with all intensity profiles filtered from the duplicates leaving the last one and sorting it as it is required by the IntensityAnalyser tool
  intensity_profile_datatable<- reactive({
    intensity_profile_dataframe<- intensitydataframe$intensity_profile_datatable()
    intensity_profile_dataframe<- intensity_profile_dataframe[,!duplicated(colnames(intensity_profile_dataframe), fromLast = TRUE)]
    intensity_profile_dataframe<- intensity_profile_dataframe%>%select(pixel,marker,sort(names(.)))
    intensity_profile_dataframe<- intensity_profile_dataframe[mixedorder(colnames(intensity_profile_dataframe))]
    intensity_profile_dataframe<- intensity_profile_dataframe%>%select(pixel,marker,everything())
    return(intensity_profile_dataframe)
  })
  
  #Making datatable with all intensity profiles filtered from the duplicates leaving the last one and sorting it as it is required by the IntensityAnalyser tool
  x_position_datatable<- reactive({
    x_position_dataframe<- x_position$x_position_datatable()
    x_position_dataframe<- x_position_dataframe[,!duplicated(colnames(x_position_dataframe), fromLast = TRUE)]
    x_position_dataframe<- x_position_dataframe%>%select(pixel,marker,sort(names(.)))
    x_position_dataframe<- x_position_dataframe[mixedorder(colnames(x_position_dataframe))]
    x_position_dataframe<- x_position_dataframe%>%select(pixel,marker,everything())
    return(x_position_dataframe)
  })
  
  #Download the table
  output$downloadExcel <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_intensityprofiles", ".xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(rbind(x_position_datatable(),intensity_profile_datatable()), file, row.names = FALSE)
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