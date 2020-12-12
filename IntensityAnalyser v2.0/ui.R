library('shinyWidgets')
library('shinyalert')
library('shinyBS')

jscode <- '
shinyjs.init = function() {
$(".nav").on("click", ".disabled", function (e) {
e.preventDefault();
return false;
});
}
'

css <- '
.disabled {
background: #eee !important;
cursor: not-allowed !important;
color: grey !important;
}
'

tagList(shinyjs::useShinyjs(), 
        tags$style(type="text/css", "body {padding-top: 70px;}",".navbar-default .navbar-brand {
                         color: #337ab7;}"),
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }"
        ),
        tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 51px;
               left: 0px;
               width: 100%;
               padding: 2px 0px 2px 0px;
               text-align: center;
               font-size: 60%;
               color: #808080;
               background-color: #EEEEEE;
               z-index: 105;
             }
          "),
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Busy...",id="loadmessage")),
        navbarPage(collapsible = TRUE, position = "fixed-top","IntensityAnalyser", id = "tabs",
           tabPanel("Upload file", value = "upload_file",
                    div(style="position:relative; margin-left: auto; margin-right: auto;",
                                  img(src="backround_tab1.png"),
                        div(style="text-align: left; position: absolute; top: 145px; fixed = TRUE; white-space: nowrap;",
                          HTML(
                            paste(
                              h5(tags$b("Upload file in a same format as example below:")),
                              h6(tags$b("1st column"), " – pixel number"),
                              h6(tags$b("2nd column"), " – marker intensities"),
                              h6(tags$b("3rd to the last column"), " – intensities of samples"),
                              h6(tags$b("1st row"), " - names of columns"),
                              h6(tags$b("2nd row"), " - pixel positions of the samples in the TRF scan")
                                 )
                              )
                           ),
                        div(style="text-align: left; position: absolute; top: 0px; fixed = TRUE; white-space: nowrap;",
                            fileInput("uploadfile", "Upload file (.xls, .xlsx), you want to analyze.",
                                      accept=c(".xlsx",".xls"), width = 320),
                            uiOutput("ui.run_program")
                            )
                    )
              ),

           tabPanel("Marker", value = "marker",
              tags$head(tags$style(HTML('#sidebar {width: 240px;}'))),
              div(style="position: absolute;",
                  numericInput("sens", "Sensitivity:", 15, width = 76)
              ),
              div(style="position: absolute;top:242px",
              selectInput("select_marker", "Select marker", width = 266, 
                          choices = list("GeneRuler™ 1kb Ladder" = 1, "GeneRuler™ 1kb Plus Ladder" = 2,
                                         "GeneRuler™ 1kb + GeneRuler™ 1kb Plus Ladder" = 3, "GeneRuler™ High Range Ladder" = 4,
                                         "GeneRuler™ High Range + GeneRuler™ 1kb Ladder" = 5, "GeneRuler™ High Range + GeneRuler™ 1kb Plus Ladder" = 6,
                                         "Lambda PFG Ladder (NEB)" = 7, "Lambda DNA Hind III Digest (NEB)" = 8
                                         ), selected = 1)
              ),
              div(style="position: absolute; left: 99px;top:86px",
                  HTML(
                    paste(
                      h6("In case of too many peak maxima,")
                    )
                  )
              ),
              div(style="position: absolute; left: 99px;top:104px",
                  HTML(
                    paste(
                      h6("increase the value and vice versa.")
                    )
                  )
              ),
              div(style="position: absolute;top: 309px;",
                  sidebarPanel(id = "sidebar",
                rhandsontable::rHandsontableOutput("markertable")
                  )
              ),
              div(style="position: absolute;left: 300px; width: 1040px;",
                plotOutput("maxima_plot")
              ),
              div(style="position: absolute; top: 470px; left: 306px;",
                  HTML(
                    paste(
                      h5(tags$b("Don't know what to do? Follow the instructions below:"))
                    )
                  )
              ),
              div(style="position: absolute; top: 505px; left: 306px; fixed = TRUE;",
                img(src="help_marker.png")
              ),
              div(style="position: absolute;left: 1380px; fixed = TRUE;",
                  img(src="warning_marker.png")
              ),
              div(style="text-align: left; position: absolute; top: 149px;left:29px; fixed = TRUE; white-space: nowrap;",
                  uiOutput("ui.proceed_without_correction")
              ),
              div(style="text-align: left; position: absolute; top: 193px;left:39px; fixed = TRUE; white-space: nowrap;",
                  uiOutput("ui.proceed_with_correction")
              )
            ),
           
           tabPanel("Marker correction", value = "marker_correction",
                    div(style="position: absolute;",
                        numericInput("sens_second_marker", "Sensitivity:", 15, width = 76)
                    ),
                    div(style="position: absolute; left: 99px;top:86px",
                        HTML(
                          paste(
                            h6("In case of too many peak maxima,")
                          )
                        )
                    ),
                    div(style="position: absolute; left: 99px;top:104px",
                        HTML(
                          paste(
                            h6("increase the value and vice versa.")
                          )
                        )
                    ),
                    div(style="position: absolute;top:144px",
                        uiOutput("slider_second_marker_option")
                    ),
                    div(style="position: absolute;top:183px;left: 335px;",
                        actionButton("forward_button_slider_for_second_marker", ">", style='padding:3px')
                    ),
                    div(style="position: absolute;top:183px;left: 316px;",
                        actionButton("backward_button_slider_for_second_marker", "<", style='padding:3px')
                    ),
                    div(style="position: absolute;left: 370px; width: 1040px;",
                        plotOutput("maxima_plot_second_marker")
                    ),
                    div(style="position: absolute;top:473px;left: 374px;",
                        HTML(
                          paste(
                            h5(tags$b("Graph of the first marker:"))
                          )
                        )
                    ),
                    div(style="position: absolute;top:512px;left: 370px; width: 1040px;",
                        plotOutput("maxima_plot_for_marker_tab")
                    ),
                    div(style="position: absolute;top:232px;",
                        actionButton("info_marker_correction_tab", icon = icon("info-circle", lib = "font-awesome"), "Don't know what to do?")
                    ),
                    div(style="position: absolute;top:232px;left: 273px;",
                        actionButton("proceed_after_second_marker", "Proceed", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                    ),
                    div(style="position: absolute;top: 286px;",
                        tags$head(tags$style(HTML('#sidebar_marker_correction {width: 306px;}'))),
                        sidebarPanel(id = "sidebar_marker_correction",
                                     rhandsontable::rHandsontableOutput("second_marker_table", width = 300)
                        )
                    )
           ),

           tabPanel("pixel/MW ratio", value = "pxmw_ratio",
                    div(style="position: absolute;left: 325px; width: 965px;",
                        plotOutput("polynom_plot")
                    ),
                    uiOutput("polynom_number_ui"),
                    div(style="position: relative; top: -10px",
                        htmlOutput("polynom_recommendation")
                    ),
                    tags$head(tags$style(HTML('#sidebar_anova {width: 244px;}'))),
                        sidebarPanel(id = "sidebar_anova",
                      tableOutput("dataframe_anova_test")
                    ),
                    div(style="position: absolute;top:532px;left: 332px;",
                      HTML(
                        paste(
                          h5(tags$b("Don't know what you should see?"))
                        )
                      )
                    ),
                    div(style="position: absolute; top: 567px;left: 332px; fixed = TRUE;",
                        img(src="help_polynom_regression.png")
                    ),
                    div(style="text-align: left; position: absolute; top: 95px;left:221px; fixed = TRUE; white-space: nowrap;",
                        uiOutput("ui.choose")
                    ),
                    div(style="position: absolute;top:539px;left:934px",
                        uiOutput("slider_pxmw_tab")
                    ),
                    div(style="position: absolute;top:578px;left: 1269px;",
                        uiOutput("forward_button_slider_pxmw_tab")
                    ),
                    div(style="position: absolute;top:578px;left: 1250px;",
                        uiOutput("backward_button_slider_pxmw_tab")
                    )
           ),
           
           tabPanel("Length calc.", value = "length_calc",
                    div(style="position: absolute;left:290px;width: 911px;top:50px",
                        plotOutput("sample_coordinates_plot", brush = brushOpts(id = "brush_sample_coordinates_plot", direction = "x", resetOnNew = TRUE))
                    ),
                    div(style="position: absolute;left: 4px;top: 505px;width: 298px;",
                        plotOutput("sample_coordinates_plot_selection")
                    ),
                    tags$style(type = "text/css", ".irs-grid-pol.small {height: 0px;}"),
                    uiOutput("slider"),
                    div(style="position: absolute;top:109px;left: 256px;",
                        actionButton("forward_button_slider", ">", style='padding:3px')
                    ),
                    div(style="position: absolute;top:109px;left: 237px;",
                        actionButton("backward_button_slider", "<", style='padding:3px')
                    ),
                    div(style="position: absolute;top:70px;left:1223px;",
                        uiOutput("threshold")
                       ),
                    div(style="position: absolute;top:324px;left: 164px;",
                        actionButton("save_sample_coordinates", "Manual Save", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                    ),
                    div(style="position: absolute;top:368px;left: 164px;",
                        actionButton("save_sample_coordinates_auto", "Auto Save", style="color: #fff; background-color: #ff4d4d; border-color: #c23a3a")
                    ),
                    div(style="position: absolute;top:412px;left: 164px;",
                        actionButton("restart_sample_coordinates", "Restart")
                    ),
                    div(style="position: absolute;top:412px;",
                        actionButton("restart_all_coordinates", "Abort Analysis", style="color: #fff; background-color: #ff4d4d; border-color: #c23a3a", width = 115)
                    ),
                    div(style="position: absolute;top:158px;",
                        HTML(
                          paste(
                            h5(tags$b("Specification of the first peak area:"))
                          )
                        )
                    ),
                    div(style="position: absolute;top:193px;left:182px",
                        actionButton("save_sample_coordinates_left_side", "Save")
                    ),
                    div(style="position: absolute;top:188px;left:30px",
                        htmlOutput("coordinates_left_peak_left_side")
                    ),
                    div(style="position: absolute;top:209px;left:30px",
                        htmlOutput("coordinates_left_peak_right_side")
                    ),
                    div(style="position: absolute;top:235px;",
                        HTML(
                          paste(
                            h5(tags$b("Specification of the last peak area:"))
                          )
                        )
                    ),
                    div(style="position: absolute;top:270px;left:182px",
                        actionButton("save_sample_coordinates_right_side", "Save")
                    ),
                    div(style="position: absolute;top:265px;left:30px",
                        htmlOutput("coordinates_right_peak_left_side")
                    ),
                    div(style="position: absolute;top:286px;left:30px",
                        htmlOutput("coordinates_right_peak_right_side")
                    ),
                    div(style="position: absolute;top:324px",
                        uiOutput("background_correction_button")
                    ),
                    div(style="position: absolute;top:368px;",
                        actionButton("restart", "Remove BGC", width = 115)
                    ),
                    div(style="position: relative;top:300px",
                        htmlOutput("info_sample")
                    ),
                    div(style="position: absolute;top:502px;",
                        HTML(
                          paste(
                            h5(tags$b("Your selected area:"))
                          )
                        )
                    ),
                    div(style="position: absolute;left:365px;top:531px;fixed = TRUE;",
                        img(src="help_length_calc.png")
                    ),
                    div(style="position: absolute;left:345px;top:491px;",
                        HTML(
                          paste(
                            h5(tags$b("Don't know what to do? Follow the instructions below:"))
                          )
                        )
                    ),
                    div(style="position: absolute;top:157px;left:1223px;",
                        img(src="explanation_button.png")
                    )
           ),
           tabPanel("Graphic result", value = "graphic_result",
                    div(style="position: absolute;top:309px",
                        prettyRadioButtons(inputId = "box_or_densityplots",
                                           label = "Plot options:",
                                           choices = c("Boxplot", "Violin plot"))
                    ),
                    div(style="position: absolute;top:389px;",
                        uiOutput("explanation_graphicresult_ui")
                    ),
                    div(style="position: absolute;top:70px",
                        numericInput("height", label = "Height (px)", value = 638, width = 77)
                    ),
                    div(style="position: absolute;top:70px;left:111px",
                        uiOutput("width_input")
                    ),
                    div(style="position: absolute;top:143px",
                        uiOutput("min_input")
                    ),
                    div(style="position: absolute;top:143px;left:111px",
                        uiOutput("max_input")
                    ),
                    div(style="position: absolute;top:216px",
                        sliderTextInput(
                          inputId = "yaxisrange",
                          label = "Labels of y-axis:",
                          choices = c("Auto","250","500","750","1000","2.5K"),
                          selected = "Auto",
                          grid = TRUE,
                          hide_min_max = TRUE,
                          width = 173
                        )
                    ),
                    div(style="position: absolute;left:43px;top:605px",
                    downloadButton("report", label = "Save report")
                    ),
                    div(style="position: absolute;top:143px;left:217px",
                        img(src="line.png")
                    ),
                    div(style="position: absolute;top:75px;left:246px",
                        actionButton("info_statistics", "How is statistics calculated?", icon = icon("info-circle", lib = "font-awesome"), width = 241)
                    ),
                    div(style="position: absolute;top:297px;left:232px",
                        sidebarPanel(id = "sidebar",
                                     rhandsontable::rHandsontableOutput("grouptable")
                        )
                    ),
                    useShinyalert(),
                    div(style="position: absolute;top:629px;left:247px;",
                        uiOutput("evaluate_button")
                    ),
                    div(style="position: absolute;top:629px;left:372px",
                        actionButton("restart_statistical_evaluation", "Cancel", width = 115, style="color: #fff; background-color: #ff4d4d; border-color: #c23a3a")
                    ),
                    div(style="position: absolute;top:233px;left:246px",
                    textInput("control_group_name", label = "", value = "Enter Control group name...", width = 241)
                    ),
                    div(style="position: absolute;top:124px;left:246px",
                        prettyRadioButtons(inputId = "variant_of_statistics",
                                           label = "Test options:",
                                           choices = c("multiple Welch's t-test (vs control)","multiple Welch's t-test (vs all)","Tamhane-Dunnett (vs control)", "Games-Howell (vs all)"))
                    ),
                    div(style="position: absolute;left:510px",
                        uiOutput("result")
                    ),
                    div(style="position: absolute;left:38px;top:651px",
                        actionButton("show_300DPI", "Show hi-res plot")
                    ),
                    tags$head(tags$style(".modal-dialog{width:4000px}")),
                    bsModal("modalplot", "Press Esc to get back", "show_300DPI", uiOutput("result_300DPI"))
                    
           )
  ),
  tags$style(css),
  tags$head(tags$style(HTML('.navbar-nav a {cursor: default}')))
)