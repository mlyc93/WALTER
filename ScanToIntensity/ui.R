if (!require('shinyWidgets')) install.packages('shinyWidgets'); library('shinyWidgets')

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
        navbarPage(collapsible = TRUE, position = "fixed-top","ScanToIntensity", id = "tabs",
                   tabPanel("Upload picture",
                            div(style="position:relative; margin-left: auto; margin-right: auto;",
                                div(style="position: absolute;top: 153px;left:1px;",
                                    img(src="doesdonts_upload.png")
                                ),
                                div(style="text-align: left; position: absolute; top: 0px; fixed = TRUE; white-space: nowrap;",
                                    fileInput("uploadpicture", "Upload picture (.bmp, .jpg) you want to transform.",
                                              accept=c(".bmp", ".jpg", ".jpeg"), width = 320),
                                    uiOutput("ui.run_program")
                                )
                            )
                   ),
                   tabPanel("Sample transformation",
                            div(style="position: absolute;left:16px;top:75px; fixed = TRUE; white-space: nowrap;",
                                img(src="help_scantointensity.png")
                            ),
                            div(style="position: absolute;left:739px;top:73px",
                                uiOutput("ui_picture_plot")
                            ),
                            div(style="position: absolute;left:405px",
                                uiOutput("size_picture_slider")
                            ),
                            div(style="position: absolute;;left:660px;top:113px",
                                prettyCheckbox(inputId = "invert",
                                               label = "Invert",
                                               icon = icon("check"),
                                               status = "primary",
                                               inline = TRUE)
                            ),
                            div(style="position: absolute;left:405px;top:163px",
                                uiOutput("slider")
                            ),
                            div(style="position: absolute;left:405px;top: 256px",
                                actionButton("save_sample_coordinates", "Calculate profile")
                            ),
                            div(style="position: absolute;left: 572px;top: 256px",
                            downloadButton("downloadExcel", label = "Intensity profiles")
                            ),
                            div(style="position: absolute;left:405px;top:529px;",
                                htmlOutput("text_history_plot")
                            ),
                            div(style="position: absolute;left:405px;top:557px",
                                plotOutput("intensity_sample_plot_history")
                            ),
                            div(style="position: absolute;left:405px;top:298px;",
                                HTML(
                                  paste(
                                    h5(tags$b("Intensity profile of current sample:"))
                                  )
                                )
                            ),
                            div(style="position: absolute;left:405px;top:326px",
                                plotOutput("intensity_sample_plot")
                            )
                   )
        ),
        tags$style(css),
        tags$head(tags$style(HTML('.navbar-nav a {cursor: default}')))
)