library(shiny)
ui <- fluidPage(titlePanel("Shiny App to compute fixed/share/private SNP count"),
                tabsetPanel(
                  tabPanel("Input and Results", fluid = TRUE,
                           sidebarLayout(
                             
                             sidebarPanel(helpText(tags$h3("Data format"),"Your data must have a specific format, a tab separated
                 table with different columns: ",br(),"scaffold, position, reference_nucleotide, 
                                      alternative_nucleotide, allele_count and number_of_sample"),
                                          "To produce this kind of file from a vcf.gz you can use the following bash command:",br(),
                                          "bcftools view -S population1_individual_list vcf_file.vcf.gz |bcftools query  -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AC\\t%AN\\n' > population1_summary.txt",br(),
                                          "The word before _summary.txt will be used as a name for each population so name your files accordingly",
                                          fileInput(
                                            inputId = "files", 
                                            label = "Choose all summary files to upload at once", 
                                            multiple = TRUE
                                          )
                             ),
                             mainPanel(fluidPage(helpText("Once you uplodaded the data, it may take several minutes to show results just under this text")#,
                                        downloadButton('downloadtable',"Download the summary table"),
                   fluidRow(column(7,dataTableOutput('dtotable')))
                                       )
                                   
                             )
                           )
                  ),
                  tabPanel("Other informations", fluid = TRUE,
                           tags$h1("How SNP are computed ?"), 
                           "within population :",br(),
                           "There is two information when we compare alleles from population A and B ",br(),
                           "The SNP is said to be \"fixed\" if all samples from a population have the same allele",br(),
                           "Between population:",br(),
                           "The SNP is said \"snp\" if at least one sample have a different allele within a population",br(),
                           "When we compare alleles present in both populations, they can have same allele content (same ref or alt or same combination of multiple alleles) then we declare them \"same\" ",br(),
                           "If allele content is different between populations, then they are declared \"diff\" ",
                           tags$h2("We can obtain these combinations"),
                           "fixed fixed	same	=>IDENTIC",br(),
                           "snp snp	same	=>SHARED",br(),
                           "fixed snp	same	=>NOT POSSIBLE",br(),
                           "snp fixed	same	=>NOT POSSIBLE",br(),
                           
                           "fixed fixed	diff	=>FIXED",br(),
                           "snp snp	diff	=>SHARED",br(),
                           "fixed snp	diff	=>PRIVATE",br(),
                           "snp fixed	diff	=>PRIVATE"
                  ),
                  #Add a tab to download the data
                  tabPanel("Download the computed dataframe",fluidPage(helpText("You won't be able to see and download the data until the input is computed"),
                    # This one is linked by the id 'download'
                    downloadButton('download',"Download the data"),
                    fluidRow(column(7,dataTableOutput('dto')))
                  ))
                  
                  
                ))


server <- function(input, output) {
  packages <- c("dplyr","stringr","tidyr")
  install.packages(setdiff(packages, rownames(installed.packages())))    
  library(stringr)
  library(dplyr)
  library(tidyr)
  
  
  options(shiny.maxRequestSize=3000*1024^2) 
  output$contents <- renderTable({
    req(head(input$files))
    
    upload = list()
    for(nr in 1:length(input$files[, 1])){
      upload[[nr]] <-  read.table(file = input$files[[nr, 'datapath']],
                                  header = F,sep = "\t",colClasses=c("character"))
    }
    names(upload)<-str_replace( input$files$name,pattern = "_summary.txt",replacement = "")# Get names of each file, supposed to be the population name
    df<-do.call(gdata::combine,upload)# Combine all file in a single dataframe with an IDpop column
    colnames(df)=c("scaffold","position","ref","alt","sites","individual_known","source")#Rename columns
    number_of_alt=str_split(df$sites,pattern = ",",simplify = T)#Need the total number of alt allele
    number_of_alt=apply(number_of_alt, 2,as.numeric)#As numeric the character matrix
    number_of_alt=apply(number_of_alt, 1,sum,na.rm=T)# Sum count of different alt allele
    df$number_of_ref=as.numeric(df$individual_known)-number_of_alt#Substract number of GT know to get the count of ref
    df$sites=paste(df$number_of_ref,df$sites,sep = ",")#Add ref count in sites column
    df=df%>%group_by(scaffold,position)%>%filter(!any(individual_known==0))#Remove SNP that are unknown in any population (dots in vcf)
    df$sites_fixed=str_replace_all(string = df$sites,pattern = "[1-9]+",replacement = "1")# Replace number higher than 0 by one (except multiple of 10)
    df$sites_fixed=str_replace_all(string = df$sites_fixed,pattern = "[1-9]0",replacement = "1")#Replace multiple of 10 by 1
    df$sites_fixed=str_replace_all(pattern = ",",replacement = "",df$sites_fixed)#Remove comma to get a code with 0 and 1
    
    
    # Get status of each snp --------------------------------------------------
    df=df%>%ungroup()%>%mutate(fixed=if_else(str_count(sites_fixed,pattern = "1")>1,"snp","fixed"))#If there is more than one genotype per population, consider it as a snp and not fixed in the population
    df$source=as.character(df$source)
    #Making pairwise population for each snp in a way that column.x is compared to column.y
    df=inner_join(df,df,by = c("scaffold"="scaffold","position"="position")) %>% filter(source.x>source.y) 
    df=df%>%mutate(diff=if_else(sites_fixed.x==sites_fixed.y,"same","diff"))#Tells if the allele distribution is the same in population x and y 
    #Compute the status of each snp base on if it is fixed or not and the same in both population or not
    df=df%>%mutate(status=case_when(diff=="diff"&fixed.x=="fixed"&fixed.y=="fixed"~"fixed",
                                    fixed.x=="snp"&fixed.y=="snp"~"shared",
                                    (fixed.x=="snp"&fixed.y=="fixed")|(fixed.x=="fixed"&fixed.y=="snp")~"private",
                                    diff=="same"&fixed.x=="fixed"&fixed.y=="fixed"~"Identical"
    ))
    df$pairwise=paste(df$source.x,"vs",df$source.y) # Get the name of the pairwise comparison
    
    #Prepare the data for download
    thedata <- reactive(df)
    output$dto <- renderDataTable({thedata()})
    output$download <- downloadHandler(
      filename = function(){"Whole_data.csv"}, 
      content = function(fname){
        write.csv(thedata(), fname)
      }
    )
    
       #Prepare the summary table for download
    thetable <- reactive(df%>%
             group_by(pairwise,status)%>%
             summarise(n=n())%>%
             spread(status, n))
    output$dtotable <- renderDataTable({thetable()})
    output$downloadtable <- downloadHandler(
      filename = function(){"Summary_table.csv"}, 
      content = function(fname){
        write.csv(thetable(), fname)
      }
    )
    
    
    return(df%>%
             group_by(pairwise,status)%>%
             summarise(n=n())%>%
             spread(status, n))
  })
}
shinyApp(ui, server)
