#' Set up Cytoscape and Necessary Apps
#'
#' This function is to set up Cytoscape and all of the necessary apps for use
#' with the BinfTools package.
#' Currently. The necessary apps include:
#' - enrichmentmap
#' - clustermaker2
#' - autoannotate
#' - wordcloud
#' - stringapp
#' - aMatReader
#' - geneMANIA
#' - cytoHubba
#' Additional apps can be added through the app argument.
#'
#'
#' @param app String with name(s) of Cytoscape apps to install. Default is NULL.
#' @export

CytoSetUp<-function(app = NULL){
  require(RCy3, quietly = TRUE)
  #list of cytoscape apps to install
  installation_responses <- c()

  #list of app to install
  cyto_app_toinstall <- c("clustermaker2", "enrichmentmap", "autoannotate", "wordcloud", "stringapp", "aMatReader",
                          "geneMANIA","cytoHubba")

  cytoscape_version <- unlist(strsplit( cytoscapeVersionInfo()['cytoscapeVersion'],split = "\\."))
  if(length(cytoscape_version) == 3 && as.numeric(cytoscape_version[1]>=3)
     && as.numeric(cytoscape_version[2]>=7)){
    for(i in 1:length(cyto_app_toinstall)){
      #check to see if the app is installed.  Only install it if it hasn't been installed
      if(!grep(commandsGET(paste("apps status app=\"", cyto_app_toinstall[i],"\"", sep="")),
               pattern = "status: Installed")){
        installation_response <-commandsGET(paste("apps install app=\"",
                                                  cyto_app_toinstall[i],"\"", sep=""))
        installation_responses <- c(installation_responses,installation_response)
      } else{
        installation_responses <- c(installation_responses,"already installed")
      }
    }
    installation_summary <- data.frame(name = cyto_app_toinstall,
                                       status = installation_responses)

    knitr::kable(list(installation_summary),
                 booktabs = TRUE, caption = 'A Summary of automated app installation'
    )
  }
  cytoscapePing()
}
