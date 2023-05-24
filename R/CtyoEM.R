#' Make an enrichment map from exported GSEA results
#'
#' Use RCy3 to make an enrichment map with positive and negative enrichment results for a single contrast
#'
#' @param GMTfile character with complete path to gmt file used for GSEA analyses
#' @param pos character with complete path to tsv file with positive enrichment scores from GSEA
#' @param neg character with complete path to tsv file with negative enrichment scores from GSEA
#' @param title character indicating title for enrichment map
#' @param parseBADER Boolean indicating if Bader Lab pathway sets were used. Default is TRUE.
#' @param PDFpre character indicating prefix/path for PDF file (if wanted). If not wanted, leave NULL.
#' @param AUTOannotate Boolean indicating if AutoAnnotate summary network should also be created. Default is FALSE.
#' @param ANNOtitle character indidcating AutoAnnotate summary network title if AUTOannotate = TRUE.
#' @return Network SUID. Two SUIDs for original and summary network, if AUTOannotate = TRUE.
#' @export

CytoEM<-function(GMTfile, pos, neg, title, parseBADER = TRUE, PDFpre = NULL, AUTOannotate = FALSE, ANNOtitle = paste0(title,"_summary")){
  require(RCy3, quietly = TRUE)
  pb<-TRUE
  if(isFALSE(parseBADER)){
    pb<-FALSE
  }
  em_command<-paste('enrichmentmap build analysisType="GSEA"',
                    "gmtFile=",GMTfile,
                    'pvalue=',"0.05",'qvalue=',"0.1",
                    'similaritycutoff=',"0.5",
                    'coefficients=',"OVERLAP",
                    'enrichmentsDataset1=',pos,
                    'enrichments2Dataset1=',neg,
                    'parseBaderlabNames=',pb,
                    sep=" ")
  nw<-commandsRun(em_command)

  renameNetwork(title, network=as.numeric(nw))
  if(isFALSE(AUTOannotate)){
    if(!is.null(PDFpre)){
      message("Exporting network to PDF: ",paste0(PDFpre,title,".PDF"))
      exportImage(paste0(PDFpre, title,".PDF"), type = "PDF")
    }
    return(nw)
  } else {
    nodetable_colnames<-getTableColumnNames(table="node", network=as.numeric(nw))
    descr_attrib<-nodetable_colnames[grep(nodetable_colnames, pattern="GS_DESCR")]
    autoannotate_url<-paste0("autoannotate annotate-clusterBoosted labelColumn=",descr_attrib," maxWords=3 ")
    current_name<-commandsGET(autoannotate_url)
    #layoutNetwork("cose", network=as.numeric(getNetworkSuid(as.numeric(nw))))
    commandsGET("autoannotate layout layout=cose network='current'")
    if(!is.null(PDFpre)){
      message("Exporting network to PDF: ",paste0(PDFpre,title,".PDF"))
      exportImage(paste0(PDFpre, title,".PDF"), type = "PDF")
    }
    commandsGET("autoannotate summary network='current'")
    summary_network_suid<-getNetworkSuid()
    renameNetwork(title=ANNOtitle,
                  network=as.numeric(summary_network_suid))

    summary_nodes<-getTableColumns(table="node", columns=c("name"))

    #Change node sizes to match number of gene sets in cluster - currently error in function if not done one at a time:
    for(i in rownames(summary_nodes)){
    clearNodePropertyBypass(i, visual.property="NODE_SIZE", network=as.numeric(summary_network_suid))
    }
    setNodeSizeMapping(table.column = "cluster node count", mapping.type="c", sizes=c(60,150),
                       style.name="EM1_Visual_Style", network=as.numeric(summary_network_suid))
    layoutNetwork("genemania-force-directed", network=as.numeric(summary_network_suid))
    if(!is.null(PDFpre)){
      message("Exporting summary network to PDF: ",paste0(PDFpre,title,"_summary.PDF"))
      exportImage(paste0(PDFpre, title,"_summary.PDF"), type = "PDF", network=as.numeric(summary_network_suid))
    }
    return(list(Original_Network=as.numeric(nw), Summary_Network=summary_network_suid))
  }
}
