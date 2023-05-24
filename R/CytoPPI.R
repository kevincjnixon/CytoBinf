#' Create a Protein-protein interaction network
#'
#' Generate a PPI network with either geneMANIA or STRING databases using
#' RNA-seq data. Nodes (genes) will be coloured by mean expression (if count table is provided) or log2FoldChange (if results table is provided).
#'
#' @param data data.frame of either normalized gene counts or DESEq2 results. Gene names must be rownames
#' @param datTYPE on of either "counts" or "res" indicating the kind of data found in data. Defaults to 'counts'.
#' @param genes string vector indicating gene names (matching rownames(data)) to which data should be subset. Leave as NULL to include all genes in table.
#' @param species string indicating organism used. Defaults to 'H. sapiens'
#' @param PPI one of either "gm" or "string" indicating which Cytoscape app to use. Default is "gm"
#' @param gmONLINE Boolean indicating if online or offline mode should be used for geneMANIA query. If query size >=1000, defaults to offline (FALSE), otherwise defaults to online (TRUE).
#' @param title string indicating name for network. Defaults to "PPI_Network".
#' @param cytoHubba Numeric indicating number of hub genes to search for using cytoHubba. Leave at 0 to skip this analysis.
#' @param PDF string indicating PDF file name for exporting generated network image. Leave NULL to not export
#' @return Network SUID
#' @export

CytoPPI<-function(data, datTYPE = "counts", genes = NULL, species = "H. sapiens", PPI = "gm", gmONLINE = TRUE, title = "PPI_Network", cytoHubba = 0, PDF = NULL){
  require(RCy3, quietly=T)
  data<-as.data.frame(data)
  if(!is.null(genes)){
    message("Subsetting data to listed genes...")
    data<-data[which(rownames(data) %in% genes),]
  }
  ppi<-NULL
  if(datTYPE=="counts"){
    message("Nodes will be coloured according to average gene expression levels")
    data$fill.col<-rowMeans(data)
    data$size.col<-rep(60, nrow(data))
    data$border.col<-rep("no", nrow(dat))
  }
  if(datTYPE=="res"){
    message("Nodes will be coloured according to log2FoldChange and sized according to average expression levels")
    message("Node borders will indicate significance")
    data$fill.col<-data$log2FoldChange
    data$size.col<-data$baseMean
    data$border.col<-ifelse(data$padj<0.05,"yes","no")
    data$border.col[which(is.na(data$border.col))]<-"no"
  }
  if(PPI=="gm"){
    message("Generating geneMANIA PPI network...")
    if(nrow(data)>=1000){
      gmONLINE<-FALSE
      message("More than 1000 genes in query, setting gmONLINE to FALSE for offline mode...")
    }
    query<-paste(rownames(data), collapse="|")
    offline<-'true'
    if(isTRUE(gmONLINE)){
      offline<-'false'
    }
    gm_command<-paste0("genemania search genes=",query," organism=",species," offline=",offline," geneLimit=0")
    ppi<-commandsRun(gm_command)
  }
  if(PPI=="string"){
    message("Generating STRINGdb PPI network...")
  }
  style.name<-"PPIstyle"
  defaults.list<-list(NODE_SHAPE="ellipse",
                      NODE_SIZE=40,
                      NODE_FILL_COLOR="#9F9F9F",
                      NODE_BORDER_PAINT="#FFFFFF",
                      NODE_BORDER_WIDTH=2.0)
  node.label.map<-mapVisualProperty('node label', 'gene name', 'p')
  createVisualStyle(style.name, defaults.list, list(node.label.map))
  setVisualStyle(style.name)
  #add data to node table
  loadTableData(data, data.key.column = "row.names", table="node", table.key.column="query term")
  newNodeTab<-getTableColumns()
  data.values<-c(min(newNodeTab$fill.col,na.rm=T),0,max(newNodeTab$fill.col,na.rm=T))
  node.colours<-c(rev(RColorBrewer::brewer.pal(length(data.values), "RdBu")))
  setNodeColorMapping("fill.col", data.values, node.colours, style.name=style.name)
  setNodeSizeMapping("size.col", mapping.type="c", sizes=c(40,100), style.name=style.name)
  setNodeBorderColorMapping("border.col", mapping.type="d", table.column.values=c("no","yes"), colors=c("black","yellow"), style.name=style.name)
  if(!is.null(PDF)){
    exportImage(PDF, type="PDF")
  }
  if(cytoHubba>0){
    stop("cytoHubba not supported by RCy3 yet...")
  }
  return(as.numeric(getNetworkSuid()))
}
