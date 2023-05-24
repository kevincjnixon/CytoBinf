#' Basic Network from node/edge tables
#'
#' Create a basic network using existing node/egde tables
#'
#' @param node data.frame with at least one column of strings named 'id'
#' @param edge data.frame with at least three columns of strings named 'source','target', and 'interaction' where 'source' and 'target' match node 'ids' from the node table.
#' @param title character indicating network name
#' @param collection character indicating network collection
#' @param png character indicating exported png file. Default is NULL if more edits need to be done before exporting.
#' @return Network SUID
#' @export

DFNetwork<-function(node, edge, title, collection, png = NULL){
  require(RCy3, quietly = TRUE)
  createNetworkFromDataFrames(nodes = node, edges = edge, title = title, collection = collection)
  if(!is.null(png)){
    exportImage(png, type = "png")
  }
  return(as.numeric(getNetworkSuid()))
}
