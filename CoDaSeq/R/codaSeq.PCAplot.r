#' Visualise and customise a PCA biplot
#'
#' This function takes an object produced by `prcomp` and produces a biplot which
#' can be customise to show sample or loading groups, plot samples or loadings as
#' text or data points, include density plots above and beside the biplot for
#' sample or loading groups, and add a legend to the biplot.
#'
#' @param pcx An object produced by `prcomp`. The input to `prcomp` must have
#'   samples in rows and features in columns.
#' @param plot.groups Logical value indicating whether or not to colour samples
#'   (i.e. `pcx$x`) by group membership. Will be drawn as text or symbols
#'   depending on `grp.sym` (default: `FALSE`).
#' @param plot.loadings Logical value indicating whether or not to plot loadings
#'   (i.e. `pcx$rotations`). Will be drawn as text or symbols depending on
#'   `loadings.sym` (default: `TRUE`).
#' @param plot.ellipses Character vector indicating whether or not to draw
#'   confidence ellipses around sample *or* loading groups. Must be one of either
#'   `"groups"` or `"loadings"`, or `NULL`.
#' @param plot.density Logical value indicating whether or not to draw density
#'   plots for the specified principal components above and beside the main
#'   biplot.
#' @param grp A list of groups where each element is a vector containing the
#'   indices of group members (default: `NULL`).
#' @param grp.col A vector of the same length as `grp` containing colours for
#'   the corresponding groups (default: `NULL`, automatically set to `"black"`
#'   if `plot.groups= TRUE` and `grp.col= NULL`).
#' @param grp.sym Either a numeric vector the same length as `grp` containing
#'   only numbers 1-25 (corresponding to the `pch` argument of many graphical
#'   functions), or `"text"`. Indicates if samples will be plotted as symbols
#'   or text on the biplot (default: `"text"`).
#' @param grp.cex A numeric value indicating the relative size of the group
#'   symbols or text to be plotted (default: `1`).
#' @param loadings.grp A list of loading groups where each element is a vector
#'   containing the indices of group members (default: `NULL`).
#' @param loadings.col A vector the same length as `loadings.grp` containing
#'   colours for the corresponding loading groups (default: `NULL`,
#'   automatically set to `rgb(0,0,0,0.05)` if `plot.loadings= TRUE` and
#'   `loadings.col= NULL`).
#' @param loadings.sym Either a numeric the same length as `loadings.grp`
#'   containing only numbers 1-25 (corresponding to the `pch` argument of many
#'   graphical functions), or  `"text"`. Indicates if loadings will be plotted
#'   as symbols or text (default: `19`).
#' @param loadings.cex A numeric value indicating the relative size of the
#'   loading symbols or text to be plotted (default: `0.5`).
#' @param PC A numeric vector of length = 2 indicating which prinicpal
#'   components should be plotted (default: `c(1,2)`)
#' @param plot.legend A character string equal to either `"groups"` or
#'   `"loadings"`, or NULL, indicating whether to draw a legend for groups,
#'   loadings or neither (default: `NULL`).
#' @param leg.position A character string equal to one of the following:
#'   `"topleft"`, `"top"`, `"topright"`, `"right"`, `"bottomright"`, `"bottom"`,
#'   `"bottomleft"`, `"left"`, or `"center"`. Cannot be called in conjunction
#'   with `leg.xy` (default: `NULL`).
#' @param leg.xy A numeric vector of length = 2, specifying the x and y
#'   coordinates at which to draw the top-right corner of the legend box. Cannot
#'   be called in conjunctionn with `leg.position` (default: `NULL`)
#' @param leg.cex A numeric value indicating the relative size of the legend
#'   text (default: `0.55`).
#' @param leg.columns A numeric calue indicating the number of columns to split
#'   legend items into (default: `1`)
#' @param title A character string containing the desired title of the biplot
#'   (default: `""`).
#'
#' @return Returns a plot showing samples and/or loadings as specified by the
#'   user, which can be saved using one of the graphics devices (e.g. `png()`).
#'
#' @export
#'
#' @seealso [prcomp()],[ALDEx2]
#'
#' @importFrom car dataEllipse
#' @importFrom ggplot2 ggplot aes geom_density xlab ylab theme theme_bw
#' @importFrom ggplot2 scale_color_manual scale_x_continuous scale_y_continuous
#' @importFrom graphics points text abline par plot.new
#' @importFrom grDevices rgb
#' @importFrom magrittr  %>%
#' @importFrom dplyr arrange bind_rows
#' @importFrom grid plotViewport pushViewport popViewport
#' @importFrom gridBase baseViewports
#' @importFrom zCompositions cmultRepl
#' @importFrom stringr str_pad
#'
#' @examples
#' #generate test data: 100 features across 30 samples (split into two groups)
#' ex.data<-as.data.frame(matrix(data = NA,ncol = 30,nrow=100))
#' colnames(ex.data)<-c(paste0("Sample",rep("_A-",15),seq(1,15,1)),paste0("Sample",rep("_B-",15),seq(1,15,1)))
#' rownames(ex.data)<-paste0("Gene_",stringr::str_pad(string = seq(1,100,1),width = 4,side = "left",pad = "0"))
#'
#' for(i in 1:30){
#'   set.seed(123+i)
#'   high<-sample(c(10000:50000),30)
#'   med<-sample(c(1000:5000),30)
#'   low<-sample(c(0:10),40,replace = TRUE,prob = c(0.5,rep(0.05,10)))
#'   if(i %in% seq(1,15,1)){
#'     ex.data[,i]<-c(high,med,low)
#'   }else{ex.data[,i]<-c(low,high,med)
#'   }
#'   rm(list = c("i","high","med","low"))
#' }
#'
#' # generate test metadata
#' metadata<-as.data.frame(cbind(id=colnames(ex.data),cond=c(rep("Group1",15),rep("Group2",15))))
#' gene.list<-as.data.frame(cbind(id=rownames(ex.data),taxa=paste0("Species",rep(1:10,each=10))))
#'
#' # get list of indices for groups and loadings
#' ex.groups<-list()
#' for(i in levels(factor(metadata$cond))){
#'   ex.groups[[i]]<-which(metadata$cond==i)
#' }
#'
#' ex.loadings<-list()
#' for(j in levels(factor(gene.list$taxa))){
#'   ex.loadings[[j]]<-which(gene.list$taxa==j)
#' }
#'
#' # define group and loadings colours
#' ex.groups.col<-c("deepskyblue","orange2")
#' ex.load.col<-colorRampPalette(c("dodgerblue2","moccasin","maroon3"),space="rgb")(10)
#'
#' # impute zeros in test dataset
#' clr.input<-cmultRepl(t(ex.data), method = "CZM",label = 0)
#'
#' # log-ratio transformation of data
#' clr.data<-as.data.frame(apply(t(clr.input), 2, function(x) log(x) - mean(log(x))))
#'
#' # perform pca
#' pca.data<-prcomp(t(clr.data))
#'
#' # plot samples as coloured symbols by groups w/ ellipses & loadings as text
#' codaSeq.PCAplot(pca.data,plot.groups = TRUE,plot.loadings = TRUE,plot.ellipses = "groups",grp = ex.groups,
#'                 grp.col = ex.groups.col,grp.sym = 19,grp.cex = 0.8, loadings.sym = "text" ,PC = c(1,2),
#'                 plot.legend = "groups",leg.position = "bottomright",leg.columns = 1,title = "Groups")
#'
#' # plot samples as text with no colour, loadings as coloured symbols by group,
#' # ellipses for loadings groups, density plots for loadings
#' codaSeq.PCAplot(pca.data, plot.groups = FALSE, plot.loadings = TRUE, plot.ellipses = "loadings",
#'                 plot.density = "loadings",loadings.grp = ex.loadings,loadings.col = ex.load.col,loadings.sym = 19,
#'                 loadings.cex = 0.8,PC = c(1,2),plot.legend = "loadings",leg.position = "bottomright",
#'                 leg.columns = 5,title = "Loadings")
codaSeq.PCAplot <- function(pcx, plot.groups=FALSE, plot.loadings=TRUE, plot.ellipses=NULL, plot.density=NULL,
                            grp=NULL, grp.col=NULL, grp.sym="text", grp.cex= 1, loadings.grp=NULL, loadings.col=NULL,
                            loadings.sym=19, loadings.cex=0.5, PC=c(1,2), plot.legend=NULL,
                            leg.position=NULL, leg.xy=NULL, leg.cex=0.55, leg.columns=1, title=""){

  if((attributes(pcx)$class == "prcomp") == FALSE) stop("please use the prcomp function for the SVD")
  if(!is.logical(plot.groups) && !plot.groups %in% c(FALSE,TRUE)) stop("'plot.groups' must be a logical vector of length = 1")
  if(!is.logical(plot.loadings) && !plot.loadings %in% c(FALSE,TRUE)) stop("'plot.loadings' must be a logical vector of length = 1")
  if(!is.null(plot.ellipses) && length(plot.ellipses)!=1) stop("'plot.ellipses' must be a character vector of length = 1")
  if(!is.null(plot.ellipses) && !plot.ellipses %in% c("groups","loadings")) stop("'plot.ellipses' must be set to EITHER 'groups' or 'loadings'")
  if(!is.null(plot.legend) && length(plot.legend)!=1) stop("'plot.legend' must be a character vector of length = 1")
  if(!is.null(plot.legend) && !plot.legend  %in%  c("groups","loadings")) stop("'plot.legend' must be set to EITHER 'groups' or 'loadings'")
  if(!is.character(title) || !length(title)==1) stop("'title' must be a be a character vector of length = 1")
  if(!is.null(plot.legend) && !plot.legend %in% c("groups","loadings")) warning("'plot.legend' supplied but not set to either 'groups' or 'loadings'")
  if(plot.groups==FALSE && plot.legend=="groups") warning("'plot.legend' set to 'groups' but 'plot.groups' set to FALSE")
  if(plot.loadings==FALSE && plot.legend=="loadings") warning("'plot.legend' set to 'loadings' but 'plot.loadings' set to FALSE")

  if(!is.null(plot.density)){
    if(length(plot.density)!=1) stop("'plot.density' must be a character vector of length = 1")
    if(!plot.density %in% c("groups","loadings")) stop("'plot.density' must be set to either 'groups' or 'loadings'")
    if(plot.density=="groups" && is.null(grp)) stop("please supply a list of sample groupings")
    if(plot.density=="loadings" && is.null(loadings.grp)) stop("please supply a list of loading vectors")
    if(plot.groups==TRUE && is.list(grp)==FALSE) stop("please ensure 'grp' is a list")

    mvar <- sum(pcx$sdev^2)
    PC1 <- paste("PC",PC[1],": ", round(sum(pcx$sdev[PC[1]]^2)/mvar, 3))
    PC2 <- paste("PC",PC[2],": ", round(sum(pcx$sdev[PC[2]]^2)/mvar, 3))

    layout.matrix<-matrix(c(2,0,1,3),ncol = 2,nrow = 2,byrow = TRUE)
    layout(mat = layout.matrix,heights = c(1,2),widths = c(4,1))
    par(mar=c(4,4,1,1))

    plot(pcx$x[,PC[1]], pcx$x[,PC[2]], pch=19, xlab=PC1, ylab=PC2,
         col=rgb(0,0,0,0.01),  main= title,
         xlim=c(max(abs(pcx$x[,PC[1]])) * -1.2, max(abs(pcx$x[,PC[1]])) * 1.2),
         ylim=c(max(abs(pcx$x[,PC[2]])) * -1.2, max(abs(pcx$x[,PC[2]])) * 1.2)
    )
    abline(h=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))
    abline(v=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))
    limits<-par("usr")

    if(plot.loadings==TRUE){
      if(is.null(loadings.col)){
        if(plot.legend=="loadings") warning("no loadings colours specified: will not plot legend for loadings")
        if(!is.numeric(loadings.cex)) stop("'loadings.cex' must be numeric")
        if(length(loadings.cex)!=1) stop("'loadings.cex' must be a single numeric value")
        if(!is.numeric(loadings.sym) && !is.character(loadings.sym)) stop("'loadings.sym' must be either numeric or a character vector equal to 'text'")
        if(is.numeric(loadings.sym) && length(loadings.sym)!=1) stop("'loadings.sym' must be a single numeric value")
        if(is.numeric(loadings.sym) && !loadings.sym  %in% seq(0,25,1)) stop("'loadings.sym' must be a numeric value between 0 and 25")
        if(is.character(loadings.sym) && length(loadings.sym)!=1 && loadings.sym!="text") stop("to plot loading labels as text, 'loadings,sym' must be a character vector of length = 1, equal to 'text'")
        scale.factor.1 <- max(abs(pcx$x[,PC[1]]))/max(abs(pcx$rotation[,PC[1]]))
        scale.factor.2 <- max(abs(pcx$x[,PC[2]]))/max(abs(pcx$rotation[,PC[2]]))

        if(is.numeric(loadings.sym)){
          points(pcx$rotation[,PC[1]]*scale.factor.1,
                 pcx$rotation[,PC[2]]*scale.factor.2, pch=loadings.sym,
                 col=rgb(0,0,0,0.05), cex=loadings.cex)
        }else if(is.character(loadings.sym) && loadings.sym=="text"){
          text(pcx$rotation[,PC[1]]*scale.factor.1,
               pcx$rotation[,PC[2]]*scale.factor.2,
               labels=rownames(pcx$rotation),col=rgb(0,0,0,0.5),
               cex=loadings.cex)
        }
      }

      if(!is.null(loadings.col)){
        if(is.null(loadings.col)) stop("please provide a vector of loadings colours")
        if(length(loadings.col) >1 && is.list(loadings.grp)==FALSE) stop("multiple loadings colours supplied, but 'loadings.grp' is not a list")
        if(length(loadings.col) >1 && length(loadings.col) != length(loadings.grp)) stop("number of loadings groups and number of loading colours must match")
        if(!is.numeric(loadings.sym) && !is.character(loadings.sym)) stop("'loadings.sym' must be either numeric or a character vector equal to 'text'")
        if(is.numeric(loadings.sym) && length(loadings.sym)!=1) stop("'loadings.sym' must be a single numeric value")
        if(is.numeric(loadings.sym) && !loadings.sym  %in% seq(0,25,1)) stop("'loadings.sym' must be a numeric value between 0 and 25")
        if(is.character(loadings.sym) && length(loadings.sym)!=1 && loadings.sym!="text") stop("to plot loading labels as text, 'loadings,sym' must be a character vector of length = 1, equal to 'text'")

        if(is.numeric(loadings.sym)){
          loadings.sym <- rep(loadings.sym, length(loadings.col))
        }
        if(is.null(loadings.grp)){
          loadings.grp=list(c(1:nrow(pcx$rotation)))
        }
        if(is.numeric(loadings.grp) && length(loadings.sym) != length(loadings.grp)) stop("number of loadings groups and number of symbols must match")

        scale.factor.1 <- max(abs(pcx$x[,PC[1]]))/max(abs(pcx$rotation[,PC[1]]))
        scale.factor.2 <- max(abs(pcx$x[,PC[2]]))/max(abs(pcx$rotation[,PC[2]]))
        if(!is.numeric(loadings.cex)) stop("'loadings.cex' must be numeric")
        if(length(loadings.cex)!=1) stop("'loadings.cex' must be a single numeric value")
        if(is.numeric(loadings.sym)){
          for(j in 1:length(loadings.grp)){
            points(pcx$rotation[,PC[1]][loadings.grp[[j]]]*scale.factor.1,
                   pcx$rotation[,PC[2]][loadings.grp[[j]]]*scale.factor.2,
                   pch=loadings.sym[j],cex=loadings.cex, col=loadings.col[j])
          }
        }else if(loadings.sym=="text"){
          for(j in 1:length(loadings.grp)){
            text(pcx$rotation[,PC[1]][loadings.grp[[j]]]*scale.factor.1,
                 pcx$rotation[,PC[2]][loadings.grp[[j]]]*scale.factor.2,
                 labels=rownames(pcx$rotation)[loadings.grp[[j]]],
                 cex=loadings.cex,col=loadings.col[j])
          }
        }

        if(!is.null(plot.ellipses) && plot.ellipses=="loadings"){
          for(l in 1:length(loadings.grp)){
            dataEllipse(pcx$rotation[,PC[1]][loadings.grp[[l]]]*scale.factor.1,
                        pcx$rotation[,PC[2]][loadings.grp[[l]]]*scale.factor.2,
                        lwd = 1, levels=c(0.67), center.cex=FALSE,
                        plot.points=FALSE, add=TRUE, col=loadings.col[l],
                        fill = TRUE, fill.alpha = 0.25)
          }
        }

        if(!is.null(plot.legend) && plot.legend=="loadings"){
          if(is.null(leg.position) & is.null(leg.xy)) stop("please specify one of either 'leg.pos' or 'leg.xy'")
          if(!is.null(leg.position) & !is.null(leg.xy)) stop("please specify ONE of either 'leg.pos' or 'leg.xy'")
          if(is.null(leg.columns)) stop("please specify the number of columns in the legend")
          if(!is.numeric(leg.columns)) stop("'leg.cols' must be numeric")
          if(is.null(names(loadings.grp))) stop("please ensure group names in your list are not NULL using names()")
          if(!is.null(leg.position)){
            legend(leg.position,legend = names(loadings.grp),
                   col=loadings.col, cex=leg.cex,ncol= leg.columns,pch = if(is.character(loadings.sym)){
                     rep("-",length(loadings.grp))
                   }else loadings.sym)
          }else {
            if(length(leg.xy)!=2) stop("'leg.xy' must be a numeric vector of length = 2, specifying x and y coordinates, respectively")
            legend(leg.xy[1],leg.xy[2],legend = names(loadings.grp),
                   col=loadings.col, cex=leg.cex,ncol= leg.columns,pch = if(is.character(loadings.sym)){
                     rep("-",length(loadings.sym))
                   }else loadings.sym)
          }
        }
      }
    }

    if (plot.groups==FALSE) {
      if(!is.null(grp)) warning("'plot.groups'= FALSE but 'grp' is supplied: did you forget to set plot.groups to TRUE?")
      if(!is.null(grp.col)) warning("'plot.groups'= FALSE but 'grp.col' is supplied: did you forget to set plot.groups to TRUE?")
      if(!is.null(plot.ellipses) && plot.ellipses =="groups") warning("'plot.groups'= FALSE but 'plot.ellipses' is supplied: did you forget to set plot.groups to TRUE?")
      if(!is.character(grp.sym) && !is.numeric(grp.sym)) stop("'grp.sym' must either be numeric or a character vector equal to 'text'")
      if(is.character(grp.sym) && grp.sym!="text") stop("to plot samples as text, 'grp.sym' must be a character vector equal to 'text'")
      if(is.character(grp.sym) && length(grp.sym)!=1) stop("to plot samples as text, 'grp.sym' must be a character vector of length = 1, equal to 'text'")
      grp.col="black"
      grp=list(c(1:nrow(pcx$x)))
      names(grp)<-"Group_Name"
      if(is.character(grp.sym)){
        for(i in 1:length(grp)){
          text(pcx$x[grp[[i]],PC[1]],pcx$x[grp[[i]],PC[2]],
               labels=rownames(pcx$x)[grp[[i]]],cex=grp.cex,col=grp.col[i])
        }
      }else if(is.numeric(grp.sym)){
        if(!grp.sym  %in% seq(0,25,1)) stop("'grp.sym' must be a numeric value between 0 and 25")
        if(length(grp.sym) !=1) stop("'plot.grp' set to FALSE but more than one group symbol supplied")
        for(i in length(grp)){
          points(pcx$x[grp[[i]],PC[1]],pcx$x[grp[[i]],PC[2]],
                 pch=grp.sym[i],cex= grp.cex,col=grp.col[i])
        }
      }
    }else if(plot.groups==TRUE){
      if(is.null(grp)) stop("please supply a list of sample groupings")
      if(is.list(grp) == FALSE) stop("please supply a list of sample groupings")
      if(is.null(grp.col)) stop("please supply a vector of group colours")
      if(length(grp)!=length(grp.col)) stop("number of groups and number of group colours must match")
      if(!is.character(grp.sym) && !is.numeric(grp.sym)) stop("'grp.sym' must either be numeric or a character vector equal to 'text'")
      if(is.character(grp.sym) && grp.sym!="text") stop("to plot samples as text, 'grp.sym' must be a character vector equal to 'text'")
      if(is.character(grp.sym) && length(grp.sym)!=1) stop("to plot samples as text, 'grp.sym' must be a character vector of length = 1, equal to 'text'")
      if(is.character(grp.sym)){
        for(i in 1:length(grp)){
          text(pcx$x[grp[[i]],PC[1]],pcx$x[grp[[i]],PC[2]],
               labels=rownames(pcx$x)[grp[[i]]],cex=grp.cex,col=grp.col[i])
        }
      }else if(is.numeric(grp.sym)){
        if(!all(grp.sym  %in% seq(0,25,1))) stop("all values in 'grp.sym' must be numeric values between 0 and 25")
        if(length(grp.sym)==1){
          grp.sym<-rep(grp.sym,length(grp))
        }
        if(length(grp.sym) >1 && length(grp.sym) != length(grp)) stop("number of group symbols and number of groups must match (and must correspond)")
        for(i in 1:length(grp)){
          points(pcx$x[grp[[i]],PC[1]],pcx$x[grp[[i]],PC[2]],
                 pch=grp.sym[i],cex= grp.cex,col=grp.col[i])
        }
      }

      if(!is.null(plot.ellipses) && plot.ellipses=="groups"){
        for(l in 1:length(grp)){
          dataEllipse(pcx$x[,PC[1]][grp[[l]]],pcx$x[,PC[2]][grp[[l]]],
                      lwd = 1, levels=c(0.67), center.cex=FALSE,
                      plot.points=FALSE, add=TRUE,col=grp.col[l],
                      fill = TRUE, fill.alpha = 0.25)
        }
      }

      if(!is.null(plot.legend) && plot.legend=="groups"){
        if(is.null(leg.position) & is.null(leg.xy)) stop("please specify one of either 'leg.pos' or 'leg.xy'")
        if(!is.null(leg.position) & !is.null(leg.xy)) stop("please specify ONE of either 'leg.pos' or 'leg.xy'")
        if(is.null(leg.columns)) stop("please specify the number of columns in the legend")
        if(!is.numeric(leg.columns)) stop("'leg.cols' must be numeric")
        if(is.null(names(grp))) stop("please ensure group names in your list are not NULL using names()")
        if(!is.null(leg.position)){
          legend(leg.position,legend = names(grp),
                 col=grp.col, cex=leg.cex,ncol= leg.columns,pch=if(is.character(grp.sym)){
                   rep("-",length(grp))
                 }else grp.sym)
        }else {
          if(length(leg.xy)!=2) stop("'leg.xy' must be a numeric vector of length = 2, specifying x and y coordinates, respectively")
          legend(leg.xy[1],leg.xy[2],legend = names(grp),
                 col=grp.col, cex=leg.cex,ncol= leg.columns,pch=if(is.character(grp.sym)){
                   rep("-",length(grp))
                 }else grp.sym)
        }
      }
    }

    if(plot.density=="groups"){
      list.group<-list()
      for(k in 1:length(grp)){
        points.pc1<-pcx$x[,1][grp[[k]]]
        points.pc2<-pcx$x[,2][grp[[k]]]
        list.group[[names(grp)[k]]]<-as.data.frame(cbind(points.pc1,points.pc2))
      }

      cols.ordered<-as.data.frame(cbind(name=names(grp),col=grp.col)) %>%
        arrange(name)
      pca.vals.group<-bind_rows(list.group,.id = "group")

      dens.pc1<-pca.vals.group%>%
        ggplot(aes(x=points.pc1,col=group))+
        geom_density(show.legend = FALSE)+
        xlab("")+ylab("Density")+
        scale_color_manual(values=cols.ordered$col)+
        scale_x_continuous(limits = limits[1:2],expand = c(0,0))+
        scale_y_continuous(expand = c(0,0))+
        theme_bw()+theme(panel.grid = element_blank())+
        theme(axis.text.y = element_blank())+theme(axis.ticks.y = element_blank())

      dens.pc2<-pca.vals.group%>%
        ggplot(aes(y=points.pc2,col=group))+
        geom_density(show.legend = FALSE)+
        xlab("Density")+ylab("")+
        scale_color_manual(values=cols.ordered$col)+
        scale_y_continuous(limits = limits[3:4],expand = c(0,0))+
        scale_x_continuous(expand = c(0,0))+
        theme_bw()+theme(panel.grid = element_blank())+
        theme(axis.text.x = element_blank())+theme(axis.ticks.x = element_blank())

      vp1<-plotViewport(c(0,1.805,0,0.43))
      plot.new()
      vps<-baseViewports()
      pushViewport(vps$figure)
      print(dens.pc1,vp=vp1)
      popViewport()
      vp2<-plotViewport(c(1.8345,0,0.45,0))
      plot.new()
      vps<-baseViewports()
      pushViewport(vps$figure)
      print(dens.pc2,vp=vp2)
      par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

    } else if(plot.density=="loadings"){
      list.loadings<-list()
      for (k in 1:length(loadings.grp)) {
        points.pc1<-pcx$rotation[,1][loadings.grp[[k]]]*scale.factor.1
        points.pc2<-pcx$rotation[,2][loadings.grp[[k]]]*scale.factor.2
        list.loadings[[names(loadings.grp)[k]]]<-as.data.frame(cbind(points.pc1,points.pc2))
      }

      cols.ordered<-as.data.frame(cbind(name=names(loadings.grp),col=loadings.col)) %>%
        arrange(name)
      pca.vals.loadings<-bind_rows(list.loadings,.id = "group")


      dens.pc1<-pca.vals.loadings%>%
        ggplot(aes(x=points.pc1,col=group))+
        geom_density(show.legend = FALSE)+
        xlab("")+ylab("Density")+
        scale_color_manual(values=cols.ordered$col)+
        scale_x_continuous(limits = limits[1:2],expand = c(0,0))+
        scale_y_continuous(expand = c(0,0))+
        theme_bw()+theme(panel.grid = element_blank())+
        theme(axis.text.y = element_blank())+theme(axis.ticks.y = element_blank())

      dens.pc2<-pca.vals.loadings%>%
        ggplot(aes(y=points.pc2,col=group))+
        geom_density(show.legend = FALSE)+
        xlab("Density")+ylab("")+
        scale_color_manual(values=cols.ordered$col)+
        scale_y_continuous(limits = limits[3:4],expand = c(0,0))+
        scale_x_continuous(expand = c(0,0))+
        theme_bw()+theme(panel.grid = element_blank())+
        theme(axis.text.x = element_blank())+theme(axis.ticks.x = element_blank())

      vp1<-plotViewport(c(0,1.805,0,0.43))
      plot.new()
      vps<-baseViewports()
      pushViewport(vps$figure)
      print(dens.pc1,vp=vp1)
      popViewport()
      vp2<-plotViewport(c(1.8345,0,0.45,0))
      plot.new()
      vps<-baseViewports()
      pushViewport(vps$figure)
      print(dens.pc2,vp=vp2)
      par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
    }

  }else if(is.null(plot.density)){

    mvar <- sum(pcx$sdev^2)
    PC1 <- paste("PC",PC[1],": ", round(sum(pcx$sdev[PC[1]]^2)/mvar, 3))
    PC2 <- paste("PC",PC[2],": ", round(sum(pcx$sdev[PC[2]]^2)/mvar, 3))

    plot(pcx$x[,PC[1]], pcx$x[,PC[2]], pch=19, xlab=PC1, ylab=PC2,
         col=rgb(0,0,0,0.01),  main= title,
         xlim=c(max(abs(pcx$x[,PC[1]])) * -1.2, max(abs(pcx$x[,PC[1]])) * 1.2),
         ylim=c(max(abs(pcx$x[,PC[2]])) * -1.2, max(abs(pcx$x[,PC[2]])) * 1.2)
    )
    abline(h=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))
    abline(v=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))

    if(plot.loadings==TRUE){
      if(is.null(loadings.col)){
        if(!is.null(plot.legend) && plot.legend=="loadings") warning("no loadings colours specified: will not plot legend for loadings")
        scale.factor.1 <- max(abs(pcx$x[,PC[1]]))/max(abs(pcx$rotation[,PC[1]]))
        scale.factor.2 <- max(abs(pcx$x[,PC[2]]))/max(abs(pcx$rotation[,PC[2]]))
        if(!is.numeric(loadings.cex)) stop("'loadings.cex' must be numeric")
        if(length(loadings.cex)!=1) stop("'loadings.cex' must be a single numeric value")
        if(!is.numeric(loadings.sym) && !is.character(loadings.sym)) stop("'loadings.sym' must be either numeric or a character vector equal to 'text'")
        if(is.numeric(loadings.sym) && length(loadings.sym)!=1) stop("'loadings.sym' must be a single numeric value")
        if(is.numeric(loadings.sym) && !loadings.sym  %in% seq(0,25,1)) stop("'loadings.sym' must be a numeric value between 0 and 25")
        if(is.character(loadings.sym) && length(loadings.sym)!=1 && loadings.sym!="text") stop("to plot loading labels as text, 'loadings,sym' must be a character vector of length = 1, equal to 'text'")
        if(is.numeric(loadings.sym)){
          points(pcx$rotation[,PC[1]]*scale.factor.1, pcx$rotation[,PC[2]]*scale.factor.2, pch=loadings.sym,
                 col=rgb(0,0,0,0.05), cex=loadings.cex)
        }else if(is.character(loadings.sym) && loadings.sym=="text"){
          text(pcx$rotation[,PC[1]]*scale.factor.1,
               pcx$rotation[,PC[2]]*scale.factor.2,
               labels=rownames(pcx$rotation),col=rgb(0,0,0,0.5),
               cex=loadings.cex)
        }
      }

      if(!is.null(loadings.col)){
        if(is.null(loadings.col)) stop("please provide a vector of loadings colours")
        if(length(loadings.col) >1 && is.list(loadings.grp)==FALSE) stop("multiple loadings colours supplied, but 'loadings.grp' is not a list")
        if(length(loadings.col) >1 && length(loadings.col) != length(loadings.grp)) stop("number of loadings groups and number of colours must match")

        if(!is.numeric(loadings.sym) && !is.character(loadings.sym)) stop("'loadings.sym' must be either numeric or a character vector equal to 'text'")
        if(is.numeric(loadings.sym) && length(loadings.sym)!=1) stop("'loadings.sym' must be a single numeric value")
        if(is.numeric(loadings.sym) && !loadings.sym  %in% seq(0,25,1)) stop("'loadings.sym' must be a numeric value between 0 and 25")
        if(is.character(loadings.sym) && length(loadings.sym)!=1 && loadings.sym!="text") stop("to plot loading labels as text, 'loadings,sym' must be a character vector of length = 1, equal to 'text'")

        if(is.numeric(loadings.sym)){
          loadings.sym <- rep(loadings.sym, length(loadings.col))
        }
        if(is.null(loadings.grp)){
          loadings.grp=list(c(1:nrow(pcx$rotation)))
        }
        if(is.numeric(loadings.sym) && length(loadings.sym) != length(loadings.grp)) stop("number of loadings groups and number of symbols must match")

        scale.factor.1 <- max(abs(pcx$x[,PC[1]]))/max(abs(pcx$rotation[,PC[1]]))
        scale.factor.2 <- max(abs(pcx$x[,PC[2]]))/max(abs(pcx$rotation[,PC[2]]))
        if(!is.numeric(loadings.cex)) stop("'loadings.cex' must be numeric")
        if(length(loadings.cex)!=1) stop("'loadings.cex' must be a single numeric value")
        if(is.numeric(loadings.sym)){
          for(j in 1:length(loadings.grp)){
            points(pcx$rotation[,PC[1]][loadings.grp[[j]]]*scale.factor.1,
                   pcx$rotation[,PC[2]][loadings.grp[[j]]]*scale.factor.2,
                   pch=loadings.sym[j],cex=loadings.cex, col=loadings.col[j])
          }
        }else if(loadings.sym=="text"){
          for(j in 1:length(loadings.grp)){
            text(pcx$rotation[,PC[1]][loadings.grp[[j]]]*scale.factor.1,
                 pcx$rotation[,PC[2]][loadings.grp[[j]]]*scale.factor.2,
                 labels=rownames(pcx$rotation)[loadings.grp[[j]]],
                 cex=loadings.cex, col=loadings.col[j])
          }
        }

        if(!is.null(plot.ellipses) && plot.ellipses=="loadings"){
          for(l in 1:length(loadings.grp)){
            dataEllipse(pcx$rotation[,PC[1]][loadings.grp[[l]]]*scale.factor.1,
                        pcx$rotation[,PC[2]][loadings.grp[[l]]]*scale.factor.2,
                        lwd = 1, levels=c(0.67), center.cex=FALSE,
                        plot.points=FALSE, add=TRUE, col=loadings.col[l],
                        fill = TRUE, fill.alpha = 0.25)
          }
        }

        if(!is.null(plot.legend) && plot.legend=="loadings"){
          if(is.null(leg.position) & is.null(leg.xy)) stop("please specify one of either 'leg.pos' or 'leg.xy'")
          if(!is.null(leg.position) & !is.null(leg.xy)) stop("please specify ONE of either 'leg.pos' or 'leg.xy'")
          if(is.null(leg.columns)) stop("please specify the number of columns in the legend")
          if(!is.numeric(leg.columns)) stop("'leg.cols' must be numeric")
          if(is.null(names(loadings.grp))) stop("please ensure group names in your list are not NULL using names()")
          if(!is.null(leg.position)){
            legend(leg.position,legend = names(loadings.grp),
                   col=loadings.col, cex=leg.cex,ncol= leg.columns, pch=if(is.character(loadings.sym)){
                     rep("-",length(loadings.grp))
                   }else loadings.sym)
          }else {
            if(length(leg.xy)!=2) stop("'leg.xy' must be a numeric vector of length = 2, specifying x and y coordinates, respectively")
            legend(leg.xy[1],leg.xy[2],legend = names(loadings.grp),
                   col=loadings.col, cex=leg.cex,ncol= leg.columns,pch=if(is.character(loadings.sym)){
                     rep("-",length(loadings.grp))
                   }else loadings.sym)
          }
        }
      }
    }

    if (plot.groups==FALSE) {
      if(!is.null(grp)) warning("'plot.groups'= FALSE but 'grp' is supplied: did you forget to set plot.groups to TRUE?")
      if(!is.null(grp.col)) warning("'plot.groups'= FALSE but 'grp.col' is supplied: did you forget to set plot.groups to TRUE?")
      if(!is.null(plot.ellipses) && plot.ellipses=="groups") warning("'plot.groups'= FALSE but 'plot.ellipses' is supplied: did you forget to set plot.groups to TRUE?")
      if(!is.character(grp.sym) && !is.numeric(grp.sym)) stop("'grp.sym' must either be numeric or a character vector equal to 'text'")
      if(is.character(grp.sym) && grp.sym!="text") stop("to plot samples as text, 'grp.sym' must be a character vector equal to 'text'")
      if(is.character(grp.sym) && length(grp.sym)!=1) stop("to plot samples as text, 'grp.sym' must be a character vector of length = 1, equal to 'text'")
      grp.col="black"
      grp=list(c(1:nrow(pcx$x)))
      names(grp)<-"Group_NoName"
      if(is.character(grp.sym)){
        for(i in 1:length(grp)){
          text(pcx$x[grp[[i]],PC[1]],pcx$x[grp[[i]],PC[2]],
               labels=rownames(pcx$x)[grp[[i]]],cex=grp.cex,col=grp.col[i])
        }
      }else if(is.numeric(grp.sym)){
        if(!grp.sym  %in% seq(0,25,1)) stop("'grp.sym' must be a numeric value between 0 and 25")
        if(length(grp.sym) !=1) stop("'plot.grp' set to FALSE but more than one group symbol supplied")
        for(i in length(grp)){
          points(pcx$x[grp[[i]],PC[1]],pcx$x[grp[[i]],PC[2]],
                 pch=grp.sym[i],cex= grp.cex,col=grp.col[i])
        }
      }
    }else if(plot.groups==TRUE){
      if(is.null(grp)) stop("please supply a list of sample groupings")
      if(is.list(grp) == FALSE) stop("please supply a list of sample groupings")
      if(is.null(grp.col)) stop("please supply a vector of group colours")
      if(length(grp) != length(grp.col)) stop("number of groups and number of group colours must match")
      if(!is.character(grp.sym) && !is.numeric(grp.sym)) stop("'grp.sym' must either be numeric or a character vector equal to 'text'")
      if(is.character(grp.sym) && grp.sym!="text") stop("to plot samples as text, 'grp.sym' must be a character vector equal to 'text'")
      if(is.character(grp.sym) && length(grp.sym)!=1) stop("to plot samples as text, 'grp.sym' must be a character vector of length = 1, equal to 'text'")
      if(is.character(grp.sym)){
        for(i in 1:length(grp)){
          text(pcx$x[grp[[i]],PC[1]],pcx$x[grp[[i]],PC[2]],
               labels=rownames(pcx$x)[grp[[i]]],cex=grp.cex,col=grp.col[i])
        }
      }else if(is.numeric(grp.sym)){
        if(!all(grp.sym  %in% seq(0,25,1))) stop("all values in 'grp.sym' must be numeric values between 0 and 25")
        if(length(grp.sym)==1){
          grp.sym<-rep(grp.sym,length(grp))
        }
        if(length(grp.sym) >1 && length(grp.sym) != length(grp)) stop("number of group symbols and number of groups must match (and must correspond)")
        for(i in 1:length(grp)){
          points(pcx$x[grp[[i]],PC[1]],pcx$x[grp[[i]],PC[2]],
                 pch=grp.sym[i],cex= grp.cex,col=grp.col[i])
        }
      }

      if(!is.null(plot.ellipses) && plot.ellipses=="groups"){
        for(l in 1:length(grp)){
          dataEllipse(pcx$x[,PC[1]][grp[[l]]],pcx$x[,PC[2]][grp[[l]]],
                      lwd = 1, levels=c(0.67), center.cex=FALSE,
                      plot.points=FALSE, add=TRUE, col=grp.col[l], fill = TRUE,
                      fill.alpha = 0.25)
        }
      }

      if(!is.null(plot.legend) && plot.legend=="groups"){
        if(is.null(leg.position) & is.null(leg.xy)) stop("please specify one of either 'leg.pos' or 'leg.xy'")
        if(!is.null(leg.position) & !is.null(leg.xy)) stop("please specify ONE of either 'leg.pos' or 'leg.xy'")
        if(is.null(leg.columns)) stop("please specify the number of columns in the legend")
        if(!is.numeric(leg.columns)) stop("'leg.cols' must be numeric")
        if(is.null(names(grp))) stop("please ensure group names in your list are not NULL using names()")
        if(!is.null(leg.position)){
          legend(leg.position,legend = names(grp),
                 col=grp.col, cex=leg.cex,ncol= leg.columns,pch=if(is.character(grp.sym)){
                   rep("-",length(grp))
                 }else grp.sym)
        }else {
          if(length(leg.xy)!=2) stop("'leg.xy' must be a numeric vector of length = 2, specifying x and y coordinates, respectively")
          legend(leg.xy[1],leg.xy[2],legend = names(grp),
                 col=grp.col, cex=leg.cex,ncol= leg.columns,pch=if(is.character(grp.sym)){
                   rep("-",length(grp))
                 }else grp.sym)
        }
      }
    }

  }
}
