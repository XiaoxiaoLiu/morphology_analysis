
# ==================================================================================
# Function: Classify using FEM .
#
#   Input: class_csv_nm  
#           -> spreadsheet file name
#           -> from spreadsheed, determine the neuron_characteristics variable.
#
#
#
# ==================================================================================
#cl_fem <- function( class_csv_nm, terminal_node_label ) #, which_categorical_labels )
#{
  library( "FisherEM" )
  library( gtools )

  class_csv_nm    <- "/data/informatics/changkyul/Ephys/RUN07012015.Extended_0.25NAImputed_SpinyASpiny_Firing/XSel.csv"  # 702 neurons, ephys only


    df              <- read.csv( class_csv_nm, stringsAsFactors=FALSE )
    M               <- data.frame( df ) 
    dims            <- dim( M )

    start_col       <- 2
    stop_col        <- dims[2]  

    # Remove specimen id.
    i_spec_id       <- which( names(M) == 'specimen_id' )
    i_use           <- setdiff( start_col:stop_col, i_spec_id )

    M               <- M[,i_use]

 
  
  b_nans      <- is.nan( as.matrix( df ))

  he          <- names( M )
  b_na                    <- is.na( as.matrix( M ))
  useM                    <- data.matrix( M )

  dims                    <- dim( useM )
  any_nan_usem            <- c()
  any_infinite_usem       <- c()
  any_na_usem             <- c()
  for( j in 1:dims[2] )
  {
    useM[,j]                <- scale( useM[,j] )
    any_nan_usem[ j ]       <- any( is.nan( scale( useM[,j] ))) 
    any_infinite_usem[ j ]  <- any( is.infinite( scale( useM[,j] ))) 
   any_na_usem[ j ]        <- any( is.na( scale( useM[,j] )))
  }

  any( is.nan( useM  ))
  any( is.infinite( useM  ))
  any( is.na( useM  ))

  #out_xsel35      <-  fem(Y = XSel, K = c(14,16), model = "all", method = "svd", init = "kmeans", nstart = 75, disp = TRUE)
  out_xsel35_K40   <-  fem(Y = useM, K = c(12,20,40), model = "DB", method = "svd", init = "kmeans", nstart = 75, disp = TRUE)
  summary( out_xsel35_K40 )

  pdf("param.xsel.K=40",width=6,height=4,paper='special')
  plot( out_xsel35_K40 )
   dev.off()
  # par(mfrow=c(1,2))
  #           plot( df$X1, out$cls )

  save( "out_xsel35_K40", file="xsel35.K=40.RData" )

  min <- min( out$D )
  max <- max( out$D )
  dD  <- dim( out$D )
  image(1:dD[2], 1:dD[1], out$D, col=ColorRamp, xlab="", ylab="", axes=FALSE, zlim=c(min,max))


  # ========================================================================================
  # Write the data .csv file and the labels to file.
  # ========================================================================================
  load( "RData/xsel.K40" )
  label_nm  <- "/data/informatics/kylel/out/labels/702_neurons_40_nodes_parametric_modelDB/label.csv"
  MM        <- c()
  MM$id     <- df[ , 1 ]
  MM$labels <- out_xsel35_K40$cls
  write.csv( MM, file=label_nm, row.names=FALSE )

  # Not working yet.
  # out_sparse <- sfem( useM, model='DB', l1=0.05, init='kmeans', method="svd"  )

#}
  
