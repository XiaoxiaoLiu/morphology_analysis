
# ==================================================================================
# Function: Classify using sparse hierarchical clustering .
#
#   Input: class_csv_nm  
#           -> spreadsheet file name
#           -> from spreadsheed, determine the neuron_characteristics variable.
#
#
# e.g. csv   ../../data/Branch0.EphysBasic.Morph.Pia.Feature.07012015.ft01.6TN.csv  
#
# ==================================================================================
cl_fem <- function( class_csv_nm, terminal_node_label ) #, which_categorical_labels )
{
  library( "sparcl" )
  library( gtools )

  # For debugging
  #terminal_node_label               <- "X6TerminalNodes"
  if( 0 )
  {
    class_csv_nm                      <- "../../data/Branch0.EphysBasic.Morph.Pia.Feature.07012015.ft01.6TN.csv"
    df                                <- read.csv( class_csv_nm, stringsAsFactors=FALSE )
    M                                 <- data.frame( df ) 
    dims                              <- dim( M )

    start_col   <- 9
    stop_col    <- dims[2]-2
    M           <- M[,start_col:stop_col]
  }
  else
  {
    class_csv_nm                      <- "../../data/experiments/Branch0.EphysBasic.Morph.Pia.Feature.07012015.ft01.6TN.e1.csv"
    df                                <- read.csv( class_csv_nm, stringsAsFactors=FALSE )
    M                                 <- data.frame( df ) 
    dims                              <- dim( M )

    start_col   <- 9
    stop_col    <- dims[2]-3
    M           <- M[,start_col:stop_col]
  }
  
  
  b_nans      <- is.nan( as.matrix( df ))

  #
  he          <- names( M )
  M$cre_line  <- match( M$cre_line, unique( mixedsort( M$cre_line )))
  M$region    <- match( M$region  , unique( mixedsort( M$region   )))
  M$layer                 <- match( M$layer   , unique( mixedsort( M$layer    )))
  M$apical_dendrite_tag   <- match( M$apical_dendrite_tag , unique( mixedsort( M$apical_dendrite_tag    )))
  M$dendrite_type_tag     <- match( M$dendrite_type_tag   , unique( mixedsort( M$dendrite_type_tag    )))
  M$cell_shape            <- match( M$cell_shape   , unique( mixedsort( M$cell_shape    )))
  M$moment1               <- scale( M$moment1 )
  M$moment2               <- scale( M$moment2 )
  M$moment3               <- scale( M$moment3 )
  M$moment4               <- scale( M$moment4 )

  b_na                    <- is.na( as.matrix( M ))
  if( any( b_na ))
  {
    #impute.mean           <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
    for( i in which(sapply(    M, is.numeric ))) 
    {
         M[is.na(   M[, i]), i] <- mean(   M[, i],  na.rm = TRUE)
    }
  }
  #which( M[ b_na ] )
  #useM                    <- data.matrix( M[ ,start_col:stop_col ] )
  useM                    <- data.matrix( M ) #[ ,start_col:stop_col ] )
  b_na                    <- is.na( useM  )
  n_cols_usem             <- dim(useM)[2]
  if( any( b_na ))
  {
    # 
    # For whatever reason, this sometimes needs to be executed by "hand".
    # For i = 52, 53, 54
    #
    for( i in n_cols_usem ) #which( is.numeric( useM )))#, is.numeric ))) 
    {
      useM[is.na( useM[, i] ), i] <- mean( useM[, i],  na.rm = TRUE)
    }
  }

  dims                    <- dim( useM )
  any_nan_usem            <- c()
  any_infinite_usem       <- c()
  any_na_usem             <- c()
  for( j in 1:dims[2] )
  {
     if( j != 52 ) # 52 is all zeros.
     {
       useM[,j]                <- scale( useM[,j] )
     }
    #else
    #{
    #  stopifnot( sum( useM[,j] ) == 0 )
    #}
    any_nan_usem[ j ]       <- any( is.nan( scale( useM[,j] ))) 
    any_infinite_usem[ j ]  <- any( is.infinite( scale( useM[,j] ))) 
    any_na_usem[ j ]        <- any( is.na( scale( useM[,j] )))
  }

    #useM[ , j ] <- scale( useM[,j] )
  #i_exclude               <- c( he[ which( he == "moment1" ) ], he[ which( he == "moment2"

  any( is.nan( data.matrix( useM )))
  any( is.infinite( data.matrix( useM )))
  any( is.na( data.matrix( useM )))

  #K                       <- 4
  ##out                     <- fem( data.matrix( useM[ 1:40, 11:21] ), K, model='all', init='kmeans' )
  #out                     <- fem( useM[ 1:40, 11:21], K, model='all', init='kmeans' )
  #K                       <- 10
  #out                     <- fem( useM[  :  ,   :  ], K, model='all', init='kmeans' )
  #
  ##K      <- 10
  #out    <- fem( useM, K, model='DB', init='kmeans' )
  #out    <- fem( useM, K=3:8, model='all', init='kmeans', disp=true )
  #out <- fem(Y = useM, K = 3:8, model = "AkB", init = "kmeans", disp=true )
  #out <-  fem(Y = useM, K = 2:10, model = "all", method = "svd", init = "kmeans", kernel = "linear", nstart = 75, disp = TRUE)
  #out_sparse <- sfem( useM, model='DB', l1=0.05, init='kmeans', method="svd"  )

  # Sparse Hierarchical Clustering
  perm.out <- HierarchicalSparseCluster.permute(useM, wbounds=c(1.5,2:6), nperms=5)
  print(perm.out)
  plot(perm.out)


  # Perform sparse hierarchical clustering
  sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists, wbound=perm.out$bestw, method="complete")
  pdf("file.pdf",width=6,height=4,paper='special')
  par(mfrow=c(1,2))
  plot(sparsehc)
  dev.off()

  p <- dims[ 2 ] 

  y <- c( rep(1,p),rep(2,p) )

  pdf("sparcl_out2.pdf",width=6,height=4,paper='special')
  plot(sparsehc$hc )#, labels=rep("", length(y)))
  dev.off()
  print(sparsehc)
  # Plot using knowledge of class labels in order to compare true class
  # labels to clustering obtained
  par(mfrow=c(1,1))
  ColorDendrogram(sparsehc$hc,y=y,main="My Simulated Data",branchlength=.007)



}
  
  # Get the classification labels (i.e. the neuron types)
 # cell_types                        <- df[[ terminal_node_label ]]

  # Initialize neuron_characteristics
 # neuron_characteristics                  = c() 
 # neuron_characteristics$n_neurons        = length( cell_types )
 # neuron_characteristics$n_types          = length( levels( as.factor( cell_types )))
 # neuron_characteristics$types            = levels( as.factor( cell_types ))

  # Convert categorical variables from strings to integers.
 # df_w_factors                      <- df 
 # which_features                    <- names( df_w_factors )[ 11:length(names(df))-2 ]
 # features                          <- df_w_factors[ which_features ]
 # i_categorical                     <- !sapply( features, is.numeric )
 # features[i_categorical]           <- lapply( features[i_categorical], as.factor )
 # category_names                    <- names( i_categorical )
 # type_feature_centers              <- c()
 # for( j in 1:length( i_categorical ) )
 # {
 #   if( i_categorical[ j ] == TRUE )
 #   {
 #     this_cat_var                      <- category_names[ j ]
 #     this_factor                       <- features[ this_cat_var ]
 #     features[this_cat_var]            <- sapply( features[this_cat_var], as.integer )

 #     for( ntype in neuron_characteristics$types )
 #     {
 #       i_neurons_this_type             <- which( cell_types == ntype )
 #       this_feature_this_type          <- features[ i_neurons_this_type, j ]

        # Save the histogram of the levels for the factor.
        #type_feature_centers$ntype[ j ] <- mean( sapply( features[this_cat_var], as.integer ), na.rm=TRUE )
 #       type_feature_centers[this_cat_var, ntype]   = this_feature_this_type #hist( this_feature_this_type, breaks=max9  ) 
 #     }
 #   }
 #   else
 #   {
      # Compute the cluster statistics for this feature.
 #     for( ntype in neuron_characteristics$types )
 #     {
 #        i_neurons_this_type              <- which( cell_types == ntype )
 #        this_feature_this_type           <- features[ i_neurons_this_type, j ]
 #        type_feature_centers$ntype[ j ]  <- mean( this_feature_this_type, na.rm=TRUE )
 #     }
 #   }
 # }

 # neuron_characteristics$prevalence_type  = 'empirical'       # flags the source of the neuron characteristics.
 # neuron_characteristics$num_features     = length( which_labels ) 
 # neuron_characteristics
#
#} # function: 

