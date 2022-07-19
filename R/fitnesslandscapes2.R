if (FALSE) {
  ## Restart session and do (in shell):
  ## tar czf fitnesslandscapes2.tar.gz fitnesslandscapes2 & R CMD INSTALL fitnesslandscapes2.tar.gz
  ## or do:
  install.packages("~/fitnesslandscapes2", repos = NULL, type = "source") # cmd + shift + 0
  library(fitnesslandscapes2)
}

require("dplyr")
require("ggplot2")
require("vegan")
require("scatterplot3d")
require("fields")

# Easy dot product on data frame
dotprod2 <- function(DF,vec) {
  out <- c()
  for (rowid in 1:nrow(DF)) {
    row <- DF[rowid,]
    out[rowid] <- sum(row*vec)
  }
  return(out)
}

# exclude -> include
INCL_EXCL <- function(DF, EXCLUDE) {
  return(setdiff(colnames(DF),EXCLUDE))
}


# Performs PPR, PCA, and LDA on a data set
DRA <- DimReduction <- function(DF=df, EXCLUDE=FALSE, INCLUDE=FALSE, TYPE="PPR", VARIABLE="Fitness", IGNORE_PREVIOUS=TRUE, VERBOSE=TRUE) { #c("Identifier")
  # Manage incorrect TYPE
  if (!(TYPE %in% c("PPR","LDA","PCA","NMDS"))) {
    cat("Error: type incorrect\n")
    return(NULL)
  }
  
  # Figure out EXCLUSION versus INCLUSION
  if (EXCLUDE != FALSE && INCLUDE != FALSE) {
    if (VERBOSE==TRUE) cat("Warning: both EXCLUDE and INCLUDE triggered (perhaps implicitly), taking set difference\n")
    return(DimReduction(DF, EXCLUDE=FALSE, INCLUDE=setdiff(INCLUDE, EXCLUDE), TYPE=TYPE, VARIABLE=VARIABLE, IGNORE_PREVIOUS=IGNORE_PREVIOUS)) 
  } else if (EXCLUDE == FALSE && INCLUDE != FALSE) {
    exclude_vector <- setdiff(INCL_EXCL(DF,INCLUDE),VARIABLE)
    if (length(exclude_vector)==0) exclude_vector <- FALSE
    return(DimReduction(DF, EXCLUDE=exclude_vector, INCLUDE=FALSE, TYPE=TYPE, VARIABLE=VARIABLE, IGNORE_PREVIOUS=IGNORE_PREVIOUS))
  } else if(EXCLUDE == FALSE && INCLUDE == FALSE) {
    EXCLUDE <- c()
    exclusion_df <- data.frame()
  } else if (EXCLUDE != FALSE && INCLUDE == FALSE) {
    exclusion_df <- DF[,EXCLUDE]
    DF[,EXCLUDE] <- NULL
  }
  
  # Flag previously generated columns
  RDA <- c("PC1","PC2","PP1","PP2","LD1","LD2","NMDS1","NMDS2","RD1","RD2")
  if (!IGNORE_PREVIOUS) {
    for (rda in RDA) {
      if (rda %in% colnames(DF) & VERBOSE==TRUE) {
        cat("Warning: ",rda," detected in column names\n", sep="")
      }
    }
  } else if (IGNORE_PREVIOUS) {
    for (rda in RDA) {
      if (rda %in% colnames(DF) & VERBOSE==TRUE) {
        cat("Warning: ",rda," detected in column names, ignoring now\n", sep="")
        DF[,rda] <- NULL
      }
    }
  }
  
  # Perform dimensionality reduction
  if (TYPE == "PPR") {
    df.norm <- as.data.frame(scale(DF[,INCL_EXCL(DF, VARIABLE)]))
    explanatory <- as.matrix(df.norm)
    response <- as.matrix(DF[,VARIABLE])
    df.ppr <- stats::ppr(explanatory,response,nterms=2,maxterms=ncol(DF))
    pprdirections <- df.ppr$alpha
    colnames(pprdirections) <- c("PP1","PP2")
    PP1 <- dotprod2(df.norm, pprdirections[,"PP1"])
    PP2 <- dotprod2(df.norm, pprdirections[,"PP2"])
    PPR_columns <- cbind(PP1,PP2)
    colnames(PPR_columns) <- c("PP1","PP2")
    output <- list(
      PPR_columns,
      pprdirections,
      df.ppr
    )
    
    names(output) <- c("columns","weights","ppr")
    return(output)
  } else if (TYPE == "LDA") {
    FORMULA <- as.formula(paste(VARIABLE,paste(INCL_EXCL(DF,c(VARIABLE)),collapse=" + "), sep=" ~ "))
    df.lda <- MASS::lda(FORMULA, DF)
    LD1 <- dotprod2(DF[,INCL_EXCL(DF,VARIABLE)], df.lda$scaling[,"LD1"])
    if (ncol(df.lda$scaling) >= 2) {
      LD2 <- dotprod2(DF[,INCL_EXCL(DF,VARIABLE)], df.lda$scaling[,"LD2"])
      LDA_columns <- cbind(LD1,LD2)
      colnames(LDA_columns) <- c("LD1","LD2")
    } else {
      LDA_columns <- data.frame(LD1=LD1)
    }
    output <- list(
      LDA_columns,
      df.lda$scaling,
      df.lda
    )
    names(output) <- c("columns","weights","lda")
    return(output)
  } else if (TYPE == "PCA") {
    df.pca <- stats::prcomp(DF[,INCL_EXCL(DF,EXCLUDE)], center = TRUE, scale. = TRUE)
    PCA_columns <- df.pca$x[,c("PC1","PC2")]
    colnames(PCA_columns) <- c("PC1","PC2")
    output <- list(
      PCA_columns,
      df.pca$rotation,
      df.pca
    )
    names(output) <- c("columns","weights","pca")
    return(output)
  } else if (TYPE == "NMDS") {
    df.metric <- as.data.frame(scale(DF[,INCL_EXCL(DF,EXCLUDE)]))
    metric_matrix <- matrix(NA,nrow(df.metric),nrow(df.metric))
    rownames(metric_matrix) <- rownames(df.metric)
    colnames(metric_matrix) <- rownames(df.metric)
    for (row_i in 1:nrow(df.metric)) {
      for (row_j in row_i:nrow(df.metric)) {
        del_vector <- df.metric[row_i,]-df.metric[row_j,]
        metric_matrix[row_i,row_j] <- metric_matrix[row_j,row_i] <- sum(del_vector*del_vector)
      }
    }
    NMDS_result <- vegan::metaMDS(metric_matrix,k=2)
    output <- list(
      NMDS_result$points,
      NMDS_result
    )
    names(output) <- c("columns","nmds")
    return(output)
  }
}


bootstrap <- PPR_replicates <- function(DF=df, EXCLUDE=c("Identifier"), INCLUDE=FALSE, VARIABLE="Fitness", IGNORE_PREVIOUS=TRUE) {
  PPR <- DimReduction(DF=DF, EXCLUDE=EXCLUDE, INCLUDE=INCLUDE, TYPE="PPR", VARIABLE=VARIABLE, IGNORE_PREVIOUS=IGNORE_PREVIOUS)$weights
  PP1_replicates <- as.data.frame(matrix(NA, nrow=0, ncol=nrow(PPR)))
  colnames(PP1_replicates) <- rownames(PPR)
  PP2_replicates <- as.data.frame(matrix(NA, nrow=0, ncol=nrow(PPR)))
  colnames(PP2_replicates) <- rownames(PPR)
  for (i in 1:nrow(DF)) {
    the.rows <- setdiff(1:nrow(DF),i)
    new.weights <- DimReduction(DF=DF[the.rows,], EXCLUDE=EXCLUDE, INCLUDE=INCLUDE, TYPE="PPR", VARIABLE=VARIABLE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, VERBOSE=FALSE)$weights
    PP1_replicates[nrow(PP1_replicates)+1,] <- new.weights[,"PP1"]
    PP2_replicates[nrow(PP2_replicates)+1,] <- new.weights[,"PP2"]
  }
  
  PP1_weights <- PPR[,"PP1"]
  PP2_weights <- PPR[,"PP2"]
  
  PP1_dots <- c()
  PP2_dots <- c()
  for (i in 1:nrow(PP1_replicates)) {
    PP1_dots[i] <- sum(PP1_weights*PP1_replicates[i,])
    PP2_dots[i] <- sum(PP2_weights*PP2_replicates[i,])
  }
  
  return(list(
    PP1_replicates=PP1_dots,
    PP1_median=median(PP1_dots),
    PP1_ratio=sum(PP1_dots>0)/length(PP1_dots),
    PP2_replicates=PP2_dots,
    PP2_median=median(PP2_dots),
    PP2_ratio=sum(PP2_dots>0)/length(PP2_dots)
  ))
}



# Returns a fitness landscape given specific parameters
landscape <- TPS_landscape <- function(DF=df, x="PP1", y="PP2", output="contour", Theta=30, Phi=30, z="Fitness", x_name=x, y_name=y, z_name=z, Lambda="default", zlim=NULL) {
  par(mar=c(5,5,2,1)+.1)
  Var1 <- DF[,x]
  Var2 <- DF[,y]
  Fitness <- DF[,z]
  tp.m <- as.matrix(data.frame(v1=Var1, v2=Var2))

  if (Lambda == "default") {
    t <- fields::Tps(x=tp.m, Y=Fitness)
  } else if (Lambda == "special") {
    t <- fields::Tps(x=tp.m, Y=Fitness, lambda=0.02691373)
  } else {
    t <- fields::Tps(x=tp.m, Y=Fitness,lambda=Lambda)
  }

  if (output=="plotly") return(plotly::plot_ly(z=~fields::predictSurface(t)$z) %>% plotly::add_surface())
  else if (output=="contour" & is.null(zlim)) return(fields::surface(t, xlab=x_name, ylab=y_name, zlab=z_name))
  else if (output=="contour" & !is.null(zlim)) return(fields::surface(t, xlab=x_name, ylab=y_name, zlab=z_name,zlim=zlim))
  else if (output=="wireframe" & is.null(zlim)) return(fields::surface(t, xlab=x_name, ylab=y_name, zlab=z_name,type="p", theta=Theta, phi=Phi))
  else if (output=="wireframe" & !is.null(zlim)) return(fields::surface(t, xlab=x_name, ylab=y_name, zlab=z_name,type="p", theta=Theta, phi=Phi,zlim=zlim))
  else if (output=="matrix") {
    p <- fields::predictSurface(t)
    return(list(x=p$x, y=p$y, z=p$z))
  }
  else if (output=="model") return(t)
  else if (output=="scatter3") return(scatterplot3d::scatterplot3d(df[,c(x,y,z)], type="h"))
  else if (output=="scatter2") {
    return(ggplot(df, aes_string(x=x,y=y,color=z))+geom_point(size=3)+theme_classic()+scale_colour_gradient(
      low = "#173F5F",
      high = "#ED553B",
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    ))
  } else if (output=="fitness") {
    return(predict(t, DF[,c(x,y)]))
  }
  else print("Error: wrong `output` type")
  
  # https://stackoverflow.com/questions/18881546/creating-a-trellised-faceted-thin-plate-spline-response-surface
}

# Creates a 2D frequency-binning
binCounts <- function(x,y,increment_x,increment_y, pdf=TRUE) {
  x_seq <- seq(min(x),max(x),increment_x)
  y_seq <- seq(min(y),max(y),increment_y)
  counts <- expand.grid(x=x_seq, y=y_seq)
  rownames(counts) <- paste(counts$x,counts$y,sep=":")
  counts$z <- 0
  for (xv in x) {
    for (yv in y) {
      xi <- floor((xv-min(x))/increment_x)*increment_x+min(x)
      yi <- floor((yv-min(y))/increment_y)*increment_y+min(y)
      counts[paste(xi,yi,sep=":"),"z"] <- counts[paste(xi,yi,sep=":"),"z"]+1
    }
  }
  colnames(counts) <- c("x","y","counts")
  if (pdf) {
    counts$counts <- counts$counts/sum(counts$counts)
  }
  return(counts)
}

# Creates a TPS density surface based on a binning
distribution <- TPS_distribution <- function(DF=df,x="PP1",y="PP2",output="contour",Theta=30,Phi=30,pdf=FALSE, Lambda="default",x_name=x,y_name=y,z_name="Frequency") { #x_divisor=2,y_divisor=2,
  x_axis <- DF[,x]
  y_axis <- DF[,y]
  if (output!="distribution") {
    return(TPS_landscape(binCounts(x_axis, y_axis, increment_x=sd(x_axis)/3, increment_y=sd(y_axis)/3, pdf=pdf),
                         "x", "y", output, x_name=x_name, y_name=y_name, z="counts", z_name=z_name, Lambda=Lambda, Theta=Theta, Phi=Phi))
  } else if (output=="distribution") {
    # return(TPS_landscape(binCounts(x_axis, y_axis, increment_x=sd(x_axis)/3, increment_y=sd(y_axis)/3, pdf=pdf),
                         # "x", "y", output="fitness", x_name=x_name, y_name=y_name, z="counts", z_name=z_name, Lambda=Lambda, Theta=Theta, Phi=Phi))
    model <- TPS_landscape(binCounts(x_axis, y_axis, increment_x=sd(x_axis)/3, increment_y=sd(y_axis)/3, pdf=pdf),
                  "x", "y", "model", x_name=x_name, y_name=y_name, z="counts", z_name=z_name, Lambda=Lambda, Theta=Theta, Phi=Phi)
    return(predict(model, DF[,c(x,y)]))
  }
}

# TPS_distribution <- function(DF=df,x="PP1",y="PP2",x_name=x,y_name=y,z_name="Frequency") { #x_divisor=2,y_divisor=2,
#   return(ggplot(df.2018, aes_string(x=x, y=y))+
#            stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")+
#            theme_classic()+
#            labs(x=x_name,y=y_name,fill=z_name)+
#            scale_fill_distiller(palette="Spectral", direction=-1))
# }


fitnesslandscape <- function(DF=df,z="Fitness",smoothing="default",output="full") {
  DF <- DF[!is.na(colnames(DF)),]
  if (ncol(DF) == 3) {
    X <- setdiff(colnames(DF),z)[1]
    Y <- setdiff(colnames(DF),z)[2]
  } else {
    pprr <- PPR_replicates(DF=DF, EXCLUDE=FALSE, INCLUDE=FALSE, VARIABLE=z)
    cat("PPR median replicates: (PP1)",pprr$PP1_median, pprr$PP2_median,"\n")
    if (pprr$PP1_median > 0 & pprr$PP2_median) {
      X <- "PP1"
      Y <- "PP2"
      method <- "PPR"
    } else {
      X <- "PC1"
      Y <- "PC2"
      method <- "PCA"
    }
    DF[,c(X,Y)] <- DimReduction(DF=DF, EXCLUDE=FALSE, INCLUDE=FALSE, TYPE=method, VARIABLE=z)$columns
  }
  if (output!="simple") {
    ggplot(DF, aes(x=Fitness))+geom_histogram()+theme_classic()
    TPS_distribution(DF=DF,x=X,y=Y,output="wireframe",Theta=30,Phi=30,pdf=FALSE, Lambda=smoothing,x_name=X,y_name=Y,z_name="Frequency")
    TPS_distribution(DF=DF,x=X,y=Y,output="contour",Theta=30,Phi=30,pdf=FALSE, Lambda=smoothing,x_name=X,y_name=Y,z_name="Frequency")
    TPS_landscape(DF=DF, x=X, y=Y, output="wireframe", Theta=30, Phi=30, z=z, x_name=X, y_name=X, z_name=z, Lambda=smoothing)
  }
  TPS_landscape(DF=DF, x=X, y=Y, output="contour", Theta=30, Phi=30, z=z, x_name=X, y_name=Y, z_name=z, Lambda=smoothing)
}



# Returns a data-frame of an ellipse-segment of a TPS model given specific parameters
ellipse_at <- function(center, radius, increment, tps_model, pdf) {
  pdf(1,1) # To check PDF function works before investing runtime
  p_x <- center[1]
  p_y <- center[2]
  radius_x <- radius[1]
  radius_y <- radius[2]
  increment_x <- increment[1]
  increment_y <- increment[2]
  min_x <- p_x-radius_x
  max_x <- p_x+radius_y
  min_y <- p_y-radius_x
  max_y <- p_y+radius_y
  x_seq <- seq(min_x, max_x, increment_x)
  y_seq <- seq(min_y, max_y, increment_y)
  ellipse <- expand.grid(x=x_seq, y=y_seq)
  ellipse$z <- NA
  for (r in 1:nrow(ellipse)) {
    # cat("sum",(ellipse[r,"x"]-p_x)^2+(ellipse[r,"y"]-p_y)^2,"radius",radius^2,"\n")
    if ((ellipse[r,"x"]-p_x)^2/radius_x^2+(ellipse[r,"y"]-p_y)^2/radius_y^2 <= 1) {
      ellipse[r,"z"] <- 0
    } else {
      ellipse[r,"z"] <- NA
    }
  }
  ellipse <- dyplr::filter(ellipse, !is.na(z))
  ellipse$z <- stats::predict(tps_model, ellipse[,c("x","y")])
  ellipse$p <- NA
  for (r in 1:nrow(ellipse)) {
    ellipse[r,"p"] <- pdf(ellipse[r,"x"]-p_x,ellipse[r,"y"]-p_y)
  }
  return(ellipse)
}

square_at <- function(min_x, max_x, increment_x, min_y, max_y, increment_y, tps_model, VERBOSE=FALSE) {
  if (VERBOSE) print(min_x)
  if (VERBOSE) print(max_x)
  x_seq <- seq(min_x, max_x, increment_x)
  # print(x_seq)
  if (VERBOSE) print(min_y)
  if (VERBOSE) print(max_y)
  y_seq <- seq(min_y, max_y, increment_y)
  # print(y_seq)
  square <- expand.grid(x=x_seq, y=y_seq)
  square$z <- NA
  square$z <- predict(tps_model, square[,c("x","y")])
  return(square)
}

EPI <- function(DF=df,x="PP1",y="PP2",z="Fitness",Lambda="default",resolution=c(20,20)) {
  X <- seq(min(from=DF[,x]),to=max(DF[,x]),length.out=resolution[1])
  Y <- seq(min(from=DF[,y]),to=max(DF[,y]),length.out=resolution[2])
  phenotype_space <- expand.grid(X,Y)
  colnames(phenotype_space) <- c(x,y)
  fitnessModel <- TPS_landscape(DF=DF,x=x,y=y,z=z,Lambda=Lambda,output="model")
  phenotype_space$fitness <- predict(fitnessModel,phenotype_space)
  phenotype_space$fitness <- phenotype_space$fitness/max(phenotype_space$fitness,na.rm=T)
  distributionModel <- TPS_distribution(DF=DF,x=x,y=y,z=z,Lambda=Lambda,output="model")
  phenotype_space$distribution <- predict(distributionModel,phenotype_space[,c(x,y)])
  phenotype_space$distribution <- phenotype_space$distribution/max(phenotype_space$distribution,na.rm=T)
  return(1-(mean((abs(phenotype_space$fitness-phenotype_space$distribution))^(2),na.rm=T)))
  # return(1-mean((abs(phenotype_space$fitness-phenotype_space$distribution))^(1),na.rm=T))
}
# return(1-mean((phenotype_space$fitness-phenotype_space$distribution)^2,na.rm=T))
# return(1-mean(abs(phenotype_space$fitness-phenotype_space$distribution),na.rm=T))

