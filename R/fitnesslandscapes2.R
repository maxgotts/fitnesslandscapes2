require("dplyr")

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


# Performs PPR (and LDA) on a data set
DimReduction <- function(DF=df, EXCLUDE=c("Identifier"), INCLUDE=FALSE, TYPE="PPR", VARIABLE="Fitness", IGNORE_PREVIOUS=TRUE) {
  # Manage incorrect TYPE
  if (!(TYPE %in% c("PPR","LDA","PCA"))) {
    cat("Error: type incorrect\n")
    return(NULL)
  }
  
  # Figure out EXCLUSION versus INCLUSION
  if (EXCLUDE != FALSE && INCLUDE != FALSE) {
    cat("Warning: both EXCLUDE and INCLUDE triggered, taking set difference\n")
    return(DimReduction(DF, FALSE, setdiff(INCLUDE, EXCLUDE), TYPE, Variable)) 
  } else if (EXCLUDE == FALSE && INCLUDE != FALSE) {
    return(DimReduction(DF, INCL_EXCL(DF,INCLUDE), FALSE, TYPE, VARIABLE))
  } else if(EXCLUDE == FALSE && INCLUDE == FALSE) {
    EXCLUDE <- c()
    exclusion_df <- data.frame()
  } else if (EXCLUDE != FALSE && INCLUDE == FALSE) {
    exclusion_df <- DF[,EXCLUDE]
    DF[,EXCLUDE] <- NULL
  }
  
  # Flag previously generated columns
  RDA <- c("PC1","PC2","PP1","PP2","LD1","LD2")
  if (!IGNORE_PREVIOUS) {
    for (rda in RDA) {
      if (rda %in% colnames(DF)) {
        cat("Warning: ",rda," detected in column names\n", sep="")
      }
    }
  } else if (IGNORE_PREVIOUS) {
    for (rda in RDA) {
      if (rda %in% colnames(DF)) {
        cat("Warning: ",rda," detected in column names, ignoring now\n", sep="")
        DF[,rda] <- NULL
      }
    }
  }
  
  # Perform dimensionality reduction
  if (TYPE == "PPR") {
    df.norm <- as.data.frame(scale(DF[,INCL_EXCL(DF, "Fitness")]))
    explanatory <- as.matrix(df.norm)
    response <- as.matrix(DF[,VARIABLE])
    df.ppr <- stats::ppr(explanatory,response,nterms=2,maxterms=5)
    pprdirections <- df.ppr$alpha
    PP1 <- dotprod2(df.norm, pprdirections[,"term 1"])
    PP2 <- dotprod2(df.norm, pprdirections[,"term 2"])
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
    LD2 <- dotprod2(DF[,INCL_EXCL(DF,VARIABLE)], df.lda$scaling[,"LD2"])
    LDA_columns <- cbind(LD1,LD2)
    colnames(LDA_columns) <- c("LD1","LD2")
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
  }
}



# Returns a fitness landscape given specific parameters
TPS_landscape <- function(DF=df, x="PP1", y="PP2", output="contour", Theta=30, Phi=30, z="Fitness", x_name=x, y_name=y, z_name=z, Lambda="special") {
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
  if (output=="contour") return(fields::surface(t, xlab=x_name, ylab=y_name, zlab=z_name))
  if (output=="wireframe") return(fields::surface(t, xlab=x_name, ylab=y_name, zlab=z_name,type="p", theta=Theta, phi=Phi))
  if (output=="matrix") return(fields::predictSurface(t)$z)
  if (output=="model") return(t)
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
TPS_distribution <- function(DF=df,x="PP1",y="PP2",output="contour",Theta=30,Phi=30,pdf=FALSE, Lambda="special") { #x_divisor=2,y_divisor=2,
  x_axis <- DF[,x]
  y_axis <- DF[,y]
  return(TPS_landscape(binCounts(x_axis, y_axis, increment_x=sd(x_axis)/3, increment_y=sd(y_axis)/3, pdf=pdf),
                       "x", "y", output, x_name=x, y_name=y, z="counts", z_name="Frequency", Lambda=Lambda))
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
