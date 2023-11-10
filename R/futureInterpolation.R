if (FALSE) {
  install.packages("~/Desktop/Research/Future\ Interpolation/futureInterpolate", repos = NULL, type = "source") # cmd + shift + 0
  library(futureInterpolate)
}

require("dplyr")
require("ggplot2")
require("vegan")
require("fields")
require("ks")
require("plotly")

# Turn exclusion into inclusion
INCL_EXCL <- function(DATA, EXCLUDE) {
  return(setdiff(colnames(DATA),EXCLUDE))
}

# Easy dot product on data frame
dot.prod <- function(DATA,vector) {
  out <- c()
  for (rowid in 1:nrow(DATA)) {
    row <- DATA[rowid,]
    out[rowid] <- sum(row*vector)
  }
  return(out)
}

# Prepare data frame
prepare.data <- function(DATA, INCLUDE=NULL, EXCLUDE=NULL, VARIABLE=NULL, IGNORE_PREVIOUS=TRUE, QUIETLY=FALSE) {
    # Limit data to relevant columns
    if (is.null(INCLUDE) & is.null(EXCLUDE)) { # Both EXCLUDE and INCLUDE are empty
        the.list <- colnames(DATA)
    } else if (is.null(INCLUDE) & !is.null(EXCLUDE)) { # Only EXCLUDE is relevant
        the.list <- INCL_EXCL(DATA,EXCLUDE)
    } else if (!is.null(INCLUDE) & is.null(EXCLUDE)) { # Only INCLUDE is relevant
        the.list <- INCLUDE
    } else if (!is.null(INCLUDE) & !is.null(EXCLUDE)) { # Both EXCLUDE and INCLUDE are relevant
        the.list <- setdiff(INCLUDE, EXCLUDE)
    }

    if (is.null(VARIABLE)) {
        data <- DATA[,the.list]
    } else if (!is.null(VARIABLE)) {
        data <- DATA[,unique(c(the.list,VARIABLE))]
    }

    # Remove and flag previous dimensionality reduction columns
    DRA.list <- c("PC1","PC2","PP1","PP2","LD1","LD2","NMDS1","NMDS2","DR1","DR2")
    if (!IGNORE_PREVIOUS & !QUIETLY) {
        # If IGNORE_PREVIOUS=False and QUIETLY=False, simply write to console
        for (dra.entry in DRA.list) {
            if (dra.entry %in% colnames(data)) {
                cat("Warning: ",dra.entry," detected in column names\n", sep="")
            }
        }
    } else if (IGNORE_PREVIOUS) {
        # If IGNORE_PREVIOUS=True, remove repeated columns
        for (dra.entry in DRA.list) {
            if (dra.entry %in% colnames(data)) {
                if (!QUIETLY) cat("Warning: ",dra.entry," detected in column names, ignoring now\n", sep="")
                data[,dra.entry] <- NULL
            }
        }
    }

    # Return output
    return(data)
}

# Perform PPR on data frame
perform.PPR <- function(DATA, INCLUDE=NULL, EXCLUDE=NULL, VARIABLE="Fitness", NTERMS=2, REPLICATE=FALSE, IGNORE_PREVIOUS=TRUE, QUIETLY=FALSE) {
    # Prepare data frame for PPR
    data <- prepare.data(DATA=DATA, INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, VARIABLE=VARIABLE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, QUIETLY=QUIETLY)

    if (!REPLICATE) {
        # Perform PPR
        df.norm <- as.data.frame(scale(data[,INCL_EXCL(data, VARIABLE)]))
        explanatory <- as.matrix(df.norm)
        response <- as.matrix(data[,VARIABLE])
        df.ppr <- stats::ppr(explanatory,response,nterms=NTERMS,maxterms=ncol(data))
        pprdirections <- df.ppr$alpha
        colnames(pprdirections) <- c("PP1","PP2")
        PP1 <- dot.prod(df.norm, pprdirections[,"PP1"])
        PP2 <- dot.prod(df.norm, pprdirections[,"PP2"])
        PPR_columns <- cbind(PP1,PP2)
        colnames(PPR_columns) <- c("PP1","PP2")

        # Produce output
        output <- list(
            PPR_columns,
            pprdirections,
            df.ppr
        )
        names(output) <- c("columns","weights","ppr")
        return(output)
    } else if (REPLICATE) {
        # Perform replicates
        PPR.full <- perform.PPR(REPLICATE=FALSE, DATA=data, INCLUDE=NULL, EXCLUDE=NULL, VARIABLE=VARIABLE, NTERMS=NTERMS, IGNORE_PREVIOUS=FALSE, QUIET=TRUE)
        PPR <- PPR.full$weights
        PP1_replicates <- as.data.frame(matrix(NA, nrow=0, ncol=nrow(PPR)))
        colnames(PP1_replicates) <- rownames(PPR)
        PP2_replicates <- as.data.frame(matrix(NA, nrow=0, ncol=nrow(PPR)))
        colnames(PP2_replicates) <- rownames(PPR)
        for (i in 1:nrow(data)) {
            the.rows <- setdiff(1:nrow(data),i)
            new.weights <- perform.PPR(REPLICATE=FALSE, DATA=data[the.rows,], INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, VARIABLE=VARIABLE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, NTERMS=NTERMS)$weights
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
            ppr=PPR.full,
            replicates=list(
                PP1_replicates=PP1_dots,
                PP1_median=median(PP1_dots),
                PP1_prop_positive=sum(PP1_dots>0)/length(PP1_dots),
                PP2_replicates=PP2_dots,
                PP2_median=median(PP2_dots),
                PP2_prop_positive=sum(PP2_dots>0)/length(PP2_dots)
            )
        ))
    }
}


# Perform PCA on data frame
perform.PCA <- function(DATA, INCLUDE=NULL, EXCLUDE=NULL, IGNORE_PREVIOUS=TRUE, QUIETLY=FALSE, VARIABLE="Ignore this parameter.") {
    # Prepare data frame for PCA
    data <- prepare.data(DATA=DATA, INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, VARIABLE=NULL, IGNORE_PREVIOUS=IGNORE_PREVIOUS, QUIETLY=QUIETLY)
    
    # Perform PCA
    df.pca <- stats::prcomp(data, center = TRUE, scale. = TRUE)
    PCA_columns <- df.pca$x[,c("PC1","PC2")]
    colnames(PCA_columns) <- c("PC1","PC2")

    # Produce output
    output <- list(
      PCA_columns,
      df.pca$rotation,
      df.pca
    )
    names(output) <- c("columns","weights","pca")
    return(output)
}


# Produce TPS plot
TPS <- function(DATA, AXES="PPR", OUTPUT=c("wireframe","contour"), POINTS=TRUE, Theta=30, Phi=30, Z="Fitness", XNAME=NULL, YNAME=NULL, ZNAME=Z, LAMBDA="default", ZLIM=NULL, COLOUR="rainbow") {
    par(mar=c(5,5,2,1)+.1)

    # Extract correct columns
    if (AXES=="PPR") {
        X <- "PP1"
        Y <- "PP2"
    } else if (AXES=="PCA") {
        X <- "PC1"
        Y <- "PC2"
    } else {
        XY <- strsplit(AXES,",")[[1]]
        X <- XY[1]
        Y <- XY[2]
    }
    Xaxis <- DATA[,X]
    Yaxis <- DATA[,Y]
    Zaxis <- DATA[,Z]
    tp.m <- as.matrix(data.frame(x=Xaxis, y=Yaxis))

    # Perform TPS
    if (LAMBDA=="default") {
        t <- fields::Tps(x=tp.m, Y=Zaxis)
    } else {
        t <- fields::Tps(x=tp.m, Y=Zaxis,lambda=LAMBDA)
    }

    # Set colour
    if (COLOUR=="rainbow") COLOUR <- fields::tim.colors(256)
    else if (COLOUR=="snow") COLOUR <- fields::snow.colors(n=256, alpha=1)
    else cat("Wrong input for variable `COLOUR`\n") 

    if (is.null(XNAME)) XNAME <- X
    if (is.null(YNAME)) YNAME <- Y
    if (is.null(ZNAME)) ZNAME <- Z

    # Produce vizualisation
    output.list <- list()
    if ("plotly" %in% OUTPUT) { 
        old.names <- names(output.list)
        if (is.null(old.names)) {
            new.names <- c("plotly")
        } else {
            new.names <- c(old.names, "plotly")
        }
        output.list <- append(output.list, function() { return(plotly::plot_ly(z=~fields::predictSurface(t)$z) %>% plotly::add_surface()) })
        names(output.list) <- new.names
    }
    if ("wireframe" %in% OUTPUT & is.null(ZLIM)) { 
        fields::surface(t, xlab=XNAME, ylab=YNAME, zlab=ZNAME,type="p", theta=Theta, phi=Phi, col=COLOUR)
    }
    if ("wireframe" %in% OUTPUT & !is.null(ZLIM)) { 
        fields::surface(t, xlab=XNAME, ylab=YNAME, zlab=ZNAME,type="p", theta=Theta, phi=Phi,zlim=ZLIM, col=COLOUR)
    }
    if ("contour" %in% OUTPUT & is.null(ZLIM)) {
        fields::surface(t, xlab=XNAME, ylab=YNAME, zlab=ZNAME, col=COLOUR)
        if (POINTS) {
            if (sum(unique(Zaxis)%in%c(0,1))==2 & length(unique(Zaxis))==2) {
                points(DATA[Zaxis==1,c(X,Y)], pch=20, col="black")
                points(DATA[Zaxis==0,c(X,Y)], pch=20, col="white")
            } else {
                points(Xaxis, Yaxis, pch=20, col=COLOUR[as.numeric(cut(Zaxis,breaks=256))])
            }
        }
    }
    if ("contour" %in% OUTPUT & !is.null(ZLIM)) { 
        fields::surface(t, xlab=XNAME, ylab=YNAME, zlab=ZNAME,zlim=ZLIM, col=COLOUR)
        if (POINTS) {
            if (sum(unique(Zaxis)%in%c(0,1))==2 & length(unique(Zaxis))==2) {
                points(DATA[Zaxis==1,c(X,Y)], pch=20, col="black")
                points(DATA[Zaxis==0,c(X,Y)], pch=20, col="white")
            } else {
                points(Xaxis, Yaxis, pch=20, col=COLOUR[as.numeric(cut(Zaxis,breaks=256))])
            }
        }
    }
    if ("matrix" %in% OUTPUT) {
        old.names <- names(output.list)
        if (is.null(old.names)) {
            new.names <- c("matrix")
        } else {
            new.names <- c(old.names, "matrix")
        }
        p <- fields::predictSurface(t)
        output.list <- append(output.list, function() { return(list(x=p$x, y=p$y, z=p$z)) })
        names(output.list) <- new.names
    }
    if ("model" %in% OUTPUT) { 
        old.names <- names(output.list)
        if (is.null(old.names)) {
            new.names <- c("model")
        } else {
            new.names <- c(old.names, "model")
        }
        output.list <- append(output.list, function() { return(t) })
        names(output.list) <- new.names
    }
    if ("scatter3" %in% OUTPUT) { 
        old.names <- names(output.list)
        if (is.null(old.names)) {
            new.names <- c("scatter3")
        } else {
            new.names <- c(old.names, "scatter3")
        }
        output.list <- append(output.list, function() { return(scatterplot3d::scatterplot3d(df[,c(X,Y,Z)], type="h")) })
        names(output.list) <- new.names
    }
    if ("scatter2" %in% OUTPUT) {
        old.names <- names(output.list)
        if (is.null(old.names)) {
            new.names <- c("scatter2")
        } else {
            new.names <- c(old.names, "scatter2")
        }
        output.list <- append(output.list, function() { return(ggplot(df, aes_string(x=X,y=Y,color=Z))+geom_point(size=3)+theme_classic()+scale_colour_gradient(
            low = "#173F5F",
            high = "#ED553B",
            space = "Lab",
            na.value = "grey50",
            guide = "colourbar",
            aesthetics = "colour"
        )) })
        names(output.list) <- new.names
    }
    if ("fitness" %in% OUTPUT) {
        old.names <- names(output.list)
        if (is.null(old.names)) {
            new.names <- c("fitness")
        } else {
            new.names <- c(old.names, "fitness")
        }
        output.list <- append(output.list, function() { return(predict(t, DATA[,c(X,Y)])) })
        # names(output.list) <- new.names
    }
    # if (length(output.list) <= 1) output.list <- unlist(output.list)
    return(output.list)
}

# Calculate new mean and standard deviation
eval_eqns <- function(DATA, alpha, lambda, INCLUDE=NULL, EXCLUDE=NULL, IGNORE_DRA=TRUE, QUIETLY=FALSE) {
    data <- prepare.data(DATA, INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, VARIABLE=NULL, IGNORE_PREVIOUS=IGNORE_DRA, QUIETLY=QUIETLY)
    
    n <- ncols(data)
    c <- matrix(NA, ncol=n, nrow=n)
    for (i in 1:n) {
        for (j in 1:i) {
            c[i,j] <- c[j,i] <- c(data[,i], data[,j])
        }
    }

    mu <- 1:n
    for (i in 1:n) mu[i] <- mean(data[,i])

    sig <- 1:n
    for (i in 1:n) sig[i] <- sqrt(c[i,i])

    mean_P = 0 # RESULT 1

    var_P <- 0 # RESULT 2
    for (i in 1:n) {
        for (j in 1:n) {
            var_P <- var_P + alpha[i]*alpha[j]*c[i,j]/(sig[i]*sig[j])
        }
    }

    mean_P_prime <- 0 # RESULT 3
    for (i in 1:n) {
        mean_P_prime <- mean_P_prime + alpha[i]*(lambda[i]-1)*mu[i]/sig[i]
    }

    var_P_prime <- 0 # RESULT 4
    for (i in 1:n) {
        for (j in 1:n) {
            var_P_prime <- var_P_prime + alpha[i]*alpha[j]*(1/(sig[i]*sig[j])*(lambda[i]*lambda[j]*(c[i,j]+mu[i]*mu[j])+(1-lambda[i]-lambda[j])*mu[i]*mu[j])-(lambda[i]-1)*(lambda[j]-1)*mu[i]*mu[j]/(sig[i]*sig[j]))
        }
    }

    return(c(
        MEAN_P=mean_P,
        SD_P=sqrt(var_P),
        MEAN_P_PRIME=mean_P_prime,
        SD_P_PRIME=sqrt(var_P_prime)
    ))
}

# Perform future interpolation
futureInterpolation <- function(DATA, ADJUSTMENTS, METHOD=NULL, VARIABLE="Fitness", INCLUDE=NULL, EXCLUDE=NULL, IGNORE_PREVIOUS=TRUE, PLOT=NULL, VIZ=NULL, QUIETLY=FALSE) {
    # Clarify METHOD=NULL
    if (is.null(METHOD)) {
        return(futureInterpolation(DATA, ADJUSTMENTS, METHOD=list(
            type="auto",
            nterms=2,
            beta.break=100, # Change this please!
            replicate.break=0 # This too!
        ), INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, PLOT=PLOT, VIZ=VIZ, QUIETLY=QUIETLY))
    } else if (!is.null(METHOD)) {
        if (typeof(METHOD)=="character") {
            if (METHOD=="PPR") {
                return(futureInterpolation(DATA, ADJUSTMENTS, METHOD=list(
                    type="PPR",
                    nterms=2
                ), INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, PLOT=PLOT, VIZ=VIZ, QUIETLY=QUIETLY))
            } else if (METHOD=="PCA") {
                return(futureInterpolation(DATA, ADJUSTMENTS, METHOD=list(
                    type="PCA"
                ), INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, PLOT=PLOT, VIZ=VIZ, QUIETLY=QUIETLY))
            }
        }
    }

    # Clarify PLOT=NULL
    if (!is.null(PLOT)) {
        return(futureInterpolation(DATA, ADJUSTMENTS, METHOD=METHOD, INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, PLOT=list(
            output=list(
                base=c("wireframe", "contour","model"),
                zoom=c("contour"),
                multiply="model",
                points=TRUE
            ),
            lambda="default"
        ), VIZ=VIZ, QUIETLY=QUIETLY))
    }

    # Clarify VIZ=NULL
    if (!is.null(VIZ)) {
        return(futureInterpolation(DATA, ADJUSTMENTS, METHOD=METHOD, INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, PLOT=PLOT, VIZ=list(
            Theta=30,
            Phi=30,
            xname=NULL,
            yname=NULL,
            zname=NULL,
            zlim=NULL,
            colour="rainbow"
        ), QUIETLY=QUIETLY))
    }

    # Casework on DRA type
    if (METHOD$type=="PPR") {
        data <- prepare.data(DATA, INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, VARIABLE=VARIABLE, QUIETLY=QUIETLY)
        dra <- perform.PPR(DATA, INCLUDE=NULL, EXCLUDE=NULL, VARIABLE=VARIABLE, NTERMS=METHOD$nterms, REPLICATE=FALSE, IGNORE_PREVIOUS=FALSE, QUIETLY=QUIETLY)
        data <- cbind(data, dra$columns)
        X <- "PP1"
        Y <- "PP2"
    } else if (METHOD$type=="PCA") {
        data <- prepare.data(DATA, INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, VARIABLE=NULL, QUIETLY=QUIETLY)
        dra <- perform.PCA(DATA, INCLUDE=NULL, EXCLUDE=NULL, IGNORE_PREVIOUS=FALSE, QUIETLY=QUIETLY)
        data <- cbind(data, dra$columns)
        X <- "PC1"
        Y <- "PC2"
    } else if (METHOD$type=="auto") {
        ppr.w.rep <- perform.PPR(DATA, INCLUDE=NULL, EXCLUDE=NULL, VARIABLE=VARIABLE, NTERMS=METHOD$nterms, REPLICATE=TRUE, IGNORE_PREVIOUS=NULL, QUIETLY=FALSE)
        beta <- ppr.w.rep$ppr$beta # check
        pp1.med <- ppr.w.rep$replicates$PP1_median
        pp2.med <- ppr.w.rep$replicates$PP2_median
        if (beta > beta.break & pp1.med > replicate.break & pp2.med > replicate.break) {
            if (QUIETLY) cat("Beta and replicate values high enough to justify using PPR\n")
            dra <- ppr.w.rep$ppr
            data <- cbind(data, dra$columns)
            X <- "PP1"
            Y <- "PP2"
        } else {
            if (QUIETLY) cat("Beta and replicate values not high enough to justify using PPR\n")
            return(futureInterpolation(DATA, ADJUSTMENTS, METHOD=list(
                type="PCA"
            ), VARIABLE=VARIABLE, INCLUDE=INCLUDE, EXCLUDE=EXCLUDE, IGNORE_PREVIOUS=IGNORE_PREVIOUS, PLOT=PLOT, VIZ=VIZ, QUIETLY=QUIETLY))
        }
    }

    # Casework on VIZ
    if (is.null(VIZ$xname)) VIZ$xname <- X
    if (is.null(VIZ$yname)) VIZ$yname <- Y
    if (is.null(VIZ$zname)) VIZ$zname <- Z
    if (is.null(VIZ$colour)) VIZ$colour <- "rainbow"

    # Casework on PLOT
    if (is.character(PLOT$output)) { PLOT$output <- list(base=PLOT$output, zoom=PLOT$output, multiply=PLOT$output, points=TRUE) }

    # Generate first figure
    first.figures <- TPS(DATA, TYPE=METHOD$type, OUTPUT=unique(c(PLOT$output$base,"model")), Theta=VIZ$Theta, Phi=VIZ$Phi, Z=VARIABLE, XNAME=VIZ$xname, YNAME=VIZ$yname, ZNAME=VIZ$zname, LAMBDA=PLOT$lambda, ZLIM=VIZ$zlim, COLOUR=VIZ$colour)
    

    TPS(DATA, TYPE=METHOD$type, OUTPUT=PLOT$output$base, Theta=VIZ$Theta, Phi=VIZ$Phi, Z=VARIABLE, XNAME=VIZ$xname, YNAME=VIZ$yname, ZNAME=VIZ$zname, LAMBDA=PLOT$lambda, ZLIM=VIZ$zlim, COLOUR=VIZ$colour)

    
    

}


