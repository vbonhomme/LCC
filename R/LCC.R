# we load required packages
library(dplyr)

rm(list = ls())

.LCC_check <- function(path=path){
  # checks for presence of data files ######
  go <- TRUE
  lf <- list.files(pattern=".csv")
  if (!"FireA.csv" %in% lf) {
    cat(" * 'FireA.csv' is missing \n")
    go <- FALSE
  }
  if (!"FireC.csv" %in% lf) {
    cat(" * 'FireC.csv' is missing \n")
    go <- FALSE
  }
  if (!"FireAsc.csv" %in% lf) {
    cat(" * 'FireAsc.csv' is missing \n")
    go <- FALSE
  }
  if (!"FireCsc.csv" %in% lf) {
    cat(" * 'FireCsc.csv' is missing \n")
    go <- FALSE
  }
  # we use this flag to go ahead or stop
  if (!go) stop()

  # imports data files #########################################################
  FireA   <<- read.csv("FireA.csv", header=TRUE, sep=";")
  FireAsc <<- read.csv("FireAsc.csv", header=TRUE, sep=";")
  FireC   <<- read.csv("FireC.csv", header=TRUE, sep=";")
  FireCsc <<- read.csv("FireCsc.csv", header=TRUE, sep=";")
  Smpl    <<- read.csv("Smpl.csv",  header=TRUE, sep=";")
  Seedle  <<- read.csv("Seedle.csv", header=TRUE, sep=";")

  # checks column names ########################################################
  # same idea as above for the ok flag
  ok <- TRUE
  # column names
  if (!identical(colnames(Seedle), c("Depth", "SdlCounts", "SdlArea"))) {
    cat(" * Seedle.csv must contain three columns named 'Depth', 'SdlCounts' and 'SdlArea'\n")
    ok <- FALSE
  }

  if (!identical(colnames(Smpl), c("Depth", "Age_calBP", "SmplCount", "SmplArea"))) {
    cat(" * Smpl.csv must contain three columns named 'Depth', 'Age_calBP', 'SmplCount' and 'SmplArea'\n")
    ok <- FALSE
  }
  if (!ok) stop()
}


.LCC_match <- function(Fire_file){
  Fire_file <- read.csv(Fire_file, header=T,sep=";")
  # we prefer this column selection rather than named columns since the .csv
  # returned by CharAnalysis is not really clean
  CC.peaks  <- Fire_file[, 19] #peak final
  CC.cpeak  <- Fire_file[, 8] #char peak
  CC.thresh <- Fire_file[, 12]  #threshold final
  CC.age    <- Fire_file[, 2] #Age
  CC.res    <- diff(CC.age[1:2]) #resolution temporelle entre la profondeur 1 et 2


  # we subset Fire_file to create a data.frame that will be used throughout the process
  CC <- select(Fire_file,
               peaks  = 19, # peaks final
               cpeak  = 8,  # char peak
               thresh = 12, # threshold final
               age    = 2)  # Age

  # we add a new column, a logical : whether cpeak is above threshold
  # §§§§ ne faudrait-il pas un >= plutôt ? je pense que les conséquences ne sont que psychologiques
  CC <- mutate(CC, above = cpeak > thresh)

  # we deduce the temporal resolution
  .res <- diff(Fire_file[1:2, 2]) # temporal resolution


  peakC.ind = which(CC$peaks==1) # where we have peaks in CC
  peakC.n   = length(peakC.ind)  # and how many

  Smpl <- select(Smpl,
                 Area = 4,  # sum of areas per sample
                 Age_calBP = 2, # sample age
                 Count = 3, # charcoal number
                 Depth = 1) # their depth


  # where to store results
  overthresh.indC <- peakCsamples.ind <- peakCsamples.count <- peakCsamples.area <- peakCsamples.age <- list()
  overthresh.intervalC <- matrix(NA, nrow=length(peakC.ind), ncol=2)
  # cosmectics
  # §§§§ au choix, on peut la lever, mais ça fait une matrice plus jolie qu'on retourne à la fin, parmi d'autres
  # tous les "cosmetics" réfèrent à ce genre de make-up. si ça te plait pas, on dégage.
  nice_names <- paste0("peak", 1:peakC.n)
  dimnames(overthresh.intervalC) <- list(nice_names, c("ends", "starts"))

  # detects contiguous TRUEs in 'set', starting from the bottom-th position and going backward. returns the top-th position
  find_top <- function(set, bottom){
    top <- bottom
    while (set[top]) {
      top <- top - 1
      if (top == 0) break # if set starts with TRUEs
    }
    return(top + 1)
  }

  # we loop over every peak to find which samples contribute to every peak
  for (i in seq(along=peakC.ind)){
    bottom.i <- peakC.ind[i]
    top.i <- find_top(CC$above, bottom.i)
    overthresh.intervalC[i, ] <- c(CC$age[top.i], CC$age[bottom.i] + .res)

    # identify all samples contributing to the peak
    peakCsamples.ind[[i]] <- which(Smpl$Age_calBP >= overthresh.intervalC[i, 1] &
                                     Smpl$Age_calBP <= overthresh.intervalC[i, 2])

    # we add (if not already included) the sample which is just above
    dy <- overthresh.intervalC[i, 1] - Smpl$Age_calBP
    peakCsamples.ind[[i]] <- union(peakCsamples.ind[[i]], which.min(dy[dy>0]))

    # we finally extract original count, area and age data
    peakCsamples.count[[i]] <- Smpl$Count[ peakCsamples.ind[[i]] ]
    peakCsamples.area[[i]]  <- Smpl$Area[  peakCsamples.ind[[i]] ]
    peakCsamples.age[[i]]   <- Smpl$Age_calBP[ peakCsamples.ind[[i]] ]
  }

  # cosmetics
  names(peakCsamples.count) <- names(peakCsamples.area) <- names(peakCsamples.age) <- names(peakCsamples.ind) <- nice_names


  ############ Calculates large and small charcoals for every depth
  area_threshold <- 0.1 # argument
  Seedle <- filter(Seedle, SdlCounts != 0) %>%                    # removes empty
    mutate(Se.area.inf = ifelse(SdlArea < area_threshold, 1, 0),  # col with 1=small, 0=large
           Se.area.sup = ifelse(SdlArea >= area_threshold, 1, 0)) # col with 0=small, 1=small

  # we merge by $Depth and remove empty depths
  dates2 <- merge(Seedle, Smpl,by="Depth",all=TRUE) %>% na.omit()
  # we count, by Age_calBP
  Samp <- summarize(group_by(dates2, Age_calBP),
                    Suminf = sum(Se.area.inf),
                    Sumsup = sum(Se.area.sup))

  # §§§§§ pas sur d'avoir compris ce délire là (lignes 135:138 de ton script originel), même si les résultats sont les mêmes
  dates <- merge(data.frame(Age_calBP = unlist(peakCsamples.age)), Samp, by="Age_calBP", all.x=FALSE)

  # we attribute a unique number to peaks, according to the number of ages they encompass
  res2 <- data.frame(pic=rep(1:peakC.n, times = sapply(peakCsamples.age, length)), Age_calBP=unlist(peakCsamples.age))
  # we merge it with dates
  dates<-merge(dates, res2, by="Age_calBP",all= TRUE)
  # we assign a number of fires to every peak, for every depth
  dates<-merge(dates, Smpl, by="Age_calBP",all.x= TRUE)  #assigne un num?ro de feux ? chaque pic sur la totalit? des  profondeurs constituant le pic
  # we replace NAs with 0s
  dates<-replace(dates,is.na(dates),0)
  # ugly but turns factors into numbers
  dates$Suminf <- as.numeric(as.character(dates$Suminf))
  dates$Sumsup <- as.numeric(as.character(dates$Sumsup))

  ############ Calculates large and small charcoals for every depth
  FireCount <- group_by(dates, pic) %>%
    summarize(Sinf=sum(Suminf), Ssup=sum(Sumsup), Sarea=sum(Area), Age_calBP=max(Age_calBP)) %>%
    mutate(Count = Sinf + Ssup) %>%
    # §§§ cosmetics (capitale)
    rename(Pic=pic) %>%
    # we calculate the FRI §§§§ on le fait là, plus logique et aussi Dmax est substitué par Age_calBP direct. ok pour ça ?
    mutate(FRI = c(diff(Age_calBP), NA)) %>%
    # §§§ cosmectics, reorder columns
    select(Pic, Age_calBP, FRI, Count, Sinf, Ssup, Sarea) %>%
    unique()

  return(FireCount)
}

.LCC_wdw <- function(file, min.charnb = 10,
                     wdw.increment = 400, wdw.width = 1000,
                     FUN = median){
  FireCount <- file
  # we import the file. if not provided, choose it interactively
  #     if (missing(file)) file <- file.choose()
  #     Firecount <- write.csv(file=file, h=TRUE, sep=",")
  #
  # we retain peaks that have at least min_charnb
  # §§§ un >= me paraitrait plus logique ici (eg 'au moins 10', plutôt que 'plus de 10'. kétenpenses ?
  peaksC <- filter(FireCount, Count > min.charnb)

  #D?finir les fen?tres coulissantes

  top_i <- seq(min(Smpl$Age_calBP), max(Smpl$Age_calBP), by=wdw.increment)
  bottom_i <- top_i + wdw.width
  wdw <- data.frame(top_i, bottom_i)

  # we loop over peaksC to retain only peaks which large charcoals
  # are above the threshold function (median by default) within every peak

  # we create a list to store results
  pic_fenetreok <- list()

  for (i in 1:nrow(wdw)){
    # a subset encompassed within the i-th frame of the window
    peaksC_i <- filter(peaksC, Age_calBP > wdw$top_i[i] & Age_calBP <= wdw$bottom_i[i])
    # the thresold value
    cut_i   <- summarize(peaksC_i, FUN(Ssup)) %>% as.numeric()
    # a subset of the subset
    pic_fenetreok[[i]] <- filter(peaksC_i, Ssup > cut_i)
  }
  # we bind them all
  pic_fenetreok <- bind_rows(pic_fenetreok) %>% unique()
  return(pic_fenetreok)}



#' Performs Large Charcoal Count method
#'
#' @param path the path where to find the required files. \code{getwd()} by default
#' @param write.csv \code{logical} whether to write \code{.csv} files
#' @param min.charnb \code{numeric} minimum number of charcoals in a peak to retain it (10 by default)
#' @param wdw.increment \code{numeric} the increment for sliding the window (400 by default)
#' @param wdw.width \code{numeric} the width of the moving window (1000 by default)
#' @param FUN \code{function} returning a scalar, to apply to the charcoal retrieved within
#' a moving window. (\code{median} by default)
#'
#' @return a list with two components : \code{$Area} and \code{$Count}, each of them
#' containing three components: \code{$raw} (raw results from CharAnalysis),
#' \code{$sc} screened results with ARCO and \code{$LCC} screening with the LCC methods.
#'
#' \code{$Area} and \code{$Count} can be passed directly to \link{plot_screening}.
#' @details
#' This functions needs several input data, in 'path', by default the current working directory.
#' Have a look to those provided available theregi: 'http://www.github/vbonhomme/LCC/data.zip'
#'
#' I. Raw data:
#'
#' I.1. \strong{Seedle.csv}: a .csv table with charcoal-particle areas.
#' Must have three columns and as many rows as the number of observations. These columns
#' must be labelled as follows:
#' \itemize{
#' \item 'Depth': depth of samples
#' \item 'SdlCounts': number of charcoal particles in each sample
#' \item 'SdlArea': charcoal-particle areas
#' }
#'
#' I.2. \strong{Smpl.csv}: a .csv table with charcoal counts and charcoal areas.
#' Must have four columns and as many rows as the number of samples. These columns
#' must be labelled as follows:
#' \itemize{
#' \item 'Depth': depth of samples
#' \item 'Age_calBP': age estimate of samples
#' \item 'SmplCount': number of charcoal particles in each sample
#' \item 'SmplArea': cumulative charcoal area in each sample
#'}
#'
#' II. Fire episode reconstructions:
#' These are output files from the
#' \itemize{
#'
#' \item CharAnalysis program by Higuera et al. (2009), which
#' is freely available there \url{http://sites.google.com/site/charanalysis/}.
#' \item ARCO by Finsinger et al. (2014) which is freely available there \url{https://github.com/wfinsinger/ARCO}
#'}
#'
#' In file names below, A stands for Area; C for counts; sc for screening with ARCO.
#'
#' \itemize{
#' \item 'FireA.csv': CharAnalysis output table from analysis of CHARCOAL AREAS
#' \item 'FireAsc.csv': CharAnalysis output table from analysis of CHARCOAL AREAS and
#'  WITH ARCO, with those parameters \code{n.boot=10000}, \code{thresh.prob=0.95},
#' \code{win.width=1000} and \code{breakage=FALSE}. Also, \code{minCountP} must be < 1.0
#' (e.g. 0.05 as in Higuera et al. 2009)
#' \item 'FireC.csv': CharAnalysis output table from analysis of CHARCOAL COUNTS
#' \item 'FireCsc.csv': CharAnalysis output table from analysis of CHARCOAL COUNTS and
#' with minimum count test, i.e. with CharAnalysis parameter \code{minCountP < 1.0}, e.g. 0.05 as in Higuera et al. (2009).
#'  }
#'
#'  Note that relevant information from all CharAnalysis is done via column number because
#'  columns names are poorly formatted. Thus, files must have columns in the expected order, specifically:
#'  \itemize{
#'  \item Column 2:  'age Top_i', the age at the top of the interpolated sample
#'  \item Column 8:  'char Peak', the Cpeak component of CHAR
#'  \item Column 12: 'thresh FinalPos', the threshold used for peak identification
#'  \item Column 19: 'peaks Final', the boolean series representing identified peaks.
#'  }
#'
#' @references
#'
#' The method is detailed in a submitted (2015) paper by Remy et al.
#'
#' Papers cited in the 'Details' section:
#'\itemize{
#' \item Finsinger, W., R. Kelly, J. Fevre, and E.K. Magyari (2014). A guide to screening
#'     charcoal peaks in macrocharcoal-area records for fire episode reconstructions.
#'     \emph{The Holocene} 24(8):1002-108.
#'     It can be found there \url{http://hol.sagepub.com/content/24/8/1002}
#'
#' \item Higuera, P.E., L.B. Brubaker, P.M. Anderson, F.S. Hu, and T.A. Brown (2009)
#'     Vegetation mediated the impacts of postglacial climate change on fire regimes in
#'     the south-central Brooks Range, Alaska. Ecological Monographs 79(2) 201:219.
#'    It can be foudn there \url{http://www.esajournals.org/doi/abs/10.1890/07-2019.1}
#'    }
#' @examples
#' \dontrun{
#' # If you have all required files in your working directory.
#' # You can find a .zip with example data there: 'http://www.github/vbonhomme/LCC/data.zip'
#' results <- LCC()
#' plot_screening(results$Area, "Area")
#' plot_screening(results$Count, "Count")
#' # Also, you should now have FireA, FireAsc, FireC, FireCsc, Seedle and Smpl
#' available in your environment.
#' # if write.csv=TRUE, you should also have .csv files in your working directory.
#' }
#' @export
LCC <- function(path=getwd(), write.csv=FALSE,
                min.charnb = 10,
                wdw.increment = 400,
                wdw.width = 1000,
                FUN = median){
  .LCC_check(path)

  FireA_raw <- .LCC_match("FireA.csv")
  FireA_sc  <- .LCC_match("FireAsc.csv")
  FireA_LCC <- .LCC_wdw(FireA_sc,
                        min.charnb = min.charnb,
                        wdw.increment = wdw.increment,
                        wdw.width = wdw.width,
                        FUN = median)

  Area <- list(raw = FireA_raw,
               sc  = FireA_sc,
               LCC = FireA_LCC)

  FireC_raw <- .LCC_match("FireC.csv")
  FireC_sc  <- .LCC_match("FireCsc.csv")
  FireC_LCC <- .LCC_wdw(FireC_sc)

  Count <- list(raw = FireC_raw,
                sc  = FireC_sc,
                LCC = FireC_LCC)


  if (write.csv) {
    write.csv(FireA_raw, file="FireA_raw.csv", row.names = FALSE)
    write.csv(FireA_sc,  file="FireA_sc.csv", row.names = FALSE)
    write.csv(FireA_LCC, file="FireA_LCC.csv", row.names = FALSE)
    write.csv(FireC_raw, file="FireC_raw.csv", row.names = FALSE)
    write.csv(FireC_sc,  file="FireC_sc.csv", row.names = FALSE)
    write.csv(FireC_LCC, file="FireC_LCC.csv", row.names = FALSE)
  }

  return(list(Area=Area, Count=Count))
}


#' Plot the screening
#'
#' Typically after an LCC() step, plots the fire detected by the three methods
#' @param list_screening an object from LCC (either $Area or $Count)
#' @param title \code{character} a title for the plot
#' @examples
#' \dontrun{
#' # see LCC examples.
#' }
#' @export
plot_screening <- function(list_screening, title=""){
  if (!identical(c("raw", "sc", "LCC"), names(list_screening)))
    stop(" * list_screening must contains '$raw', '$sc' and '$LCC' components")
  Fire_raw <- list_screening$raw
  Fire_sc  <- list_screening$sc
  Fire_LCC <- list_screening$LCC
  FireALL  <- bind_rows(Fire_raw, Fire_sc, Fire_LCC)
  add_pts <- function(x, y, pch, ...){points(x$Age_calBP, rep(y, nrow(x)), pch=pch, ...)}

  par(mar=c(3, 1, 1, 1))
  ma <- max(FireALL$Age_calBP)
  plot(NA, xlim=c(max(FireALL$Age_calBP), 0), ylim=c(1, 6), ann=FALSE, frame=FALSE, axes=FALSE)
  add_pts(Fire_raw, 4, 19, col="grey50")
  add_pts(Fire_sc, 3, 19, col="grey20")
  add_pts(Fire_LCC, 2, 17, col="grey20")
  axis(1, at=c(seq(0, ma, by=500), ma))
  legend("topright", pch = c(19, 19, 17), col=c("grey50", "grey20", "grey20"),
         legend = c("CharAnalysis", "Screened", "LCC"), bty="n")
  title(title, line=-2)
}



######






