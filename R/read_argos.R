pushTimes <- function (time, id, pushincr = 1)
{
  dups <- unlist(tapply(time, id, duplicated), use.names = FALSE)
  if (any(dups)) {
    time[dups] <- time[dups] + pushincr
    time <- Recall(time, id)
  }
  time
}

getit <- function(lab) as.character(x  %>% xml_find_all(sprintf("//%s", lab)) %>% xml_contents())

#' Read Argos XML records
#'
#'
#' @param x xml Argos files
#'
#' @return a dplyr tbl_df 
#' @export
#'
read_argosxml <- function(x) {
  al <- vector("list", length(x))
  for (ifile in seq_along(x)) {
    xd <- read_xml(x[ifile])
    tags <- "platformHexId"
    tags1 <- c("locationDate", "latitude", "longitude", "locationClass")
    tags2 <- c("longitude2", "latitude2", "errorRadius", "semiMajor", "semiMinor", "orientation", "hdop")
    tag <- as.character(xd  %>% xml_find_all(sprintf("//%s", tags)) %>% xml_contents())
    locs  <- setNames(lapply(tags1, getit) , tags1)
    diag <- setNames(lapply(tags2, getit), tags2)
    locs$locationClass <- ordered(locs$locationClass, c("Z", "B", "A", "0", "1", "2", "3"))
    locs$locationDate <- as.POSIXct(strptime(locs$locationDate, "%Y-%m-%dT%H:%M:%S"), tz = "UTC")
    locs$latitude <- as.numeric(locs$latitude)
    locs$longitude <- as.numeric(locs$longitude)
    locs$hdop <- as.integer(diag$hdop)
    diag$hdop <- NULL
    diag <- lapply(diag, as.numeric)
    d <- data_frame(platformHexId = rep(tag, length(locs$longitude)))
    this <- try(d %>% bind_cols(as_data_frame(locs), as_data_frame(diag)), silent = TRUE)
    if (inherits(this, "try-error")) {
      warning(sprintf("cannot parse this file: %s", x[ifile]))
      print(attr(this, 'condition')$message)
    } else {
      al[[ifile]] <- this
    }
  }
  do.call(bind_rows, al)
}


#' Read Argos PRV records.
#'
#' Read data from Argos PRV.
#'
#' See \code{\link{track_validity}} for details on the tests and fixes applied to records.
#' @param x Connection to PRV records
#' @param mungefix logical, should all track-validity errors be fixed.
#' @param dtFormat
#' @param tz
#' @param duplicateTime_eps
#' @param p4
#' @param verbose
#'
#' @return A dplyr tbl_df
#' @export
#' @importFrom dplyr data_frame mutate distinct %>%
read_argos <- function (x, mungefix = TRUE, dtFormat = "%Y-%m-%d %H:%M:%S",
                        tz = "UTC", duplicateTime_eps = 0.01, p4 =  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0",
                        verbose = FALSE) {
  dout <- vector("list", length(x))
  for (icon in seq_along(x)) {
    old.opt <- options(warn = -1)
    dlines <- strsplit(readLines(x[icon]), "\\s+", perl = TRUE)
    options(old.opt)
    loclines <- sapply(dlines, length) == 12
    if (any(loclines)) {
      dfm <- matrix(unlist(dlines[sapply(dlines, length) ==
                                    12]), ncol = 12, byrow = TRUE)
      if (dfm[1, 7] == "LC") {
        msg <- paste(" appears to be a diag file, skipping.",
                     "Use read_diag to obtain a dataframe. \n\n")
        cat("file ", icon, msg)
        next
      }
      df <- vector("list", 12)
      names(df) <- c("prognum", "ptt", "nlines", "nsensor",
                     "satname", "class", "date", "time", "latitude",
                     "longitude", "altitude", "transfreq")
      for (i in c(1:4, 9:12)) df[[i]] <- as.numeric(dfm[, i])
      for (i in 5:6) df[[i]] <- factor(dfm[, i])
      for (i in 7:8) df[[i]] <- dfm[, i]
      df <- as_data_frame(df) %>% mutate(utc = as.POSIXct(strptime(paste(df$date, df$time), dtFormat), tz))

      dout[[icon]] <- df
    }
    else {
      cat("Problem with file: ", x[icon], " skipping\n")
    }
  }
  if (all(sapply(dout, is.null)))
    stop("No data to return: check the files")
  dout <- do.call(bind_rows, dout)
  if (mungefix) {
    ##dout <- dout[order(dout$ptt, dout$gmt), ]
    ##dout <- dout[!duplicated(dout), ]
    dout <- dout %>% arrange(ptt, utc) %>% distinct()
    dt.by.id <- unlist(tapply(dout$utc, dout$ptt, function(x) c(-1,
                                                                diff(x))))
    dup.by.eps <- which(abs(dt.by.id) < duplicateTime_eps)
    if (length(dup.by.eps) >= 1) {
      if (verbose) {
        cat("Adjusting duplicate times\n.....\n")
        for (i in dup.by.eps) {
          ind <- i + (-2:1)
          print(cbind(dout[ind, c("ptt", "utc", "class")],
                      row.number = ind))
        }
      }
      dout$utc <- pushTimes(dout$utc, dout$ptt)
      if (verbose) {
        cat("\n  Adjusted records now: \n\n")
        for (i in dup.by.eps) {
          ind <- i + (-2:1)
          print(cbind(dout[ind, c("ptt", "utc", "class")],
                      row.number = ind))
        }
      }
    }
    if (any(dout$longitude > 180) & !grepl("+over", p4)) {
      msg <- paste("\nLongitudes contain values greater than 180,",
                   "assuming proj.4 +over\n\n")
      cat(msg)
      p4 <- sprintf("%s +over", p4)
    }
    dout$class <- ordered(dout$class, levels = c("Z", "B", "A", "0", "1", "2", "3"))
    #     coordinates(dout) <- c("longitude", "latitude")
    #     proj4string(dout) <- CRS(p4)
    #     test <- try(dout <- trip(dout, c("gmt", "ptt")))
    #     if (!is(test, "trip")) {
    #       cat("\n\n\n Data not validated: returning object of class ",
    #           class(dout), "\n")
    #       return(dout)
    #     }
    #     cat("\n\n\n Data fully validated: returning object of class ",
    #         class(dout), "\n")
    #     return(dout)
    #   }
    #   cat("\n\n\n Data not validated: returning object of class ",
    #       class(dout), "\n")
  }
    dout
  }


#' Read Argos DIAG streams
#'
#' Read data from Argos PRV.
#'
#' DIAG contains two sets of alternate coordinates,
#' @param x DIAG connect (or character vector of files)
#'
#' @return dplry tbl_df
#' @export
read_diag <- function (x)
{
  data <- vector("list", length(x))
  for (ifile  in seq_along(x)) {
    fl <- x[ifile]
    d <- readLines(fl)
    locs <- d[grep("LON1", d, ignore.case = TRUE)]
    tms <- d[grep("DATE", d, ignore.case = TRUE)]
    bad <- (grep("\\?", locs))
    if (length(bad) > 0) {
      if (length(locs[-bad]) == 0) {
        warning(paste("no valid locations in:", fl, "\n ...ignoring"))
        next
      }
      locs <- locs[-bad]
      tms <- tms[-(bad)]
    }
    dlines <- paste(locs, tms)
    dlines <- strsplit(dlines, "\\s+", perl = TRUE)
    reclen <- length(dlines[[1]])
    dfm <- matrix(unlist(dlines[sapply(dlines, length) ==
                                  reclen]), ncol = reclen, byrow = TRUE)
    lonlat <- dfm[, c(4, 7, 10, 13)]
    dic <- dfm[, c(14, 17, 18, 21, 24), drop = FALSE]
    id <- dic[, 1]
    gmt <- as.POSIXct(strptime(paste(dic[, 2], dic[, 3]),
                               "%d.%m.%y %H:%M:%S"), tz = "GMT")
    lq <- dic[, 4]
    iq <- dic[, 5]
    ll <- as.vector(lonlat)
    ll[grep("S", ll)] <- paste("-", ll[grep("S", ll)], sep = "")
    ll <- gsub("S", "", ll)
    ll[grep("N", ll)] <- paste("", ll[grep("N", ll)], sep = "")
    ll <- gsub("N", "", ll)
    ll[grep("E", ll)] <- paste("", ll[grep("E", ll)], sep = "")
    ll <- gsub("E", "", ll)
    ll[grep("W", ll)] <- paste("-", ll[grep("W", ll)], sep = "")
    ll <- gsub("W", "", ll)
    ll <- matrix(as.numeric(ll), ncol = 4)
    lon <- ll[, 2]
    lon2 <- ll[, 4]
    lq <- factor(lq, ordered = TRUE, levels = c("Z", "B",
                                                "A", "0", "1", "2", "3"))
    data[[ifile]] <- data_frame(lon1 = lon, lat1 = ll[,
                                                         1], lon2 = lon2, lat2 = ll[, 3], gmt = gmt, id = id,
                                   lq = lq, iq = iq)
  }
  do.call(bind_rows, data)
}
