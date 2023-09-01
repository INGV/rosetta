# rosetta cycle function
#

exe.rosetta.cycle <- function(eqFile) {
  # first things first: check if there are enough MDPs in the data file ----
  eqraw <- read.csv(eqFile,
                    header = TRUE,
                    sep = "|")
  
  eqraw2 <-
    select(eqraw,
           ReferenceLongitude,
           ReferenceLatitude,
           ExpectedIntensity)
  
  eqNMO <- nrow(eqraw) # MDP number
  
  if (eqNMO >= eqIntNumb55Min) {
    # loading and parsing header to obtain eq metadata ---------------------
    inputID <-
      unlist(strsplit(read.csv(
        eqFile, header = T, sep = "|"
      )$X.EventID[1], split = "/"))[3]
    
    eqheadfile <-
      paste (downloadPath, inputID, "_dbmihead.txt", sep = "")
    
    download.file(paste(asmiServUrl, "?eventid=", inputID, "&format=text", sep = ""),
                  destfile = eqheadfile)
    
    
    eqheadraw <- read.csv(eqheadfile,
                          header = F,
                          sep = "|")
    
    eqName <- inputID
    eqID <- inputID
    eqDate <- eqheadraw$V2[2]
    eqEpiLon <- as.numeric(eqheadraw$V4[2])
    eqEpiLat <- as.numeric(eqheadraw$V3[2])
    eqDepDef <- as.numeric(eqheadraw$V5[2])
    
    
    eqMw <- as.numeric(eqheadraw$V11[2]) # Registered/estimanted Mw
    eqEpiArea <- eqheadraw$V13[2]
    
    # non-numeric to numeric transformation for MDPs -----------------------
    eqraw2$MCS <-
      with(eqraw2, ifelse(
        ExpectedIntensity == "11-12",
        11.5,
        ifelse(
          ExpectedIntensity == "10-11",
          10.5,
          ifelse(
            ExpectedIntensity == "9-10",
            9.5,
            ifelse(
              ExpectedIntensity == "8-9",
              8.5,
              ifelse(
                ExpectedIntensity == "7-8",
                7.5,
                ifelse(
                  ExpectedIntensity == "6-7",
                  6.5,
                  ifelse(
                    ExpectedIntensity == "5-6",
                    5.5,
                    ifelse(
                      ExpectedIntensity == "4-5",
                      4.5,
                      ifelse(
                        ExpectedIntensity == "3-4",
                        3.5,
                        ifelse(
                          ExpectedIntensity == "2-3",
                          2.5,
                          ifelse(
                            ExpectedIntensity == "F",
                            3.9,
                            ifelse(
                              ExpectedIntensity == "SF",
                              2.9,
                              ifelse(
                                ExpectedIntensity == "NF",
                                1,
                                ifelse(
                                  ExpectedIntensity == "HF",
                                  5.1,
                                  ifelse(
                                    ExpectedIntensity == "SD",
                                    5.6,
                                    ifelse(
                                      ExpectedIntensity == "D",
                                      6.4,
                                      ifelse(ExpectedIntensity == "HD", 8.6, ExpectedIntensity)
                                    )
                                  )
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      ))
    
    eqraw3 <-
      select(eqraw2, ReferenceLongitude, ReferenceLatitude, MCS)
    
    colnames(eqraw3) <- c("Lon", "Lat", "MCS")
    
    # MCS are still recognized as "char"
    eqraw3 <- transform(eqraw3, MCS = as.numeric(MCS))
    
    # calculating distances of municipalities from epicenter in km
    distances <-
      sapply(1:nrow(eqraw2), function(i)
        round(distm(eqraw2[i, 1:2], c(
          eqEpiLon , eqEpiLat
        )) / 1000, digits = 1))
    
    eqraw3$Dist_epi <- distances
    
    # set a variable to enable/disable subsequent calculations
    doNotCalc <- FALSE
    
    ## Magnitude filter function ----------------------------------------
    
    # Establishing base filter parameters to select earthquakes
    eq <- select(eqraw3, MCS, Dist_epi)
    eqUnder55 <-
      nrow(subset(eq, eq$Dist_epi <= 55)) # number of MDP between epicenter
    # and 55 km
    if (eqUnder55 > 0) {
      ## rosetta algorithm calculations ---------------------------------
      
      # renaming header columns to something more practical
      names(eq)[1] <- "Int"
      names(eq)[2] <- "Dist"
      
      # adding nMDP vs mobile windows histogram plots
      eqMDPs <- paste(eqID, "MDPs", sep = "")
      
      # cleaning data file from fake (0) values
      eq <- eq[!(eq$Int == 0),]
      
      # creating subset of intensities within ranges of 10 km
      # moving away from epicenter
      v0_10 <- subset(eq[1], eq[2] >= 0 & eq[2] < 10)
      v5_15 <- subset(eq[1], eq[2] >= 5 & eq[2] < 15)
      v10_20 <- subset(eq[1], eq[2] >= 10 & eq[2] < 20)
      v15_25 <- subset(eq[1], eq[2] >= 15 & eq[2] < 25)
      v20_30 <- subset(eq[1], eq[2] >= 20 & eq[2] < 30)
      v25_35 <- subset(eq[1], eq[2] >= 25 & eq[2] < 35)
      v30_40 <- subset(eq[1], eq[2] >= 30 & eq[2] < 40)
      v35_45 <- subset(eq[1], eq[2] >= 35 & eq[2] < 45)
      v40_50 <- subset(eq[1], eq[2] >= 40 & eq[2] < 50)
      v45_55 <- subset(eq[1], eq[2] >= 45 & eq[2] < 55)
      
      # create dataframe with mean distance and MDP number to plot in histogram
      eq_MDPs_number <- cbind.data.frame(
        AvgDist = seq(5, 50, by = 5),
        nMDP = c(
          nrow(v0_10),
          nrow(v5_15),
          nrow(v10_20),
          nrow(v15_25),
          nrow(v20_30),
          nrow(v25_35),
          nrow(v30_40),
          nrow(v35_45),
          nrow(v40_50),
          nrow(v45_55)
        )
      )
      
      # calculating averages for distance intervals
      eq_avgs <- cbind.data.frame(
        Dist = seq(5, 50, by = 5),
        Int = c(
          mean(v0_10$Int),
          mean(v5_15$Int),
          mean(v10_20$Int),
          mean(v15_25$Int),
          mean(v20_30$Int),
          mean(v25_35$Int),
          mean(v30_40$Int),
          mean(v35_45$Int),
          mean(v40_50$Int),
          mean(v45_55$Int)
        )
      )
      
      # Consider only earthquakes with minimum mobile avgs and minimum avg
      # between 5 and 15 km from epicenter as established
      mobAvgNumb <- nrow(subset(eq_avgs, !is.nan(eq_avgs$Int)))
      if (!is.na(mean(v0_10$Int)) &&
          mean(v0_10$Int) > meanv0_10IntMin &&
          mobAvgNumb >= mobAvgNumbMin) {
        eq_IntNumbTot <- length(eq[[1]]) # MDPs total
        eq_IntNumb55 <-
          length(subset(eq[1], eq[2] < 55)[[1]]) # MDPs within 55 km
        
        # calculating earthquake depth using Rosetta algorithm
        cAng <- lm(eq_avgs$Int ~ eq_avgs$Dist)$coefficient[2]
        slope <- abs(as.numeric(cAng))
        eq_Dnext <- exp((slope - 0.087) / (-0.018))
        
        # Calculating intercept
        slopeInt <-
          summary(lm(eq_avgs$Int ~ eq_avgs$Dist))$coefficient[1]
        
        # calculating intercept Mw
        MwI <- 0.18 * log(eq_Dnext) + 0.56 * slopeInt + 1.44
        
        # set to NA parameters that might not be created if the following
        # "if" statement will not be performed (MwI < 6.75)
        slopeIntCorr <- NA
        MwICorr <- NA
        Rwc <- NA
        
        ## W & C recalculation for MwI >= 6.75 ---------------------------
        # If MwI is >= of 6.75 re-calculate everything taking into account
        # fault radius as for Wells and Coppersmith (1994)
        if (MwI >= 6.75) {
          # source R script for fault parameters calculation through W&C1994
          source("wcCalcFun.R")
          wcCalc(MwI)
          load("faultMetaWC.RData")
          Rwc <- round(faultRadius, digits = 0)
          
          ### recalculate earthquake parameters ---------------------------
          startRange <- Rwc
          endRange <- Rwc + 10
          
          # create vector for number of MDPs to be used to create histograms of
          # distance vs number of MDPs
          number_of_mdps <- c()
          
          # first two ranges (fixed values)
          vRwc <- subset(eq, eq$Dist >= 0 & eq$Dist <= Rwc)
          vRwcplus <-
            subset(eq, eq$Dist > startRange & eq$Dist <= endRange)
          
          # calculates the average intensity between epicenter and Rwc and add the first
          # row to the dataframe
          eq_avgs <-
            data.frame(Rwc / 2, round(mean(vRwc$Int), digits = 1))
          
          # add number of MDPs to vector
          number_of_mdps <- append(number_of_mdps, nrow(vRwc))
          
          # gives a name to columns
          names(eq_avgs) <- c("Dist", "Int")
          
          # calculates the average intensity between Rwc and Rwc+10km (2nd window) and
          # add it to the dataframe
          eq_avgs[nrow(eq_avgs) + 1, ] = c((Rwc + 5),
                                           round(mean(vRwcplus$Int),
                                                 digits = 1))
          number_of_mdps <- append(number_of_mdps, nrow(vRwcplus))
          
          # move forward
          startRange <- startRange + 5
          endRange <- startRange + 10
          
          # subsetting remnant eq intensities following above directions; calculating
          # average intensity for every subset and storing results in the dataframe
          while (endRange <= 55) {
            nam <- paste("v", startRange, endRange, sep = "")
            assign(nam, subset(eq, eq[2] > startRange &
                                 eq[2] <= endRange))
            
            #transform nam character in variable name
            namX <- eval(parse(text = nam))
            
            # add numbero of MDPs to vector to plot histograms
            number_of_mdps <- append(nrow(namX), number_of_mdps)
            
            # insert mean in eq_avgs dataframe
            eq_avgs[nrow(eq_avgs) + 1, ] = c(startRange + 5,
                                             round(mean(namX$Int), digits = 1))
            
            # move forward
            startRange <- startRange + 5
            endRange <- startRange + 10
          }
          
          eq_MDPs_number <-
            cbind.data.frame(eq_avgs$Dist, number_of_mdps)
          
          names(eq_MDPs_number)[1] <- "AvgDist"
          names(eq_MDPs_number)[2] <- "nMDP"
          
          # plot of histograms
          if (histYES) {
            hist <- ggplot(eq_MDPs_number,
                           aes(x = AvgDist, y = nMDP)) +
              geom_bar(stat = "identity",
                       color = "blue",
                       fill = "transparent") +
              scale_x_continuous(
                name = "Epicentral distance (km)",
                limits = c(0, 55),
                breaks = round(eq_MDPs_number$AvgDist, digits = 0)
              ) +
              scale_y_continuous(name = "Number of MDPs") +
              labs(title = paste(
                eqID,
                " - ",
                eqEpiArea,
                "\n# MDP vs epicentral distance",
                sep = ""
              )) +
              geom_text(aes(y = nMDP - nMDP * 10 / 100, label = nMDP))
            
            ggsave(
              filename = paste(outdir, "/", eqID, "WC_hist.png", sep = ""),
              hist,
              device = "png"
            )
          }
          
          # Consider only earthquakes with minimum mobile avgs and minimum avg
          # between 5 and 15 km from epicenter as established
          mobAvgNumb <- nrow(subset(eq_avgs, !is.nan(eq_avgs$Int)))
          if (!is.na(eq_avgs$Int[1]) &&
              eq_avgs$Int[1] > meanv0_10IntMin &&
              mobAvgNumb >= mobAvgNumbMin) {
            eq_IntNumbTot <- nrow(eq) # MDPs total
            eq_IntNumb55 <-
              nrow(subset(eq, eq$Dist <= 55)) # MDPs within 55 km
            
            # calculating earthquake depth using Rosetta algorithm
            cAng <- lm(eq_avgs$Int ~ eq_avgs$Dist)$coefficient[2]
            slope <- abs(as.numeric(cAng))
            eq_Dnext <- exp((slope - 0.087) / (-0.018))
            
            # If eq_Dnext is shallower than 5 km it will be set to 5 km
            if (eq_Dnext < 5) {
              eq_Dnext <- 5
            }
            
            # if eq_Dnext is deeper than 73 km it will be set to 73 km
            if (eq_Dnext > 73) {
              eq_Dnext <- 73
            }
            
            # Calculating intercept
            slopeInt <-
              summary(lm(eq_avgs$Int ~ eq_avgs$Dist))$coefficient[1]
            # Intercept correction due to W&C fault radius
            slopeIntCorr <- slopeInt - slope * Rwc
            
            # (re-)calculating intercept Mw and corrected one
            MwI <- 0.18 * log(eq_Dnext) + 0.56 * slopeInt + 1.44
            MwICorr <-
              0.18 * log(eq_Dnext) + 0.56 * slopeIntCorr + 1.44
          } else {
            doNotCalc <- TRUE
          }
        }
        
        if (!doNotCalc) {
          # calculating stdError
          stdError <-
            summary(lm(eq_avgs$Int ~ eq_avgs$Dist))$coefficients[4]
          roundStdError <- round(stdError, digits = 3)
          
          # considerate earthquake with established conditions only
          if (coef(lm(eq_avgs$Int ~ eq_avgs$Dist))[[2]] < 0 &&
              eq_IntNumb55 >= eqIntNumb55Min &&
              summary(lm(eq_avgs$Int ~ eq_avgs$Dist))$coefficients[4] <= maxStdError) {
            ## analysis on macroseismic field homogeneity on the basis of ----
            # intensities distribution along azimuth section (azEqInt algorithm)
            
            # create two dataframes with MDPs and their azimuthal direction
            # with respect to epicenter
            dut <- 1
            outPut <-
              data.frame(
                intensity = integer(),
                azimuth = double(),
                distance = double()
              )
            outPutDeg <-
              data.frame(
                intensity = integer(),
                azimuth = double(),
                distance = double()
              )
            
            while (dut <= nrow(eqraw3)) {
              outPut[nrow(outPut) + 1,] = c(eqraw3$MCS[dut],
                                            deg2rad(bearing(
                                              c(eqEpiLon, eqEpiLat),
                                              c(eqraw3$Lon[dut],  eqraw3$Lat[dut])
                                            )),
                                            eqraw3$Dist_epi[dut])
              outPutDeg[nrow(outPutDeg) + 1,] = c(eqraw3$MCS[dut],
                                                  bearing(
                                                    c(eqEpiLon, eqEpiLat),
                                                    c(eqraw3$Lon[dut],  eqraw3$Lat[dut])
                                                  ),
                                                  eqraw3$Dist_epi[dut])
              dut <- dut + 1
            }
            
            # distinguish MDPs in relation to their azimuth
            outPutDeg3 <- subset(outPutDeg, intensity <= 3)
            outPutDeg5 <-
              subset(outPutDeg, intensity > 3 & intensity <= 5)
            outPutDeg7 <-
              subset(outPutDeg, intensity > 5 & intensity <= 7)
            outPutDeg9 <-
              subset(outPutDeg, intensity > 7 & intensity <= 9)
            outPutDeg11 <- subset(outPutDeg, intensity > 9)
            
            outPutDegSel <-
              subset(outPutDeg, distance >= 10 & distance <= 55)
            
            outPutDegSel3 <- subset(outPutDegSel, intensity <= 3)
            outPutDegSel5 <-
              subset(outPutDegSel, intensity > 3 & intensity <= 5)
            outPutDegSel7 <-
              subset(outPutDegSel, intensity > 5 & intensity <= 7)
            outPutDegSel9 <-
              subset(outPutDegSel, intensity > 7 & intensity <= 9)
            outPutDegSel11 <- subset(outPutDegSel, intensity > 9)
            
            # evaluate if there is at least 1 intensity per interval
            # of 10 degrees 18 times (i.e. the intensities must be distributed
            # for a minimum of 180 azimuth degrees wherever located in the windrose)
            outCount <-
              data.frame(
                azimuthInterval = character(),
                azimuthCount = integer(),
                stringsAsFactors = FALSE
              )
            det <- 0
            while (det <= 170) {
              outCount[nrow(outCount) + 1, ] = list(paste(det, ">", det + 10, sep = ""), length(
                which(
                  outPutDegSel$azimuth >= det & outPutDegSel$azimuth < (det + 10)
                )
              ))
              det <- det + 10
            }
            detinv <- 0
            while (detinv >= -170) {
              outCount[nrow(outCount) + 1, ] = list(paste(detinv, "<", detinv - 10, sep = ""),
                                                    length(
                                                      which(
                                                        outPutDegSel$azimuth < detinv &
                                                          outPutDegSel$azimuth >= (detinv - 10)
                                                      )
                                                    ))
              detinv <- detinv - 10
            }
            resultAz <- length(which(outCount$azimuthCount > 0))
            
            # plotting ---------------------------------------------------
            # plot histogram of number of MDPs vs distance from epicenter
            if (histYES) {
              hist <- ggplot(eq_MDPs_number,
                             aes(x = AvgDist, y = nMDP)) +
                geom_bar(stat = "identity",
                         color = "blue",
                         fill = "transparent") +
                scale_x_continuous(
                  name = "Epicentral distance (kms)",
                  limits = c(0, 55),
                  breaks = seq(5, 50, by = 5)
                ) +
                scale_y_continuous(name = "Number of MDP") +
                labs(title = paste(
                  eqID,
                  " - ",
                  eqEpiArea,
                  "\n# MDP vs epicentral distance",
                  sep = ""
                )) +
                geom_text(aes(y = nMDP - nMDP * 10 / 100, label = nMDP))
              
              ggsave(
                filename = paste(outdir, "/", eqID, "_hist.png", sep = ""),
                hist,
                device = "png"
              )
            }
            if (resultAz >= resultAzMin) {
              if (liPlotto) {
                pdfOUT0_50 <- ggplot(eq_avgs, aes(Dist, Int)) +
                  ggtitle(paste("DBMI_", eqID, "|", eqDate, "-", eqEpiArea)) +
                  geom_point(aes(y = Int, colour = "Intensities")) +
                  scale_x_continuous(
                    name = "Epicentral distance (km)",
                    limits = c(0, 50),
                    breaks = seq(0, 50, by = 5)
                  ) +
                  theme(legend.position = "none") +
                  scale_y_continuous(
                    name = "Avg. intensities (MCS)",
                    limits = c(-10, 51),
                    breaks = seq(1, max(eq_avgs$Int, na.rm = TRUE) + 1)
                  ) +
                  coord_cartesian(xlim = c(0, 50),
                                  ylim = c(1, max(eq_avgs$Int, na.rm = TRUE) + 1)) +
                  geom_smooth(
                    method = "lm",
                    fullrange = TRUE,
                    se = TRUE,
                    color = "blue"
                  ) +
                  annotate(
                    "rect",
                    xmin = 0,
                    xmax = 33,
                    ymin = 1,
                    ymax = 2,
                    alpha = 0.6,
                    color = "grey",
                    fill = "white"
                  ) +
                  annotate(
                    geom = "text",
                    x = 1,
                    y = 1.6,
                    hjust = 0,
                    size = 3.5,
                    label = paste0(
                      "\nSteepness: ",
                      abs(round(cAng, digits = 3)),
                      "\nExpected depth: ",
                      round(eq_Dnext, digits = 1),
                      "\nSteepness std. error: ",
                      roundStdError,
                      "\nMDPs within 55 km: ",
                      eq_IntNumb55
                    )
                  ) +
                  annotate(
                    geom = "text",
                    x = 16,
                    y = 1.6,
                    hjust = 0,
                    size = 3.5,
                    label = paste0(
                      "\nIntercept value: ",
                      round(slopeInt, digits = 2),
                      "\nCorrected Intercept value: ",
                      round(slopeIntCorr, digits = 2),
                      "\nMw y-intercept (this work): ",
                      round(MwI, digits = 2),
                      "\nCorrected Mw y-intercept (this work): ",
                      round(MwICorr, digits = 2)
                    )
                  )
              }
              
              if (liSalvo) {
                ggsave(
                  filename =  paste(eqID,
                                    "-0_50",
                                    ".pdf",
                                    sep = ""),
                  plot = pdfOUT0_50,
                  device = "pdf",
                  path = outdir
                )
              }
              
              # exporting eq parameters to ouput csv file -------------------
              newline <-
                list(
                  eqID,
                  eqEpiArea,
                  eqDate,
                  eqMw,
                  eqDepDef,
                  eq_IntNumbTot,
                  eq_IntNumb55,
                  mobAvgNumb,
                  eqEpiLon,
                  eqEpiLat,
                  round(coef(
                    lm(eq_avgs$Int ~ eq_avgs$Dist)
                  )[[2]], digits = 3),
                  round(eq_Dnext, digits = 2),
                  round(summary(
                    lm(eq_avgs$Int ~ eq_avgs$Dist)
                  )$coefficients[4], digits = 4),
                  round(slopeInt, digits = 3),
                  round(slopeIntCorr, digits = 3),
                  round(MwI, digits = 2),
                  round(MwICorr, digits = 2),
                  resultAz,
                  Rwc,
                  paste("file://", outdir, "/", eqID, "-0_50.pdf", sep = "")
                )
              write.table(
                newline,
                file = outFileName,
                append = TRUE,
                sep = ",",
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE
              )
            }
          }
        }
      }
    }
  }
}