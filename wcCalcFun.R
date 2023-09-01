# Fault parameters calculation as for Wells and Coppersmith (1994)
# modified by Pierfrancesco (Paco) Burrato

wcCalc <- function(eqMag) {
  faultArea <- 10 ** (-3.49 + (0.91 * eqMag))
  faultRadius <- sqrt((faultArea/pi))
  faultLength <- 10 ** (-3.22 + (0.69 * eqMag))
  faultWidth <- faultArea/faultLength
  faultProjWidth <- faultWidth/(sqrt(2))
  faultProjArea <- faultLength*(faultWidth/(sqrt(2)))
  faultProjRadius <- sqrt(faultProjArea/pi)
  save(
    faultArea,
    faultRadius,
    faultLength,
    faultWidth,
    faultProjWidth,
    faultProjArea,
    faultProjRadius,
    file = "faultMetaWC.RData"
  )
}