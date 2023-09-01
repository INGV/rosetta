## ABOUT THIS WORK

The code in tis repository is an R (<https://cran.r-project.org>) implementation of the *"Rosetta"* algorithm described by Sbarra et al. 2023 (<https://doi.org/10.5194/nhess-23-1007-2023>). It calculates an estimation of hypocentral depth and magnitude of earthquakes based on Macroseismic Data Points (MDP) reported in the Database Macrosismico Italiano (DBMI) (<https://doi.org/10.13127/DBMI/DBMI15.4>).

It takes a CSV file with a list of DBMI IDs as input, downloads data from the internet, and calculates every earthquake's estimated depth and magnitude plotting mobile averaged macroseismic intensities for an epicentral distance of 0 to 55 km.

To increase its readability, the code is extensively commented.

**N.B.** The values shown in the S1 table of the article mentioned above were calculated based on DBMI version 15.2 (<https://doi.org/10.13.127/DBMI/DBMI15.2>). Instead, the script downloads the last available version from DBMI web services on the fly. There may be slight differences between versions, resulting in minor differences in the calculated values.

## BASIC USAGE

In main.R modify the parameters

`eqEventsFile` to point to input CSV file

and

`outdir` to indicate where to store results (CSV file and plots)

then run the code in R (or RStudio) environment.

The input file is a txt file with list of earthquakes to download File must be a single-column txt file with DBMI IDs one per row; see file `eventListTemplate.csv` for an example.

## ADVANCED USAGE {#sec-advanced-usage}

In addition to what is described in the previous paragraph, users can modify requirement filter values (file `main.R`, rows 40 -- 46). Defaults are set as reported in the article mentioned above.

*Minimum number of MDPs within 55 km from the epicenter (default = 30)*

`eqIntNumb55Min <- 30`

*Maximum standard error admitted for estimated attenuation steepness (default = 0.01)*

`maxStdError <- 0.01`

*The attenuation steepness must be calculated based on six or more averaged points; thus, at least 6 of the 10 rings must contain suitable MDPs (default = 6)*

`mobAvgNumbMin <- 6`

*The minimum average value of MDPs between 0 and 10 km (or fault R for Mw \>= 6.75) from the epicenter (default = 4)*

`meanv0_10IntMin <- 4`

*The MDPs falling 10--55 km from the epicenter must be distributed in an azimuthal range ≥180 degrees to consider only earthquakes with a "homogenous" macroseismic field. To this aim, the macroseismic field has been divided into 36 "slices" of 10 degrees, and the following parameter tells us how many of them must contain at least 1 MDP to be included in the count. (default = 18)*

`resultAzMin <- 18`

## OUTPUT DESCRIPTION

For each event that meets the stated requirements in the output directory, you'll find a plot, a histogram, and a line in a CSV file with the calculated parameters. Events that do not meet the requirements will be ignored and they will not generate output at all.

The plot shows the average values of the Macroseismic Data Points (MDPs) for every mobile window of distance from the epicenter (0 km - 10 km, 5 km - 15 km, 10 km - 20 km,...), the related linear regression, and a box with several calculated parameters.

The histogram shows the number of MDPs used in the calculation vs. epicentral distance windows.

### Description of CSV fields

| Field name                 | Description                                                                                    |
|-------------------|-----------------------------------------------------|
| DBMI Id                    | DBMI Id                                                                                        |
| Epicentral area            | Epicentral area name as for DBMI                                                               |
| Date                       | Date and time                                                                                  |
| Mw                         | Registered/estimanted local magnitude (as for DBMI)                                            |
| Instrumental depth         | Registered hypocenter depth (if any)                                                           |
| Total intensities          | Total number of reported MDPs                                                                  |
| Intensities within 55 km   | Number of MDPs in a radius of 55 km from epicenter                                             |
| Mobile avgs within 55 km   | Mobile averages calculated in a radius of 55 km from epicenter                                 |
| Epicenter Lon              | Epicenter longitude                                                                            |
| Epicenter Lat              | Epicenter latitude                                                                             |
| Angular coefficient        | Steepness of the attenuation curve                                                             |
| Calculated Depth           | Estimated depth (this work)                                                                    |
| Std. error                 | standard error of the estimated attenuation steepness                                          |
| Y-axis intercept           | Expected intensity at the epicenter                                                            |
| Y-axis corrected intercept | Expected intensity at the epicenter for earthquakes with a Mw \>= 6.75                         |
| Mw_Y intercept             | Expected Mw at the epicenter (this work)                                                       |
| Mw_I corrected intercept   | Expected Mw at the epicenter for earthquakes with a Mw \>= 6.75 (this work)                    |
| N. azimuth slices          | N. of azimuth slices of 10° with MDPs (see [Advanced Usage](#sec-advanced-usage), last point). |
| Fault radius (W&C)         | Seismogenic fault estimated radius following Wells & Coppersmith (1994)                        |
| Plot file                  | Location in the file system of the related plot file                                           |

For a detailed description of all the parameters shown in the CSV file, please refer to Sbarra et al. 2023 (<https://doi.org/10.5194/nhess-23-1007-2023>).

If you experience any problems, have any doubts or suggestions, please, [open an issue](https://github.com/INGV/rosetta/issues) on this GitHub repository.

## TO DO

-   Add implementation for HSIT (<https://doi.org/10.13127/hsit>) and CFTI (<https://doi.org/10.6092/ingv.it-cfti5>) data.

-   Implement as RScript

## References

1.  Guidoboni E., Ferrari G., Mariotti D., Comastri A., Tarabusi G., Sgattoni G., Valensise G. (2018). CFTI5Med, Catalogo dei Forti Terremoti in Italia (461 a.C.-1997) e nell'area Mediterranea (760 a.C.-1500) (Version 5). Istituto Nazionale di Geofisica e Vulcanologia (INGV). <https://doi.org/10.6092/ingv.it-cfti5>

2.  Locati M., Camassi R., Rovida A., Ercolani E., Bernardini F., Castelli V., Caracciolo C.H., Tertulliani A., Rossi A., Azzaro R., D'Amico S., Antonucci A. (2019). Database Macrosismico Italiano (DBMI15), versione 2.0. Istituto Nazionale di Geofisica e Vulcanologia (INGV). <https://doi.org/10.13127/DBMI/DBMI15.2>

3.  Locati M., Camassi R., Rovida A., Ercolani E., Bernardini F., Castelli V., Caracciolo C.H., Tertulliani A., Rossi A., Azzaro R., D'Amico S., Antonucci A. (2022). Database Macrosismico Italiano (DBMI15), versione 4.0. Istituto Nazionale di Geofisica e Vulcanologia (INGV). <https://doi.org/10.13127/DBMI/DBMI15.4>

4.  Sbarra, P., Burrato, P., De Rubeis, V., Tosi, P., Valensise, G., Vallone, R., and Vannoli, P.: Inferring the depth and magnitude of pre-instrumental earthquakes from intensity attenuation curves, Nat. Hazards Earth Syst. Sci., 23, 1007--1028, <https://doi.org/10.5194/nhess-23-1007-2023,> 2023

5.  Tosi P., De Rubeis V., Sbarra P., Sorrentino D. (2007). Hai Sentito Il Terremoto (HSIT). Istituto Nazionale di Geofisica e Vulcanologia (INGV). <https://doi.org/10.13127/hsit>

6.  Donald L. Wells, Kevin J. Coppersmith; New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. *Bulletin of the Seismological Society of America* 1994;; 84 (4): 974--1002. doi: <https://doi.org/10.1785/BSSA0840040974>
