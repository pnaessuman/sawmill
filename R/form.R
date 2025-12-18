##' Calculate stem volume based on height, diameter and species
##' 
##' This function uses the approximation function of Franz et al. (1973) to
##' calculate the stem volume based on the species specific form number. It is
##' estimated in two steps: a first level of functions is derived by species and
##' diameter, then this set of functions enters the final approximation function
##' which depends on height alone.
##' 
##' At the moment, only beech and spruce are implemented.
##' @title Stem volume
##' @param d numeric vector with stem diameters in 1.3 m height in cm
##' @param h numeric vectors with height in m
##' @param species character
##' @return numeric vector with volume of the trees in cubic metres
##' @examples
##' volume(30:35, 20:25, "Fagus sylvatica")
##' @export
volume <- function(d, h, species = "Fagus sylvatica") {

  if (!is.numeric(d) | !is.numeric(h)) {
    stop("d and h have to be numeric")
  }

  if (length(d) != length(h)) {
    stop("d and h have to have the same length")
  }

  species_list <- c(fagus_sylvatica = "Fagus sylvatica",
                    picea_abies = "Picea abies",
                    pinus_sylvestris = "Pinus sylvestris",
                    abies_alba = "Abies alba",
                    larix_decidua = "Larix decidua",
                    pseudotsuga_menziesii = "Pseudotsuga menziesii",
                    quercus_robur = "Quercus_robur",
                    alnus_glutinosa = "Alnus glutinosa",
                    fraxinus_excelsior = "Fraxinus excelsior",
                    pinus_cembra = "Pinus cembra",
                    nonexplicit_broadleaved = "Sonst. LH",
                    nonexplicit_conifer = "Sonst. NH")
  
  if (!any(species_list == species)) {
    stop("unknown species")
  }
  
  species_coefs <- list(
    picea_abies = c(
      k11  = -3.5962
      ,k12 = 1.8021
      ,k13 = -0.2882
      ,k21 = 1.0625
      ,k22 = -0.1290
      ,k23 = 0.0353
      ,k31 = 0.1423
      ,k32 = -0.0583
      ,k33 = 0.0046)
    ,fagus_sylvatica = c(
      k11  = -2.7284
      ,k12 = 0.8376
      ,k13 = -0.1058
      ,k21 = 1.6228
      ,k22 = -0.2148
      ,k23 = 0.0289
      ,k31 = -0.0880
      ,k32 = 0.0326
      ,k33 = -0.0045
      )
    ,pinus_sylvestris = c(
      k11  = -5.80915
      ,k12 = 3.387
      ,k13 = -0.494392
      ,k21 = 3.67116 
      ,k22 = -1.83211
      ,k23 = 0.273999 
      ,k31 = -0.459282
      ,k32 = 0.29989
      ,k33 = -0.0444931
    )
    ,abies_alba = c(
      k11  = -7.41365
      ,k12 = 3.33667
      ,k13 = -0.426419
      ,k21 = 4.00998 
      ,k22 = -1.39533
      ,k23 = 0.165198
      ,k31 = -0.321612
      ,k32 = 0.14401
      ,k33 = -0.0165461
    )
    ,larix_decidua = c(
      k11  = -9.26182 
      ,k12 = 4.75438
      ,k13 = -0.672495
      ,k21 = 5.17159 
      ,k22 = -2.27654 
      ,k23 = 0.311633
      ,k31 = -0.555379
      ,k32 = 0.302799 
      ,k33 = -0.041251
    )
    , pseudotsuga_menziesii = c(
      k11  = -12.5017
      ,k12 = 6.62441 
      ,k13 = -0.911185
      ,k21 = 7.27277 
      ,k22 = -3.58346 
      ,k23 = 0.489149
      ,k31 = -0.87715
      ,k32 = 0.515586
      ,k33 = -0.0714395
    )
  ,quercus_robur = c(
    k11  = -3.06118
    ,k12 = 1.45506
    ,k13 = -0.19992 
    ,k21 = 1.93898
    ,k22 = 0.112653
    ,k23 = -0.165102
    ,k31 =  0.120127
    ,k32 =-0.0202543
    ,k33 =-0.639727   
    )
  ,alnus_glutinosa = c(
    k11  =  -5.98031
      ,k12 = 2.85905 
      ,k13 =  -0.3374
      ,k21 = 3.78395 
      ,k22 = 0.133661
      ,k23 =- 0.540955
      ,k31 =0.296957
      ,k32 = -0.0385165 
      ,k33 =-1.47316 
    )
  
  cl2 <- species_coefs[[names(species_list)[which(species_list == species)]]]

  cl1f <- function(d) {
    unname(c(
    cl2["k11"] + cl2["k12"] * log(d) + cl2["k13"] * log(d)^2
    ,cl2["k21"] + cl2["k22"] * log(d) + cl2["k23"] * log(d)^2
    ,cl2["k31"] + cl2["k32"] * log(d) + cl2["k33"] * log(d)^2
    ))
  }

  cl1 <- lapply(d, cl1f)

  formhf <- function(h, coef) {
    exp(coef[1] + coef[2] * log(h) + coef[3] * log(h)^2)
  }

  formh <- mapply(formhf, h, cl1)

  d^2 * pi/40000 * formh
}
