#' @title Example Data Frame for adsorption, desorption, and supply parameter
#' @description User is advised to prepare the data as suggested in the example to fit adsorption or desorption data pertaining to a particular nutrient, e.g., phosphorus (as phosphate), sulphur (as sulphate), micronutrient cations or anions, silicon (as silicate), etc. in soil to the functions namely, Freundlich_A, Freundlich_D, Langmuir, and SP. However, to use these functions, researchers should carry out the adsorption-desorption study in the method described as follows: W g (e.g., 2 g) of processed soil samples are taken in a series of 50 mL centrifuge tubes (polypropylene) each containing V mL (e.g., 20 mL) 0.1 M sodium chloride solution with different levels (e.g., 10, 20, 40, 60, … mg/L) of the nutrient under consideration (e.g., P). The tubes are shaken continuously for a defined period (e.g., 24 h) using a mechanical shaker. Immediately after that, each of them is centrifuged (e.g., at 10000 rpm) for sufficient time (e.g., 10-12 min), and a certain volume (e.g., 15 mL) of clear aliquot is pipetted out and filtered through a Whatman No. 42 filter paper. Then the nutrient content (e.g., P) in the filtrate is estimated by standard procedures. The difference between the quantities of the nutrient in the bathing solution before (initial concentration) and after equilibration (equilibrium or final concentration) is taken as the amount of the nutrient adsorbed by the soil from bathing solution. For desorption study, once supernatant is removed after completion of the adsorption step, a certain volume (e.g., 15 mL) of 0.1 M sodium chloride solution is added in the same sample to reach a final volume (Vfinal). For convenience in calculation, Vfinal should be equal to V (as taken initially). Generally, soil sample with the highest level of added nutrient from the adsorption study is used in desorption experiment. With the help of mechanical shaker, the soil is re-suspended and equilibrated for 12 h. After shaking, the tubes are centrifuged (e.g., at 10000 rpm) for sufficient time (e.g., 10–12 min), and the same volume (as in the previous step, e.g., 15 mL) of clear supernatant solution is pipetted out for determination of the concentration of the nutrient. This process can be repeated multiple times depending on the objective and convenience of the researcher. The difference between equilibrium concentration at a particular desorption step and equilibrium concentration previous to that desorption step, multiplied by solution:soil ratio (e.g., 20:2) is considered as desorbed amount of the nutrient. The adsorbed amount remaining after each desorption step can be computed by subtracting the desorbed amount at a particular step from the amount of adsorbed nutrient present before that desorption step.
#' @usage df_sordes
#' @format Write the following notations on the spreadsheet:
#' Initial_conc for Initial concentration (mg/L) of the added element, e.g., phosphorus
#' Equilibrium_conc for Final or equilibrium concentration (mg/L) after adsorption of the same element
#' Cf_Des for Final or equilibrium concentration (mg/L) after desorption of the same element

#' @export
df_sordes = read.table(text = "Initial_conc	Equilibrium_conc	Cf_Des
10	0.619589386	13.62847222
20	1.633333333	7.477678571
40	17.2781449	6.051587302
60	35.18963333	4.836309524", header = TRUE)

#' @title Freundlich isotherm fitted to adsorption data
#'
#' @description The linear form of Freundlich adsorption isotherm can be fitted to adsorption data to find out the empirical constants of Freundlich equation (Freundlich, 1926).
#' @importFrom graphics abline legend
#' @importFrom stats coef lm
#'
#' @usage Freundlich_A(W = W, V = V, Ci = Ci, Cf = Cf,...)
#'
#' @param W Mass of soil sample (g)
#' @param V Volume of extractant solution (mL)
#' @param Ci Initial concentration (mg/L) of the added element, e.g., phosphorus
#' @param Cf Final or equilibrium concentration (mg/L) after adsorption of the same element
#' @param ... Any other argument that can be passed to base plot
#' @return a - empirical constant (unitless)
#'  1/n - empirical constant (unitless)
#'
#' @references Freundlich, H., 1926. Colloid and Capillary Chemistry. London: Methuen.
#' @details Freundlich equation/isotherm is used to study the adsorption behavior of any element in soils (Freundlich, 1926). It is of the general form: x/m= ac^(1/n), where ‘x’ is the amount of the adsorbate (e.g., P, Ni, etc.) adsorbed on ‘m’ amount of the adsorbent (e.g., soil), and ‘c’ is the equilibrium concentration of the adsorbate. This adsorption model helps to understand the relationship between quantity of any element adsorbed per unit soil weight and their concentration in soil solution. Freundlich equation does not predict or include maximum adsorption capacity of soil, so better suited to dilute solutions (of the element under consideration) in contact with the adsorbent. The ‘a’ and ‘1/n’ are two empirical-constants, sensitive to the given adsorbent-adsorbate system and temperature. Freundlich ‘a’ indicates the binding affinity of any element in soil. The ‘1/n’ is the exponent of the equilibrium concentration term in the Freundlich equation, and is always < 1 indicating that the energy of adsorption decreases logarithmically as the fraction of adsorbent-surface covered by the adsorbate increases with increasing equilibrium concentration.
#' @export
#' @examples
#' with(data = df_sordes, Freundlich_A(W = 2, V = 20, Ci = Initial_conc,
#' Cf = Equilibrium_conc))

Freundlich_A <- function(W = W, V = V, Ci = Ci, Cf = Cf,...){

  p_adsorbed <- (Ci - Cf )*(V/W)
  logC <- log10(Cf)
  logX <- log10(p_adsorbed)

  myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
  plot(x, y, xlab=xlab, ylab=ylab, ...)
}
myplot(x = logC, y = logX, ...)
abline(lm(logX ~ logC), col = "blue")
fit.slr <- lm(logX ~ logC)

##Round the coefficients for better output
cf <- round(coef(fit.slr), 3)

##Make the regression equation
eq <- paste0("y = ", cf[1],
             ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")

#Get the R2 value
Rsq = format(summary(fit.slr)$r.squared,digits=3)

##Printing of the equation and R2 on the plot
legend("topleft", legend = c(eq, as.expression(bquote(R^2 == .(Rsq)))), bty = "n")
a <- 10^coef(fit.slr)[[1]]
onebyn <- coef(fit.slr)[[2]]
return(list(a = a, onebyn = onebyn))
}

#' @title Freundlich isotherm for desorption data in soil
#'
#' @description This function fits Freundlich adsorption isotherm to data pertaining to desorption of an adsorbed nutrient, e.g., phosphorus (as phosphate), sulphur (as sulphate), micronutrient cations or anions, silicon (as silicate), etc. in soil.
#' @importFrom graphics abline legend
#' @importFrom stats coef lm
#' @usage Freundlich_D(W = W, V = V, Ci = Ci, Cf = Cf, Cf_Des = Cf_Des, Vres= Vres,
#' Vfinal = Vfinal, onebyn_ads = onebyn_ads, ...)
#'
#' @param W Mass of soil sample (g)
#' @param V Volume of extractant solution (mL)
#' @param Ci Initial concentration (mg/L) of the added element, e.g., phosphorus
#' @param Cf Final or equilibrium concentration (mg/L) after adsorption of the same element
#' @param Cf_Des Final or equilibrium concentration (mg/L) after each desorption step of the same element
#' @param Vres Volume of carried over solution from previous adsorption or desorption step (mL)
#' @param Vfinal Final volume of the solution for each desorption step (mL)
#' @param onebyn_ads Values of ‘1/n’ obtained from Freundlich_A
#' @param ... Any other argument that can be passed to base plot
#'
#' @return a - empirical constant (unitless)
#' 1/n - empirical constant (unitless)
#'
#' @references Barman, M., Datta, S.P., Rattan, R.K., Meena, M.C., 2013. Sorption and desorption of nickel in soils in relation to its availability to plants. Agrochimica LVII (3), 235–249.
#'
#' @export
#' @examples
#' with(data = df_sordes, Freundlich_D(W = 2, V = 20, Ci = Initial_conc, Cf = Equilibrium_conc,
#' Cf_Des = Cf_Des, Vres= 5, Vfinal = 20, onebyn_ads = 0.2056263))
#'
Freundlich_D <- function(W = W, V = V, Ci = Ci, Cf = Cf, Cf_Des = Cf_Des, Vres= Vres,
               Vfinal = Vfinal, onebyn_ads = onebyn_ads, ...){
  ##Initialize the calculation
  Ads <- (Ci - Cf)*(V/W)
  Ads_Des <- numeric(length(Cf_Des))
  Ci_x <- max(Cf, na.rm=T)/(Vfinal/Vres)
  Ads_Des[[1]] <- Ads_x <- max(Ads, na.rm=T) - (Cf_Des[[1]] - Ci_x)*(Vfinal/W)
  ##For loop
  for (i in seq_along(Cf_Des)[-1]) {
    Des_x <- (Cf_Des[[i]] - Cf_Des[[i - 1]]/(Vfinal/Vres))*(Vfinal/W)
    Ads_Des[[i]] <- Ads_x <- Ads_x - Des_x
  }
  ## processing
  log_Cf_Des <- log10(Cf_Des)
  log_Ads_Des <- log10(Ads_Des)

  myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
    plot(x, y, xlab=xlab, ylab=ylab, ...)
  }
  myplot(x = log_Cf_Des, y = log_Ads_Des, ...)
  abline(lm(log_Ads_Des ~ log_Cf_Des), col = "blue")
  fit.slr <- lm(log_Ads_Des ~ log_Cf_Des)

  ##Round the coefficients for better output
  cf <- round(coef(fit.slr), 4)

  ##Make the regression equation
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")

  #Get the R2 value
  Rsq = format(summary(fit.slr)$r.squared, digits=3)

  ##Printing of the equation and R2 on the plot
  legend("topleft", legend = c(eq, as.expression(bquote(R^2 == .(Rsq)))), bty = "n")
  a <- 10^coef(fit.slr)[[1]]
  onebyn <- coef(fit.slr)[[2]]
  DI <- onebyn_ads/onebyn
  return(list(a = a, onebyn = onebyn, DI = DI))
}

#' @title Langmuir isotherm fitted to adsorption data
#' @description The linear form of Langmuir adsorption isotherm can be fitted to adsorption data to find out adsorption maxima, affinity coefficient, and maximum buffering capacity (Langmuir, 1918).
#' @importFrom graphics abline legend
#' @importFrom stats coef lm
#' @usage Langmuir(W = W, V = V, Ci = Ci, Cf = Cf,...)
#' @param W Mass of soil sample (g)
#' @param V Volume of extractant solution (mL)
#' @param Ci Initial concentration (mg/L) of the added element, e.g., phosphorus
#' @param Cf Final or equilibrium concentration (mg/L) after adsorption of the same element
#' @param ... Any other argument that can be passed to base plot
#' @return b - Maximum monolayer adsorption or adsorption maxima (mg/kg)
#' k - Constant related to binding energy or the affinity coefficient for the nutrient, e.g., P (L/mg)
#' MBC - Maximum buffering capacity of soil for the nutrient under consideration (L/kg)
#'
#' @references Havlin, J.L., Tisdale, S.L., Nelson, W.L., Beaton, J.D., 2016. Phosphorus. In: Soil fertility and fertilizers. Pearson Education India, pp.160–198.
#' Holford, I.C.R., Mattingly, G.E.G., 1976. Phosphate adsorption and plant availability of phosphate. Plant and Soil 44, 377–389. https://doi.org/10.1007/BF00015889
#' Langmuir, I., 1918. The adsorption of gases on plane surfaces of glass, mica, and platinum. Journal of American Chemical Society 40, 1361–1403. https://doi.org/10.1021/ja02242a004
#' Shirvani, M., Shariatmadari, H., Kalbasi, M., 2005. Phosphorus buffering capacity indices as related to soil properties and plant uptake. Journal of Plant Nutrition, 28, 537–550. https://doi.org/10.1081/PLN-200049235
#' @export
#' @examples
#' with(data = df_sordes, Langmuir(W = 2, V = 20, Ci = Initial_conc, Cf = Equilibrium_conc))

Langmuir <- function(W = W, V = V, Ci = Ci, Cf = Cf,...){

  p_adsorbed <- (Ci - Cf)*(V/W)
  CbyX <- Cf/p_adsorbed

  myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
    plot(x, y, xlab=xlab, ylab=ylab, ...)
  }
  myplot(x = Cf, y = CbyX, ...)
  abline(lm(CbyX ~ Cf), col = "blue")
  fit.slr <- lm(CbyX ~ Cf)

  ##Round the coefficients for better output
  cf <- round(coef(fit.slr), 4)

  ##Make the regression equation
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")

  #Get the R2 value
  Rsq = format(summary(fit.slr)$r.squared, digits=3)

  ##Printing of the equation and R2 on the plot
  legend("topleft", legend = c(eq, as.expression(bquote(R^2 == .(Rsq)))), bty = "n")
  b <- 1/coef(fit.slr)[[2]]
  k <- coef(fit.slr)[[2]]/coef(fit.slr)[[1]]
  MBC <- b*k
  return(list(b = b, k = k, MBC = MBC))
}

#' @title Supply parameter of phosphorus in soil
#' @description This function generates the supply parameter (SP) of phosphorus in soil as described by Khasawneh and Copeland (1973).
#' @importFrom graphics abline legend
#' @importFrom stats coef lm
#' @usage SP(W = W, V = V, Ci = Ci, Cf = Cf,...)
#' @param W Mass of soil sample (g)
#' @param V Volume of extractant solution (mL)
#' @param Ci Initial concentration (mg/L) of the added element, e.g., phosphorus
#' @param Cf Final or equilibrium concentration (mg/L) after adsorption of the same element
#' @param ... Any other argument that can be passed to base plot
#'
#' @return SP - Supply parameter (mg^0.5)/(kg L)^0.25
#'
#' @references Khasawneh, F.E., 1971. Solution ion activity and plant growth. Soil Science Society of America Proceedings 35, 426–436.
#' Khasawneh, F.E., Copeland, J.P., 1973. Cotton root growth and uptake of nutrients: relation of phosphorus uptake to quantity, intensity, and buffering capacity. Soil Science Society of America Proceedings 37, 250–254.
#' @export
#' @examples
#' with(data = df_sordes, SP(W = 2, V = 20, Ci = Initial_conc, Cf = Equilibrium_conc,
#' col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5), pch = 16, cex = 1))

SP <- function(W = W, V = V, Ci = Ci, Cf = Cf,...){

  p_adsorbed <- (Ci - Cf)*(V/W)
  onebyC <- 1/Cf
  onebyX <- 1/p_adsorbed

myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
  plot(x, y, xlab=xlab, ylab=ylab, ...)
}
myplot(x = onebyC, y = onebyX, ...)
abline(lm(onebyX ~ onebyC), col = "blue")
fit.slr <- lm(onebyX ~ onebyC)

##Round the coefficients for better output
cf <- round(coef(fit.slr), 4)

##Make the regression equation
eq <- paste0("y = ", cf[1],
             ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")

#Get the R2 value
Rsq = format(summary(fit.slr)$r.squared, digits=3)

##Printing of the equation and R2 on the plot
legend("topleft", legend = c(eq, as.expression(bquote(R^2 == .(Rsq)))), bty = "n")

k1 <- 1/coef(fit.slr)[[1]]
k2 <- k1*coef(fit.slr)[[2]]
SP <- (k1*k2)^(1/4)*sqrt(Cf*p_adsorbed/(k1*k2))
return(list(SP = SP))
}
