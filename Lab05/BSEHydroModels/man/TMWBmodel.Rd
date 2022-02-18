\name{TMWBmodel}
\alias{TMWBmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TMWBmodel(TMWB = TMWB, fcres = 0.25, SFTmp = 0, bmlt6 = 2.5, bmlt12 = 1, Tlag = 0.5, AWCval = 200, Slope = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TMWB}{
%%     ~~Describe \code{TMWB} here~~
}
  \item{fcres}{
%%     ~~Describe \code{fcres} here~~
}
  \item{SFTmp}{
%%     ~~Describe \code{SFTmp} here~~
}
  \item{bmlt6}{
%%     ~~Describe \code{bmlt6} here~~
}
  \item{bmlt12}{
%%     ~~Describe \code{bmlt12} here~~
}
  \item{Tlag}{
%%     ~~Describe \code{Tlag} here~~
}
  \item{AWCval}{
%%     ~~Describe \code{AWCval} here~~
}
  \item{Slope}{
%%     ~~Describe \code{Slope} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (TMWB = TMWB, fcres = 0.25, SFTmp = 0, bmlt6 = 2.5, 
    bmlt12 = 1, Tlag = 0.5, AWCval = 200, Slope = 0) 
{
    attach(TMWB)
    SNO_Energy = SnowMelt(date, P, MaxTemp - 3, MinTemp - 3, 
        myflowgage$declat, slope = Slope, aspect = 0, tempHt = 1, 
        windHt = 2, groundAlbedo = 0.25, SurfEmissiv = 0.95, 
        windSp = 2, forest = 0, startingSnowDepth_m = 0, startingSnowDensity_kg_m3 = 450)
    detach(TMWB)
    TMWB$SNO = SNO_Energy$SnowWaterEq_mm
    TMWB$SNOmlt = SNO_Energy$SnowMelt_mm
    attach(TMWB)
    TMWB$Albedo = 0.23
    TMWB$Albedo[TMWB$SNO > 0] = 0.95
    PET = PET_fromTemp(Jday = (1 + as.POSIXlt(date)$yday), Tmax_C = MaxTemp, 
        Tmin_C = MinTemp, lat_radians = myflowgage$declat * pi/180) * 
        1000
    TMWB$PET = PET
    detach(TMWB)
    rm(list = c("PET"))
    TMWB$AWC = AWCval
    TMWB$dP = 0
    TMWB$ET = 0
    TMWB$AW = 0
    TMWB$Excess = 0
    attach(TMWB)
    for (t in 2:length(AW)) {
        ET[t] = min(AW[t - 1], PET[t])
        ET[t] = (AW[t - 1]/AWC[t - 1]) * PET[t]
        if (AvgTemp[t] >= SFTmp) {
            dP[t] = P[t] - ET[t] + SNOmlt[t]
        }
        else {
            dP[t] = ET[t]
        }
        if (dP[t] <= 0) {
            values <- soildrying(AW[t - 1], dP[t], AWC[t])
        }
        else if ((dP[t] > 0) & (AW[t - 1] + dP[t]) <= AWC[t]) {
            values <- soilwetting(AW[t - 1], dP[t], AWC[t])
        }
        else {
            values <- soil_wetting_above_capacity(AW[t - 1], 
                dP[t], AWC[t])
        }
        AW[t] <- values[1]
        Excess[t] <- values[2]
    }
    TMWB$AW = AW
    TMWB$Excess = Excess
    TMWB$dP = dP
    TMWB$ET = ET
    detach(TMWB)
    rm(list = c("AW", "Excess", "dP", "ET"))
    TMWB$Qpred = NA
    TMWB$Qpred[1] = 0
    TMWB$S = NA
    TMWB$S[1] = 0
    attach(TMWB)
    for (t in 2:length(Qpred)) {
        S[t] = S[t - 1] + Excess[t]
        Qpred[t] = fcres * S[t]
        S[t] = S[t] - Qpred[t]
    }
    TMWB$S = S
    TMWB$Qpred = Qpred
    detach(TMWB)
    rm(list = c("S", "Qpred"))
    return(TMWB)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
