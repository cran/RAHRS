\name{svd.calibration}
\alias{svd.calibration}
\title{
calibration coefs by Merayo's technique
}
\description{
Performs the estimation of the calibration coefs by complementary EKF 
using Merayo technique with a non iterative algoritm
}
\usage{
svd.calibration(X)
}
\arguments{
  \item{X}{Calibration data, recorded while rotating corresponding sensor in 3D
}
}
\value{
      coefs[1x12] - vector of sensor's calibration coeffs
      X_ - calibrated data
}
\references{
using Merayo technique with a non iterative algoritm
J.Merayo et al. "Scalar calibration of vector magnemoters"
Meas. Sci. Technol. 11 (2000) 120-132.

Alain Barraud, Suzanne Lesecq 2008
MgnCalibration - magnetometer calibration from a measurements data set.
\url{http://www.mathworks.com/matlabcentral/fileexchange/23398-magnetometers-calibration/content/MgnCalibration.m}
}
\author{
Jose Gama
}
\keyword{ ~Programming }
