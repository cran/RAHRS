\name{ekf.calibration.indirect}
\alias{ekf.calibration.indirect}
\title{
Calibration coefs by complementary EKF
}
\description{
Estimation of the calibration coefs by complementary EKF
}
\usage{
ekf.calibration.indirect(m, initMean = NA)
}
\arguments{
  \item{m}{ Calibration data, recorded while rotating corresponding sensor in 3D
}
  \item{initMean}{ Initial guess for coefs
}
}
\value{
      coefs[1x12]    vector of sensor's calibration coeffs
      m_    calibrated data
      tr_    Covariance matrix trace
}
\references{
Vlad Maximov, 2012
Scalar Calibration of Vector accelerometers and magnetometers,
GyroLib documentation
}
\author{
Jose Gama
}
\keyword{ ~Programming }
