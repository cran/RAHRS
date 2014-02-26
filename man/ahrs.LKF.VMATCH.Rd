\name{ahrs.LKF.VMATCH}
\alias{ahrs.LKF.VMATCH}
\title{
Quaternion estimation with vector matching and Kalman filter
}
\description{
Attitude quaternion estimation by means of complementary
Kalman filter.
}
\usage{
ahrs.LKF.VMATCH(Filter, Sensors, q, Parameters)
}
\arguments{
  \item{Filter}{data structure for Linear Kalman Filter
   Filter.x   State vector [3x1]
   Filter.P   Covariance matrix [3x3]
   Filter.Q   System noise matrix [3x3]
   Filter.R   Measurement noise matrix [6x6]
}
  \item{Sensors}{sensors data structure
   Sensors.w    current calibrated gyroscope measurement [3x1], rad/sec
   Sensors.a    current calibrated accelerometer measurement [3x1], g
   Sensors.m    current calibrated magnetometer measurement [3x1], |m| = 1
}
  \item{q}{ quaternion}
  \item{Parameters}{AHRS Parameters
   Parameters.mn      Magnetic Field Vector In Navigation Frame [3x1], |m|
   = 1
   Parameters.an      Acceleration vector In Navigation Frame [3x1], g
   Parameters.dt      Sampling period, 1/Hz
}
}
\value{
	\item{Filter}{data structure for Linear Kalman Filter}
	\item{Q}{Correct quaternion}
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
