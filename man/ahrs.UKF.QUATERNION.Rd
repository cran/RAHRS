\name{ahrs.UKF.QUATERNION}
\alias{ahrs.UKF.QUATERNION}
\title{
Estimate the current quaternion attitude by UKF-based AHRS algorithm
}
\description{
Estimates the current quaternion attitude with an implementation of the 
UKF-based AHRS algorithm based on measurements from three-component 
accelerometer with orthogonal axes, vector magnetometer 
and three-axis gyroscope.
}
\usage{
ahrs.UKF.QUATERNION(Filter, Sensors, Parameters)
}
\arguments{
  \item{Filter}{data structure for Unscented Kalman Filter
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
  \item{Parameters}{AHRS Parameters
   Parameters.mn      Magnetic Field Vector In Navigation Frame [3x1], |m|
   = 1
   Parameters.an      Acceleration vector In Navigation Frame [3x1], g
   Parameters.dt      Sampling period, 1/Hz
}
}
\value{data structure for Unscented Kalman Filter
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
