\name{ahrs.EKF.QUATERNION}
\alias{ahrs.EKF.QUATERNION}
\title{
EKF-based AHRS algorithm
}
\description{
Implementation of the EKF-based AHRS algorithm based on measurements from 
three-component accelerometer with orthogonal axes, vector magnetometer and 
three-axis gyroscope. Estimates the current quaternion attitude.
}
\usage{
ahrs.EKF.QUATERNION(Filter, Sensors, Parameters)
}
\arguments{
  \item{Filter}{ data structure for Extended Kalman Filter }
  \item{Sensors}{ calibrated gyroscope, accelerometer and magnetometer 
  	  measurements
    Sensors.w    current calibrated gyroscope measurement [3x1], rad/sec
    Sensors.a    current calibrated accelerometer measurement [3x1], g
    Sensors.m    current calibrated magnetometer measurement [3x1], |m| = 1
  }
  \item{Parameters}{ AHRS Parameters
    Parameters.mn      Magnetic Field Vector In Navigation Frame [3x1], |m|= 1
    Parameters.an      Acceleration vector In Navigation Frame [3x1], g
    Parameters.dt      Sampling period, 1/Hz
  }
}
\value{
   Q    estimated attitude quaternion [1x4]
   Filter    data structure for an Extended Kalman Filter
   dw    estimated gyroscopes bias [1x3]
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
