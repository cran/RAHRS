\name{ahrs.LKF.QUATERNION}
\alias{ahrs.LKF.QUATERNION}
\title{
LKF-based AHRS algorithm
}
\description{
Implementation of the LKF-based AHRS algorithm based on measurements from 
three-component accelerometer with orthogonal axes, vector magnetometer and 
three-axis gyroscope. Estimates the current quaternion attitude.
}
\usage{
ahrs.LKF.QUATERNION(Filter, Sensors, q, Parameters, dw)
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
  \item{dw}{ angular rate}
}
\value{
	\item{Filter}{ Data structure for Linear Kalman Filter}
	\item{Q}{ Correct quaternion}
	\item{dw}{ Correct angular rate}
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
