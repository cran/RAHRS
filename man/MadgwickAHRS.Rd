\name{MadgwickAHRS}
\alias{MadgwickAHRS}
\title{
Madgwick's AHRS algorithm
}
\description{
Implementation of Madgwick's AHRS algorithm.
}
\usage{
MadgwickAHRS(MSamplePeriod, MBeta, q, Gyroscope, Accelerometer, Magnetometer)
}
\arguments{
  \item{MSamplePeriod}{
Sample Period.
}
  \item{MBeta}{
Algorithm gain beta.
}
  \item{q}{
Quaternion.
}
  \item{Gyroscope}{
Gyroscope measurement.
}
  \item{Accelerometer}{
Accelerometer measurement.
}
  \item{Magnetometer}{
Magnetometer measurement.
}
}
\value{
Quaternion describing the Earth relative to the sensor.
}
\references{
S. O. H. Madgwick, An efficient orientation filter for inertial and 
inertial/magnetic sensor arrays, Technical report, 
University of. Bristol University, UK, 2010.
}
\author{
Jose Gama
}
\keyword{ ~Programming }
