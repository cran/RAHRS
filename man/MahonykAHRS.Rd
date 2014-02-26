\name{MahonykAHRS}
\alias{MahonykAHRS}
\title{
Mahony's AHRS algorithm
}
\description{
Implementation of Mahony's AHRS algorithm.
}
\usage{
MahonykAHRS(MSamplePeriod, Kp = 2, Ki = 0.005, q, Gyroscope, Accelerometer, Magnetometer)
}
\arguments{
  \item{MSamplePeriod}{
Sample Period.
}
  \item{Kp}{
Constant Kp.
}
  \item{Ki}{
Constant Ki.
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

R. Mahony, T. Hamel, and J.-M. Pflimlin. Nonlinear complementary filters on the
special orthogonal group. Automatic Control, IEEE Transactions on, 53(5):1203 –1218,
june 2008.
}
\author{
Jose Gama
}
\keyword{ ~Programming }
