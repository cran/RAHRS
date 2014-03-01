\name{MadgwickAHRSupdate}
\alias{MadgwickAHRSupdate}
\alias{MadgwickAHRSupdate2}
\alias{MadgwickAHRSupdateDbl}
\title{
Madgwick's AHRS update method
}
\description{
Madgwick's algorithm for AHRS update method.
}
\usage{
MadgwickAHRSupdate(gxi, gyi, gzi, axi, ayi, azi, mxi, myi, mzi, sampleFreqi, 
betai, q0i, q1i, q2i, gz)
}
\arguments{
  \item{gxi}{
Gyroscope x axis measurement in radians/s.
}
  \item{gyi}{
Gyroscope y axis measurement in radians/s.
}
  \item{gzi}{
Gyroscope z axis measurement in radians/s.
}
  \item{axi}{
Accelerometer x axis measurement in any calibrated units.
}
  \item{ayi}{
Accelerometer y axis measurement in any calibrated units.
}
  \item{azi}{
Accelerometer z axis measurement in any calibrated units.
}
  \item{mxi}{
Magnetometer x axis measurement in any calibrated units.
}
  \item{myi}{
Magnetometer y axis measurement in any calibrated units.
}
  \item{mzi}{
Magnetometer z axis measurement in any calibrated units.
}
  \item{sampleFreqi}{
Sample frequency.
}
  \item{betai}{
Algorithm gain beta.
}
  \item{q0i}{
Quaternion[0]
}
  \item{q1i}{
Quaternion[1]
}
  \item{q2i}{
Quaternion[2]
}
  \item{gz}{
Quaternion[3]
}
}
\value{
output quaternion describing the Earth relative to the sensor.
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
