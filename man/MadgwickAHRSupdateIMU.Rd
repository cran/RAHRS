\name{MadgwickAHRSupdateIMU}
\alias{MadgwickAHRSupdateIMU}
\alias{MadgwickAHRSupdateIMU2}
\alias{MadgwickAHRSupdateIMUDbl}
\title{
Madgwick's IMU update method
}
\description{
Madgwick's algorithm for IMU update method.
}
\usage{
MadgwickAHRSupdateIMU(gxi, gyi, gzi, axi, ayi, azi, sampleFreqi, betai, 
q0i, q1i, q2i, sampleFreq)
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
  \item{sampleFreq}{
Sample frequency.
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
