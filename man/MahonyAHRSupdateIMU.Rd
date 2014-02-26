\name{MahonyAHRSupdateIMU}
\alias{MahonyAHRSupdateIMU}
\alias{MahonyAHRSupdateIMU2}
\alias{MahonyAHRSupdateIMUDbl}
\title{
Mahony's IMU update method
}
\description{
Mahony's algorithm for IMU update method.
}
\usage{
MahonyAHRSupdateIMU(gxi, gyi, gzi, axi, ayi, azi, sampleFreqi, twoKpi, 
twoKii, integralFBxi, integralFByi, integralFBzi, q0i, q1i, q2i, halfex)
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
  \item{twoKpi}{
Kp constant passed as a parameter
}
  \item{twoKii}{
Ki constant passed as a parameter
}
  \item{integralFBxi}{
Integral FB[0]
}
  \item{integralFByi}{
Integral FB[1]
}
  \item{integralFBzi}{
Integral FB[2]
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
  \item{halfex}{
Half 
}
}
\value{
output quaternion describing the Earth relative to the sensor.
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
