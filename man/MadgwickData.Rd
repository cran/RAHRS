\name{MadgwickData}
\alias{MadgwickData}
\title{Test data for Sebastian Madgwick's IMU and AHRS sensor fusion algorithms}
\usage{MadgwickData}
\description{\code{MadgwickData} contains test data created by Sebastian Madgwick for his IMU and AHRS sensor fusion algorithms.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Time}{ Time of reading}
    \item{Gx}{ x-value of gyroscope data}
    \item{Gy}{ y-value of gyroscope data}
    \item{Gz}{ z-value of gyroscope data}
    \item{Ax}{ x-value of accelerometer data}
    \item{Ay}{ y-value of accelerometer data}
    \item{Az}{ z-value of accelerometer data}
    \item{Mx}{ x-value of magnetometer data}
    \item{My}{ y-value of magnetometer data}
    \item{Mz}{ z-value of magnetometer data}
 }
}
\source{
S. O. H. Madgwick, An efficient orientation filter for inertial and 
inertial/magnetic sensor arrays, Technical report, 
University of. Bristol University, UK, 2010.
}
\references{
S. O. H. Madgwick, An efficient orientation filter for inertial and 
inertial/magnetic sensor arrays, Technical report, 
University of. Bristol University, UK, 2010.
}
\author{Jose Gama}
\examples{
data(MadgwickData)
MadgwickData
}
\keyword{datasets}
