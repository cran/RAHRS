\name{anglesGyroLib}
\alias{anglesGyroLib}
\alias{calibrationGyroLib}
\title{Test recording and calibration recording data}
\usage{anglesGyroLib}
\description{\code{anglesGyroLib} and \code{calibrationGyroLib} contain test recording and calibration recording data read by Vladimir Maximov for the project GyroLib. Data W, A and M should be divided by 3000 to get the actual values.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{sincIMU}{ IMU synch value}
    \item{Wx}{ x-value of gyroscope data, angular rate in rad/sec}
    \item{Wy}{ y-value of gyroscope data, angular rate in rad/sec}
    \item{Wz}{ z-value of gyroscope data, angular rate in rad/sec}
    \item{Ax}{ x-value of accelerometer data, acceleration in g}
    \item{Ay}{ y-value of accelerometer data, acceleration in g}
    \item{Az}{ z-value of accelerometer data, acceleration in g}
    \item{Mx}{ x-value of magnetometer data, magnetic field, uncalibrated, anti-clockwise positive}
    \item{My}{ y-value of magnetometer data, magnetic field, uncalibrated, anti-clockwise positive}
    \item{Mz}{ z-value of magnetometer data, magnetic field, uncalibrated, anti-clockwise positive}
    \item{temperature}{ temperature}
 }
}
\source{
Vlad Maximov, 2012
Scalar Calibration of Vector accelerometers and magnetometers,
GyroLib documentation.
}
\references{
Vlad Maximov, 2012
Scalar Calibration of Vector accelerometers and magnetometers,
GyroLib documentation.
}
\author{Jose Gama}
\examples{
data(anglesGyroLib)
anglesGyroLib
}
\keyword{datasets}
