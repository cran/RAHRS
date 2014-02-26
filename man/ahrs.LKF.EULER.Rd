\name{ahrs.LKF.EULER}
\alias{ahrs.LKF.EULER}
\title{
Euler angles estimation with vector matching and Kalman filter
}
\description{
Attitude Euler angles estimation by means of complementary
Kalman filter.
}
\usage{
ahrs.LKF.EULER(Sensors, State, Parameters)
}
\arguments{
  \item{Sensors}{ calibrated gyroscope, accelerometer and magnetometer 
  	  measurements
    Sensors.w    current calibrated gyroscope measurement [3x1], rad/sec
    Sensors.a    current calibrated accelerometer measurement [3x1], g
    Sensors.m    current calibrated magnetometer measurement [3x1], |m| = 1
  }

  \item{State}{ previous state
    State.q
    State.dB
    State.dG
    State.dw
    State.P
}
  \item{Parameters}{ AHRS Parameters
    Parameters.mn      Magnetic Field Vector In Navigation Frame [3x1], |m|= 1
    Parameters.an      Acceleration vector In Navigation Frame [3x1], g
    Parameters.dt      Sampling period, 1/Hz
  }
}
\value{
   Attitude    estimated attitude Euler angles [1x4]
   State    estimated current state
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
