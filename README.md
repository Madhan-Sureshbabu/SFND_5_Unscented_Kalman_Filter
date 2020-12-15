# Udacity Nanodegree : Sensor Fusion 

## Unscented Kalman Filter 

<img src="media/ukf_highway_tracked.gif" width="700" height="400" />

The objective of this project is to implement an Unscented Kalman Filter to estimate the state of multiple cars on a highway using noisy lidar and radar measurements. 

The red spheres above cars represent the (x,y) lidar detection and the purple lines show the radar measurements with the velocity magnitude along the detected angle. The Z axis is not taken into account for tracking, so you are only tracking along the X/Y axis.

```
```

1. CTRV (Constant Turn Rate and Velocity) motion model has been used for predicting the state of the vehicle. This prediction step uses the time delta_t from the last update to account for uncertainity. Sigma points are generated from the state mean and covariance matrix augmented with the process noise parameters. These are propagated through the prediction step to obtain the predicted sigma points. 
2. The measurement update step with the radar uses the predicted sigma points to propagate through the non-linear measurement function before updating the state.
3. The measurement update step with the lidar uses the predicted mean and covariance matrix calculated from the predicted sigma points for updating the state.

Project pipeline :  
<img src="media/block_diagram.png" width="700" height="400" />


### Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
 * PCL 1.2

### Basic Build Instructions

1. Clone this repo.  
2. Create a build directory:  
```
mkdir build && cd build
```
3. Compile:  
```
cmake .. && make
```
4. Run it:  
```
./ukf_highway
```


