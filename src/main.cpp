#include <fstream>
#define _USE_MATH_DEFINES

#include <math.h>
#include "uWS/uWS.h"
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"
#include "json.hpp"

#include "spline.h"
#include <map>

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Support globals
std::map<int, double[3]> cars;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
stringstream hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_last_of("]");
  if (found_null != string::npos) {
    return stringstream();
  } 
  else if (b1 != string::npos && b2 != string::npos) {
	  stringstream tmp = stringstream();
	  tmp.str(s.substr(b1, b2 - b1 + 1));
	  return tmp;
  }
  return stringstream();
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(unsigned int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}


/* ============================== Start Support functions ============================== */
// Implementation by Jordi Tudela

// DEFINITION OF MODEL PARAMETERS
// maximum driving speed [m/s]
const double _MAX_SPEED_ = 20;

// Time between updates [s]
const double _DELTA_TIME_ = 0.02;

// Prediction span [s] - how many seconds into the future we are predicting 
const double _TIME_PRED_ = 1;

// Number of data points we will be generating
const int _POINTS_ = (int)(_TIME_PRED_ / _DELTA_TIME_);

// s prediction [s]: we will use a rough predicted s to get smoother paths
const double _FUTURE_S_ = 10;

// Safety distance [m] -	the car will actively try not to get closer than the safe distance to any 
//							car in front, nor will try to change lanes if in the target lane there is a car
//							behind at less than the safe distance.
const double _SAFE_DISTANCE_ = 15;

// Maximum angle variation between steps
const double _MAX_THETA_ = 0.008; // 0.006;

// Index of the most left and most right lanes
const int _MIN_LANE_ = 0;
const int _MAX_LANE_ = 2;

// Maximum acceleration
const double _MAX_ACCEL_ = 5;

// Maximum speed variation within a step
const double _ACCEL_STEP_ = _MAX_ACCEL_ * _DELTA_TIME_;

// List of states
enum Actions { FOLLOW, PREPARE, CHANGE };

// Struct to define persistently the current state
struct State {
	Actions state;
	double target_d;
	unsigned int follow_id;

	State() {
		state = FOLLOW;
		target_d = -1;
		follow_id = 0;
	}
};

// Computes the module of two numbers
double GetModule(double x, double y) { return sqrt(x*x + y*y); }

// Get the lane index of a particular coordinate d
// @ d: frenet coordinate
int GetLane(double d) { return static_cast <int> (std::floor(d / 4.0)); }

// Given an s parameter, find the closes car in that is beyond that parameter
// within the lane that encloses d
// @ s: frenet coordinate s
// @ d: frenet coordinate d
int FindNextCarInLane(double s, double d) {
	int id = -1;
		
	double other_s = std::numeric_limits<double>::infinity();

	for (auto const& x : cars) {

		// Check if the current car is on the same lane
		//if (abs(x.second[1] - car_d) > _LANE_DIST_) continue;
		if (GetLane(x.second[1]) != GetLane(d)) continue;

		// Check if the car is in front of me
		// and Check if it is closer than the previous closest one
		if (s < x.second[0] && x.second[0] < other_s) {
			// update current closest one
			other_s = x.second[0];
			id = x.first;
		}
	}

	return id;
}

// Given the current speed and an objective, update the current speed
// to get close to that objective within the acceleration limitations
// @ speed: current speed >> updated speed 
// @ objective: objective speed
void UpdateSpeed(double& speed, double objective) 
{		
	// we need to control the variation of speed to avoid spikes in acceleration / jerk
	if (speed + _ACCEL_STEP_ < objective) speed += _ACCEL_STEP_;
	else if (speed - _ACCEL_STEP_ > objective) speed -= _ACCEL_STEP_;
	else speed = objective;
}

// Get the next allowed position following a given direction
// @ xy: current position >> next position
// @ direction: direction vector 
// @ speed: speed value to compute the distance between the current position and the next
// @ yaw: current yaw angle of the car
void GetNextStep(vector<double>& xy, vector<double> direction, double speed, double& yaw) 
{
	// get the angle of the direction vector and normalize it [0, 2*PI)
	double theta = atan2(direction[1], direction[0]);
	if (theta < 0) theta += 2 * pi();
	else if (theta >= 2*pi()) theta -= 2 * pi();

	// angle difference between the current yaw and the angle of the direction vector
	// normalized to [-PI, PI]
	double diff = yaw - theta;
	while (diff >   pi()) diff -= 2 * pi();
	while (diff <= -pi()) diff += 2 * pi();

	// Compute the next (x,y) pair with the given speed and angle.
	// If the difference is greater than the allowed angle change, use the maximum
	// angle to compute the next position. If not, use the direction vector's angle.
	if (diff > _MAX_THETA_)
	{
		xy[0] += (speed * _DELTA_TIME_)*cos(yaw - _MAX_THETA_);
		xy[1] += (speed * _DELTA_TIME_)*sin(yaw - _MAX_THETA_);
		yaw -= _MAX_THETA_;
	}
	else if (diff <  -_MAX_THETA_)
	{
		xy[0] += (speed * _DELTA_TIME_)*cos(yaw + _MAX_THETA_);
		xy[1] += (speed * _DELTA_TIME_)*sin(yaw + _MAX_THETA_);
		yaw += _MAX_THETA_;
	}
	else {
		xy[0] += (speed * _DELTA_TIME_)*cos(theta);
		xy[1] += (speed * _DELTA_TIME_)*sin(theta);
		yaw = theta;
	}
}

// Update the next_x and next_y vextors with a new (x,y) pair
// @ xy: new position
// @ index: position within the next_x and next_y vectors
// @ next_x: vector of next x values to feed to the car
// @ next_y: vector of next y values to feed to the car
void UpdateValue(vector<double> xy, unsigned int index, vector<double>& next_x, vector<double>& next_y) 
{
	// Ratio
	// To improve smoothness, we are going to interpolate between the previous path
	// and the path we are currently generating using a cubic formula
	double r = 1 - (0.000016 * index * index * index - 0.0012 * index * index + 0.0008 * index + 1);
	if (r > 1) r = 1;
	else if (r < 0) r = 0;

	// If index does not correspond to a valid position on the next_x and next_y vectors
	// feed the position values as is. Else use a weighted average between the previous
	// value and the new value
	if (next_x.size() <= index)
	{
		next_x.push_back(xy[0]);
		next_y.push_back(xy[1]);
	}
	else
	{
		next_x[index] = (r * xy[0] + (1 - r) * next_x[index]);
		next_y[index] = (r * xy[1] + (1 - r) * next_y[index]);
	}

}

// Drive following the current lane using the following policy:
//	- if there is no car in front: tend to go at max speed
//	- if there is a car in front: the speed is proportional to the gap in front
//								  we try to match speeds but if the gap gets narrower
//								  we reduce speed to widen the gap and avoid undesired collisions
// @ state: current state of the car
// @ car_x: cartesian coordinate x of the car
// @ car_y: cartesian coordinate y of the car
// @ car_s: frenet coordinate s of the car
// @ car_d: frenet coordinate d of the car
// @ car_speed: current speed of the car
// @ yaw: current yaw of the car
// @ next_x: vector of next x values to feed to the car
// @ next_y: vector of next y values to feed to the car
// @ map_waypoints_s: list of s waypoints
// @ map_waypoints_x: list of x waypoints
// @ map_waypoints_y: list of y waypoints
void FollowLane(State& state, double car_x, double car_y, double car_s, double car_d, double car_speed, double yaw, vector<double>& next_x, vector<double>& next_y,
	const vector<double>& map_waypoints_s, const vector<double>& map_waypoints_x, const vector<double>& map_waypoints_y)
{
	// INITIAL VALUES
	// Define the target d corrdinate
	double d = state.target_d;

	// Initally we will set our target speed to the maximum allowed
	double objSpeed = _MAX_SPEED_;

	// The initial speed is the current speed of the car
	double speed = car_speed;

	// The initial angle is the current yaw of the car
	double theta = yaw;

	// The initial (x,y) pair is the current (x,y) position of the car
	vector<double> xy = { car_x, car_y };

	// UPDATE CAR STATE
	// First, we find a car in front of us within the current lane
	// If our state is "follow" the lane, we check the whole car list
	// If our state is "prepare" to change lane, we are currently following a car
	// in our lane, so get its id.
	int id;
	if (state.state == FOLLOW) id = FindNextCarInLane(car_s, d);
	else if (state.state == PREPARE) id = state.follow_id;

	// Find the distance between our car and the one in front
	double diff = cars[id][0] - car_s;
	if (diff <= _SAFE_DISTANCE_)
	{
		// if the car is too close, set the target speed proportionally to the
		// one of the car in front
		objSpeed = (diff / _SAFE_DISTANCE_) * cars[id][2];

		// set the current state of the car as "prepare" to change lanes
		state.state = PREPARE;
		state.follow_id = id;
	}
	// if the car in front is not too close, set the current state as "follow" the lane
	else if (diff > 1.0 * _SAFE_DISTANCE_) state.state = FOLLOW;

	// We define our starting parameter s to be a certain distance in front of the
	// current position of the car. We will be using this parameter to get the
	// direction vector that the car needs to follow and having it in front of us
	// gives us smoother results. 
	// Note that if we place it near the current position, it will create a noisy path
	// but if we place it far away, it won't properly follow lane markings in turns
	double s = car_s + _FUTURE_S_ * speed * _DELTA_TIME_;

	// MAIN LOOP
	for (int i = 0; i < _POINTS_; ++i)
	{
		// First we update the speed at this step to tend to our objective speed within acceleration limitations
		UpdateSpeed(speed, objSpeed);

		// Get the next s position: x(t) = x(0) + v * t
		 s += _DELTA_TIME_ * speed;

		 // Get a new (x,y) rough pair for this step
		vector<double> xy_f = getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

		// Get the direction vector as the difference between the previous (x,y) pair
		// and the one the we have just computed
		vector<double> direction = { (xy_f[0] - xy[0]), (xy_f[1] - xy[1]) };

		// Get the new (x,y) position by apply speed and angle constraints
		GetNextStep(xy,direction,speed,theta);

		// Update the values of the output vectors
		UpdateValue(xy, i, next_x, next_y);

	}
}

// Check side lanes for opportunities to perform a changing lane manouver
void CheckSideLanes(State& state, double car_s, double car_d)
{
	//INITIAL VALUES
	// index of the target lane
	int lane = GetLane(state.target_d);

	// Set the initial values for the cars in front of each lane to -1
	// (i.e. no car found in front)
	int left_id = -1;
	int right_id = -1;

	// Set the default values for frenet s coordinates at -Infinity
	// (i.e. never try to get to that lane)
	double left_s = -std::numeric_limits<double>::infinity();
	double right_s = -std::numeric_limits<double>::infinity();

	// CHECK SIDE LANES
	// if the target lane is valid, we try to find the next car in that lane that
	// would be in front of us. 
	// Then we try to determine if it would be safe to change lanes by setting the s coordinate a 
	// safe distance behind us and trying to find the car in front. If the id's of both cars match,
	// it means that there is no other car in that space in we store its s coordinate. If the id's
	// don't match, there is another car and we keep the -Infinity value.
	if (lane > _MIN_LANE_)
	{
		// Check the left lane (if its a valid lane)
		left_id = FindNextCarInLane(car_s, 2 + 4 * (lane - 1));
		if (FindNextCarInLane(car_s - _SAFE_DISTANCE_, 2 + 4 * (lane - 1)) == left_id) left_s = cars[left_id][0];
	}

	if (lane < _MAX_LANE_)
	{
		// Check the right lane (if its a valid lane)
		right_id = FindNextCarInLane(car_s, 2 + 4 * (lane + 1));
		if (FindNextCarInLane(car_s - _SAFE_DISTANCE_, 2 + 4 * (lane + 1)) == right_id) right_s = cars[right_id][0];
	}

	// CHECK FOR MOST PROFITABLE ALTERNATIVE
	// To commit to change lanes, we want to chose the one with the most free space
	// in front of us, given that this space is enough.
	// In other words, if we had two cars side by side, we wouldn't want to keep changing lanes.
	// *Note that if a lane score set to -Infinity would never be better than the current lane score
	//
	// If we decide to change lanes, we update the car's state and we set the new target d coordinate
	if (right_s > left_s && right_s > car_s + 2 * _SAFE_DISTANCE_)
	{
		// move to right lane
		state.state = CHANGE;
		state.target_d =  2 + 4 * (lane + 1);
	}
	else if (left_s > car_s + 2 * _SAFE_DISTANCE_)
	{
		// move to left lane
		state.state = CHANGE;
		state.target_d = 2 + 4 * (lane - 1);
	}

}

// Get the (x,y) coordinates following a Changing Lane manouver
// To change lanes, match the velocity of the car in front and slightly turn until you reaching the desired d position
// @ state: current state of the car
// @ car_x: cartesian coordinate x of the car
// @ car_y: cartesian coordinate y of the car
// @ car_s: frenet coordinate s of the car
// @ car_d: frenet coordinate d of the car
// @ car_speed: current speed of the car
// @ yaw: current yaw of the car
// @ next_x: vector of next x values to feed to the car
// @ next_y: vector of next y values to feed to the car
// @ map_waypoints_s: list of s waypoints
// @ map_waypoints_x: list of x waypoints
// @ map_waypoints_y: list of y waypoints
void ChangeLane(State& state, double car_x, double car_y, double car_s, double car_d, double car_speed, double yaw, vector<double>& next_x, vector<double>& next_y,
	const vector<double>& map_waypoints_s, const vector<double>& map_waypoints_x, const vector<double>& map_waypoints_y)
{
	// INITIAL VALUES
	// The initial target speed is the speed of the car that is currently in front of us
	double objSpeed = cars[state.follow_id][2];

	// The initial speed is the current speed of the car
	double speed = car_speed;

	// The initial d values are the current d values of the car
	double d = car_d;

	// We define our starting parameter s to be a certain distance in front of the
	// current position of the car. We will be using this parameter to get the
	// direction vector that the car needs to follow and having it in front of us
	// gives us smoother results. 
	// Note that if we place it near the current position, it will create a noisy path
	// but if we place it far away, it won't properly follow lane markings in turns
	double s = car_s + _FUTURE_S_ * speed * _DELTA_TIME_;

	// the initial theta is the current yaw of the car
	double theta = yaw;

	// define the angle and z variables for later use
	double angle;
	double z; 

	// Set the initial (x,y) position to the current (x,y) position of the car
	vector<double> xy = { car_x, car_y};

	// UPDATE STATE
	// Check if there are cars comming in the target lane
	int id = FindNextCarInLane(car_s - 0.5 * _SAFE_DISTANCE_, state.target_d);
	if (car_s - cars[id][0] > 0) 
	{
		if (state.target_d > car_d) state.target_d = 2 + 4 * (GetLane(state.target_d) - 1);
		else state.target_d = 2 + 4 * (GetLane(state.target_d) + 1);
	}
	// if we are resonably close to our target d position, we change states to "follow" lane  
	if (abs(state.target_d - d) < 0.1) state.state = FOLLOW;

	// MAIN LOOP
	for (int i = 0; i < _POINTS_; ++i)
	{
		// First we update the speed at this step to tend to our objective speed within acceleration limitations
		UpdateSpeed(speed, objSpeed);
		
		// we define z as the distance to center of the current lane
		if (GetLane(d) == GetLane(state.target_d)) z = state.target_d - d;
		else z = d - 2 - 4 * GetLane(d);
			
		// We set the current angle in the frenet base using a cuadratic formula
		// This formula has been manually tuned to get a smooth turning and a regular
		// value when tending to 2 (the edge of a lane)
		angle = 0.01 + 0.65 * abs(z) - 0.0525 * z * z;
		
		// We bound that angle to avoid unnecessary wobbling
		if (angle > 0.1) angle = 0.1;
		else if (angle < -0.1) angle = -0.1;

		// Update s coordinate
		s += speed * _DELTA_TIME_ * cos(angle);

		// Update d, depending on its position within the lane
		// if it gets close to the center, just follow the center to avoid wobbling
		if (state.target_d - d > 0.1) d += speed * _DELTA_TIME_ * sin(angle);
		else if (state.target_d - d < -0.1) d -= speed * _DELTA_TIME_ * sin(angle);
		else d = state.target_d;

		// get the new (x,y) rough objective for this step
		auto xy_f = getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

		// get the direction of movement for this step
		vector<double> v = { xy_f[0] - xy[0], xy_f[1] - xy[1] };

		// get the new (x,y) coordinates bound by speed and angle limitations
		GetNextStep(xy, v, speed, theta);

		// Pass the new (x,y) pair to the output vectors
		UpdateValue(xy, i, next_x, next_y);		
	}
}

/* ============================== End Support functions ============================== */


int main() {
  uWS::Hub h;
  State state;

  // Set a default car for id -1 (ie. there is no car)
  // I am placing the "car" at s = infinity. d and speed are irrelevant in this case.
  cars[-1][0] = std::numeric_limits<double>::infinity();
  cars[-1][1] = 0;
  cars[-1][2] = 0;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &state](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(string(data));

      if (s.str() != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
		  // Main car's localization Data
			double car_x = j[1]["x"];
			double car_y = j[1]["y"];
			double car_s = j[1]["s"];
			double car_d = j[1]["d"];
			double car_yaw = j[1]["yaw"];
			double car_speed = 0.44704 * (double)j[1]["speed"];

			// Previous path data given to the Planner
			auto previous_path_x = j[1]["previous_path_x"];
			auto previous_path_y = j[1]["previous_path_y"];
			// Previous path's end s and d values 
			double end_path_s = j[1]["end_path_s"];
			double end_path_d = j[1]["end_path_d"];

			// Sensor Fusion Data, a list of all other cars on the same side of the road.
			auto sensor_fusion = j[1]["sensor_fusion"];


          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	// Define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
			// Implementation by Jordi Tudela
			// ==========================================================================================
			try {
				// Initialize the state's target_d coordinate with the center position of the current lane
				if (state.target_d == -1) state.target_d = 2 + 4 * GetLane(car_d);

				// I pass all the previous path points to the new path vectors becuse we are
				// doing some averaging when computing the path in this iteration to try to get
				// a smoother driving style
				for (unsigned int i = 0; i < previous_path_x.size(); i++)
				{
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
				}

				// Update the dictionary/map of other cars with their updated
				// frenet coordinates and absolute speed
				for (unsigned int i = 0; i < sensor_fusion.size(); i++)
				{
					int id = sensor_fusion[i][0];

					cars[id][0] = sensor_fusion[i][5]; // Frenet Coordinate - s
					cars[id][1] = sensor_fusion[i][6]; // Frenet Coordinate - d
					cars[id][2] = GetModule(sensor_fusion[i][3], sensor_fusion[i][4]); // Total Velocity - v
				}

				// I have defined to movement functions: FollowLane and ChangeLane
				// FollowLane will tend to go at the maximum speed possible within the current lane and when it find a car in front,
				//			  the car will slow down to avoid collision (it will match speeds at a safe distance and it will slow
				//			  down more it the distance is not safe to create a wider space between them)
				// ChangeLane will match speeds with the car in front and it will manouver to change lane
				if (state.state == CHANGE)
					ChangeLane(state, car_x, car_y, car_s, car_d, car_speed, deg2rad(car_yaw), next_x_vals, next_y_vals, map_waypoints_s, map_waypoints_x, map_waypoints_y);
				else
					FollowLane(state, car_x, car_y, car_s, car_d, car_speed, deg2rad(car_yaw), next_x_vals, next_y_vals, map_waypoints_s, map_waypoints_x, map_waypoints_y);

				// The "Prepare" state triggers when the cars fins another car in front and it triggers CheckSideLanes
				// CheckSideLanes looks the immediate viable side lanes of the lane the car is currently drivind to 
				// see if it would be possible to change lanes.
				if (state.state == PREPARE) CheckSideLanes(state, car_s, car_d);
			}
			catch (const std::exception&) {
				std::cout << "ERROR: There has been an error in the normal flow of the algorithm. Using previouly computed path." << std::endl;
 				for (unsigned int i = 0; i < previous_path_x.size(); ++i) {
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
				}
			}
			// ==========================================================================================


          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen("0.0.0.0", port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}