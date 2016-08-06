
#include "stdafx.h"

using namespace std;
using namespace cimg_library;
using namespace boost::numeric::odeint;




//==================//
// Type Definitions //
//==================//

typedef double num_type;
typedef boost::array< num_type, 4 > state_type;
typedef vector< vector< num_type > > matrix_type;
typedef runge_kutta_dopri5< state_type > stepper_type;

enum mode {mode_a, mode_b, mode_a_and_b};
enum event_type {type_a, type_b, type_a_and_b, type_a_or_b};




//=====================//
// Function Prototypes //
//=====================//

// calculates the state and time of the specified number of events
void findEvent(state_type &x, num_type &t, event_type type);

// calculates the normal between an initial state and final state
num_type calculateEventNorm(const state_type x0);

// generates a matrix of calculated normals
matrix_type generateData();

// saves data as a csv file
bool saveCSV(string filename, const matrix_type data);

// saves data as a bmp file
bool saveBMP(string filename, const matrix_type data);




//===========//
// Variables //
//===========//

// ode function to integrate
void (*ode) (const state_type x, state_type & dxdt, const num_type t);

// event to track
num_type (*event_a) (const state_type x, const num_type t);
num_type (*event_b) (const state_type x, const num_type t);

// integration tolerance
num_type tolerance;

// number of events to track
int num_events_a;
int num_events_b;

// which events to track
mode event_mode;

// initial conditions
num_type x_pos = 0;
num_type y_pos = 0;

// image parameters
int width;
int height;
num_type xv0;
num_type yv0;
num_type dv;




//==============================//
// Main Function and Parameters //
//==============================//

// object mass
num_type m = 1.0;

// spring constant
num_type k1 = 1.00;
num_type k2 = 1.00;

// spring natural length
num_type L1 = 2.0;
num_type L2 = 2.0;

// position of spring on x-axis
num_type p1 = -2.0;
num_type p2 = 2.0;

// degree of relation between force and position F=k*(d-L)^n
num_type deg = 0.0;

// variable degree force toward points ode
void degree_force_ode(const state_type x, state_type &dxdt, const num_type t) {
	num_type c1 = sqrt(pow(p1 - x[0], 2) + pow(x[1], 2)) - L1;
	c1 = ((c1 > 0) - (c1 < 0)) * pow(abs(c1), deg);
	num_type c2 = sqrt(pow(p2 - x[0], 2) + pow(x[1], 2)) - L2;
	c2 = ((c2 > 0) - (c2 < 0)) * pow(abs(c2), deg);

	dxdt[0] = x[2];
	dxdt[1] = x[3];
	dxdt[2] = (k1 * c1 * (p1 - x[0]) + k2 * c2 * (p2 - x[0])) / m;
	dxdt[3] = -(k1 * c1 + k2 * c2) * x[1] / m;
}

// double spring system ode
void springs_ode(const state_type x, state_type &dxdt, const num_type t) {
	num_type c1 = 1 - L1 / sqrt(pow(p1 - x[0], 2) + pow(x[1], 2));
	num_type c2 = 1 - L2 / sqrt(pow(p2 - x[0], 2) + pow(x[1], 2));

	dxdt[0] = x[2];
	dxdt[1] = x[3];
	dxdt[2] = (k2 * c2 * (p2 - x[0]) + k1 * c1 * (p1 - x[0])) / m;
	dxdt[3] = -(k2 * c2 + k1 * c1) * x[1] / m;
}

// constant force toward points ode
void constant_force_ode(const state_type x, state_type &dxdt, const num_type t) {
	num_type c1 = 1 / sqrt(pow(p1 - x[0], 2) + pow(x[1], 2));
	num_type c2 = 1 / sqrt(pow(p2 - x[0], 2) + pow(x[1], 2));
	
	dxdt[0] = x[2];
	dxdt[1] = x[3];
	dxdt[2] = (k2 * c2 * (p2 - x[0]) + k1 * c1 * (p1 - x[0])) / m;
	dxdt[3] = -(k2 * c2 + k1 * c1) * x[1] / m;
}

int main() {
	ode = degree_force_ode;
	event_a = [](const state_type x, const num_type t) {
		return x[1];
	};
	event_b = [](const state_type x, const num_type t) {
		return x[0];
	};



	// number of X axis crossings
	num_events_a = 6;

	// number of Y axis crossings
	num_events_b = 0;

	// which events to track
	event_mode = mode::mode_a;

	// integration tolerance
	tolerance = 1e-10;

	// image parameters
	width = 400;
	height = 200;
	xv0 = -4;
	yv0 = 0;
	dv = 0.02;

	// directory
	string directory = "C:\\Users\\Curtis\\Documents\\C++\\OdeExplorer\\Data\\";

	// name of file to save to
	string base_name = "img";

	// number of images
	int num_imgs = 21;

	// change the attributes for each image
	void(*func)(int i) = [](int i) {
		deg += 0.1;
	};



	// generate the images
	
	int frame = 0;
	int group = 0;
	base_name = directory + base_name;

	while (true) {
		string f = base_name + "_" + to_string(group) + "_" + to_string(frame);
		
		if (!ifstream(f + ".bmp") && !ifstream(f + ".png"))
			break;
		group++;
	}
	
	base_name = base_name + "_" + to_string(group) + "_";
	for (int i = 0; i < num_imgs; i++) {
		func(i);

		// generate a matrix of data
		matrix_type data = generateData();

		// try saving to a new bmp file
		saveBMP(base_name + to_string(frame), data);
		frame++;

		cout << endl;
	}

}




//===========================//
// General Purpose Functions //
//===========================//

// calculates the state and time of the specified number of events
void findEvent(state_type &x, num_type &t, event_type type) {

	// determine the initial dt time step
	state_type x_temp = state_type(x);
	num_type t_temp = num_type(t);
	num_type dt = 0.1;
	make_controlled(tolerance, tolerance, stepper_type())
		.try_step(ode, x_temp, t_temp, dt);

	// initialize vars
	state_type x_prev;
	num_type t_prev;
	int count_a = num_events_a;
	int count_b = num_events_b;

	// operation to use
	bool(*operation) (bool lhs, bool rhs);
	switch (type) {
	case event_type::type_a:
		operation = [](bool lhs, bool rhs) { return lhs; };
		break;
	case event_type::type_b:
		operation = [](bool lhs, bool rhs) { return rhs; };
		break;
	case event_type::type_a_and_b:
		operation = [](bool lhs, bool rhs) { return lhs || rhs; };
		break;
	case event_type::type_a_or_b:
		operation = [](bool lhs, bool rhs) { return lhs && rhs; };
		break;
	default:
		operation = [](bool lhs, bool rhs) { return lhs; };
		break;
	}

	// create the stepper
	controlled_runge_kutta< stepper_type > stepper;
	stepper = make_controlled(tolerance, tolerance, stepper_type());

	// take an initial step
	stepper.try_step(ode, x, t, dt);
	t += dt;

	int final = 0;
	
	// integrate up to the last event
	while (operation(count_a > 0, count_b > 0)) {
		x_prev = x;
		t_prev = t;
		stepper.try_step(ode, x, t, dt);
		t += dt;

		final = 0;
		if (count_a > 0 && event_a(x_prev, t_prev) * event_a(x, t) <= 0) {
			count_a--;
			final++;
		}
		if (count_b > 0 && event_b(x_prev, t_prev) * event_b(x, t) <= 0) {
			count_b--;
			final++;
		}
	}

	// determine the final event to track
	num_type(*event) (const state_type x, const num_type t);

	if (type == event_type::type_a) {
		event = event_a;
	} else if (type == event_type::type_b) {
		event = event_b;
	} else if (event_a(x_prev, t_prev) * event_a(x, t) <= 0) {
		//
		// a bit problematic code here:
		// causes bugs when both events occur close together
		//
		if (event_b(x_prev, t_prev) * event_b(x, t) <= 0) {
			event = event_b;
		} else {
			event = event_a;
		}
	} else {
		event = event_b;
	}

	num_type min_dt = 0;
	num_type max_dt = dt;

	// interpolate to higher precision solution
	while (abs(event(x, t)) > tolerance) {

		x = x_prev;
		t = t_prev;
		dt = (min_dt + max_dt) / 2;
		stepper_type().do_step(ode, x, t, dt);

		if (event(x, t) * event(x_prev, t_prev) <= 0) {
			max_dt = dt;
		} else {
			min_dt = dt;
		}
	}

}

// calculates the squared normal between the initial state and the calculated event state
num_type calculateEventNorm(const state_type x0) {
	
	// copy the initial state and time
	num_type t = 0;
	state_type x = x0;
	num_type sum = 0;

	if (event_mode == mode::mode_a) {
		findEvent(x, t, event_type::type_a);
		for (int i = 0; i < x0.size(); i++) {
			sum += pow(x0[i] - x[i], 2);
		}

	} else if (event_mode == mode::mode_b) {
		findEvent(x, t, event_type::type_b);
		for (int i = 0; i < x0.size(); i++) {
			sum += pow(x0[i] - x[i], 2);
		}

	} else {
		findEvent(x, t, event_type::type_a_and_b);
		for (int i = 0; i < x0.size(); i++) {
			sum += pow(x0[i] - x[i], 2);
		}
		t = 0;
		x = x0;
		findEvent(x, t, event_type::type_a_or_b);
		for (int i = 0; i < x0.size(); i++) {
			sum += pow(x0[i] - x[i], 2);
		}
		sum /= 2;
	}
	
	return sum;
}

// generates a matrix of calculated normals
matrix_type generateData() {
	// start the timer
	int runtime = time(NULL);

	// initialize data matrix
	matrix_type data;

	// main calculation loop
	int p = -1;
	for (int y = 0; y < height; y++) {

		// output the percent completion
		if (y * 10 / height > p) {
			p = y * 10 / height;
			cout << (10 * p) << "% ";
		}

		// generate another row
		data.push_back(vector<num_type>());
		for (int x = 0; x < width; x++) {
			data[y].push_back(calculateEventNorm({ x_pos, y_pos, xv0 + dv * x, yv0 + dv * y }));
		}
	}
	cout << "100%" << endl;

	// output total runtime
	runtime = time(NULL) - runtime;
	cout << "Total runtime = " << (runtime / 60) << ":";
	cout << setfill('0') << setw(2) << (runtime % 60) << endl;

	return data;
}

bool saveData(string filename) {

}

// saves data as a csv file
bool saveCSV(string filename, const matrix_type data) {
	if (ifstream(filename + ".csv")) {
		return false;
	}

	ofstream out(filename + ".csv", ofstream::out);
	if (!out.good()) {
		out.close();
		return false;
	}

	// output each row of data
	for (int x = 0; x < data.size(); x++) {
		out << data[x][0];
		for (int y = 1; y < data[x].size(); y++) {
			out << "," << data[x][y];
		}
		out << endl;
	}
	out.close();

	cout << filename << ".csv" << endl;
	return true;
}

// saves data as a bmp file
bool saveBMP(string filename, const matrix_type data) {
	if (ifstream(filename + ".bmp")) {
		return false;
	}
	if (ifstream(filename + ".png")) {
		return false;
	}

	try {
		unsigned int width = data[0].size();
		unsigned int height = data.size();

		CImg<num_type> img(width, height);

		for (unsigned int x = 0; x < width; x++) {
			for (unsigned int y = 0; y < height; y++) {
				img(x, y) = data[y][x];
			}
		}

		// edit and save log image
		img.mirror('y');
		img += tolerance;
		img.log();
		img -= log(tolerance);
		img /= img.max();
		img *= 255;
		img.save_bmp((filename + ".bmp").c_str());
		cout << filename << ".bmp" << endl;

		//CImgDisplay disp;
		//disp.display(img_log);
		//disp.show();
		//while (!disp.is_closed());

		return true;

	}
	catch (...) {
		return false;
	}
}
