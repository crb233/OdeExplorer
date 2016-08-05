#pragma once

#include "stdafx.h"
#include "ode_explorer.h"

using namespace std;
using namespace cimg_library;
using namespace boost::numeric::odeint;




//==============//
// Constructors //
//==============//

template <typename num_type, typename state_type>
ode_explorer<num_type, state_type>::ode_explorer(num_type tolerance) {
	this.tolerance = tolerance;
}




//===========//
// Functions //
//===========//

// calculates the state and time of the specified number of events
template <typename num_type, typename state_type>
void ode_explorer<num_type, state_type>::find_event(state_type &x, num_type &t) {

	// stepper type used for integration
	typedef runge_kutta_dopri5< state_type > stepper_type;
	
	// determine the initial dt time step
	state_type x_temp = state_type(x);
	num_type t_temp = num_type(t);
	num_type dt = 0.1;
	make_controlled(tolerance, tolerance, stepper_type()).try_step(ode, x_temp, t_temp, dt);

	// initialize vars
	state_type x_prev;
	num_type t_prev;
	int count = num_events;

	// create the stepper
	controlled_runge_kutta< stepper_type > stepper;
	stepper = make_controlled(tolerance, tolerance, stepper_type());

	// take an initial step
	stepper.try_step(ode, x, t, dt);
	t += dt;

	// integrate up to the last event
	while (count > 0)) {
		x_prev = x;
		t_prev = t;
		stepper.try_step(ode, x, t, dt);
		t += dt;

		if (event(x_prev, t_prev) * event(x, t) <= 0) {
			count_a--;
		}
	}

	num_type min_dt = 0;
	num_type max_dt = dt;

	// interpolate to higher precision solution
	while (abs(event(x, t)) > tolerance) {

		x = x_prev;
		t = t_prev;
		dt = (min_dt + max_dt) / 2;
		stepper_type().do_step(ode, x, t, dt);

		if (event(x_prev, t_prev) * event(x, t) <= 0) {
			max_dt = dt;
		} else {
			min_dt = dt;
		}
	}

}

// calculates the normal between an initial state and final state
template <typename num_type, typename state_type>
num_type ode_explorer<typename num_type, typename state_type>::
calculate_event_norm(const state_type x0) {

}

// generates a matrix of calculated normals
template <typename num_type, typename state_type>
vector<vector<num_type>> ode_explorer<typename num_type, typename state_type>::
generate_composite_data(int width, int height) {

}

template <typename num_type, typename state_type>
CImg<int> ode_explorer<typename num_type, typename state_type>::
generate_composite_image(int width, int height) {

}

template <typename num_type, typename state_type>
vector<CImg<int>> ode_explorer<typename num_type, typename state_type>::
generate_composite_animation(int width, int height, int num_frames) {

}
