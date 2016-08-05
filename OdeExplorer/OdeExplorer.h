#pragma once

#include "stdafx.h"

using namespace std;
using namespace cimg_library;
using namespace boost::numeric::odeint;

template <typename num_type>
template <typename state_type>
class OdeExplorer {

public:
	
	//==================//
	// Type Definitions //
	//==================//

	// matrix type used for storing data
	typedef vector< vector< num_type > > matrix_type;

	// stepper type used for integration
	typedef runge_kutta_dopri5< state_type > stepper_type;




	//===========//
	// Variables //
	//===========//

	// integration tolerance
	num_type tolerance;

	// initial conditions
	state_type initial_state;
	
	// image parameters
	int width = 0;
	int height = 0;

	// number of events to track
	int num_events = 0;

	// functions
	template <typename num_type,  state_type>
	num_type (*event) (const state_type x, const num_type t);

	template <typename num_type, state_type>
	state_type (*pixel_to_state) (int x, int y);




	//==============//
	// Constructors //
	//==============//

	public OdeExplorer();




	//========================//
	// Pure Virtual Functions //
	//========================//

	// ode function to integrate
	private virtual void ode(const state_type x, state_type & dxdt, const num_type t) = 0;

	// determines event to track
	// private virtual num_type event(const state_type x, const num_type t) = 0;

	// generates initial integration state
	// private virtual state_type pixel_to_state(int x, int y) = 0;

	//




	//===========//
	// Functions //
	//===========//

	// calculates the state and time of the specified number of events
	private void findEvent(state_type &x, num_type &t, event_type type);

	// calculates the normal between an initial state and final state
	private num_type calculateEventNorm(const state_type x0);

	// generates a matrix of calculated normals
	private matrix_type generate_composite_data(int width, int height);

	public CImg generate_composite_image(int width, int height);

	public vector<CImg> generate_composite_animation(int width, int height, int num_frames);

};
