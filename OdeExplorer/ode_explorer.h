#pragma once

#include "stdafx.h"

using namespace std;
using namespace cimg_library;
using namespace boost::numeric::odeint;




//==========================//
// ode_explorer Class Header //
//==========================//

template <typename num_type, typename state_type>
class ode_explorer {




	//==================//
	// Type Definitions //
	//==================//

	private:

	// matrix type used for storing data
	typedef vector< vector< num_type > > matrix_type;




	//===========//
	// Variables //
	//===========//

	public:
	
	// integration tolerance
	num_type tolerance;

	// initial conditions
	state_type initial_state;
	
	// image parameters
	int width = 0;
	int height = 0;
	int num_images = 0;

	// number of events to track
	int num_events = 0;




	//==============//
	// Constructors //
	//==============//

	public:

	ode_explorer(num_type tolerance);




	//========================//
	// Pure Virtual Functions //
	//========================//

	private:

	// ode function to integrate
	virtual void ode(const state_type x, state_type & dxdt, const num_type t) = 0;

	// set parameters depending on the pixel and image number
	virtual void set_pixel_params(int x, int y, int n) = 0;

	// determines event to track
	virtual num_type event(const state_type x, const num_type t) = 0;

	


	//===========//
	// Functions //
	//===========//

	public:

	// calculates the state and time of the specified number of events
	void find_event(state_type &x, num_type &t);

	// calculates the normal between an initial state and final state
	num_type calculate_event_norm(const state_type x0);
	
	// generates a matrix of calculated normals
	vector<vector<num_type>> generate_composite_data(int width, int height);

	CImg<int> generate_composite_image(int width, int height);

	vector<CImg<int>> generate_composite_animation(int width, int height, int num_frames);

};
