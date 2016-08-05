#pragma once

#include "stdafx.h"
#include "OdeExplorer.h"

using namespace std;
using namespace cimg_library;
using namespace boost::numeric::odeint;

template <typename num_type>
template <typename state_type>

//==================//
// Type Definitions //
//==================//

// matrix type used for storing data
//typedef vector< vector< num_type > > matrix_type;

// stepper type used for integration
//typedef runge_kutta_dopri5< state_type > stepper_type;




//==============//
// Constructors //
//==============//

template <typename num_type>
template <typename state_type>
OdeExplorer::OdeExplorer() {

}




//===========//
// Functions //
//===========//

// calculates the state and time of the specified number of events
template <typename num_type>
template <typename state_type>
void OdeExplorer::findEvent(state_type &x, num_type &t, event_type type) {

}

// calculates the normal between an initial state and final state
template <typename num_type>
template <typename state_type>
num_type OdeExplorer::calculateEventNorm(const state_type x0) {

}

// generates a matrix of calculated normals
template <typename num_type>
template <typename state_type>
matrix_type OdeExplorer::generate_composite_data(int width, int height) {

}

template <typename num_type>
template <typename state_type>
CImg OdeExplorer::generate_composite_image(int width, int height) {

}

template <typename num_type>
template <typename state_type>
vector<CImg> OdeExplorer::generate_composite_animation(int width, int height, int num_frames) {

}
