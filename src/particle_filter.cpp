/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#include "helper_functions.h"
#include "map.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).



	num_particles=10;
	default_random_engine gen;

	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);
	
	for(int i=0; i< num_particles; i++)
	{
		Particle p;
		p.id=i;
		p.x= dist_x(gen);  
		p.y= dist_y(gen);
		p.theta=  dist_theta(gen);
		p.weight= 1.0;
		particles.push_back(p);
		weights.push_back(1.0);
	}

	is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	normal_distribution<double> dist_pos_x(0, std_pos[0]);
	normal_distribution<double> dist_pos_y(0, std_pos[1]);
	normal_distribution<double> dist_pos_theta(0, std_pos[2]);	
	default_random_engine gen;

	for(int i=0; i< num_particles; i++)
	{

		double temp_yaw_rate=yaw_rate;
		if(fabs(yaw_rate)<0.00001)
		{
			temp_yaw_rate=0.00001;
		} 
		double theta= particles[i].theta;
		particles[i].x += velocity / temp_yaw_rate * (sin(theta + temp_yaw_rate*delta_t) - sin(theta)) + dist_pos_x(gen);
		particles[i].y += velocity / temp_yaw_rate * (cos(theta) - cos(theta + temp_yaw_rate*delta_t)) + dist_pos_y(gen);
		particles[i].theta += temp_yaw_rate * delta_t + dist_pos_theta(gen);

		// particles[i].x+= dist_pos_x(gen);
		// particles[i].y+= dist_pos_y(gen);
		// particles[i].theta+= dist_pos_theta(gen);

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	
	for(int i=0; i<observations.size(); i++)
	{
		double o_x= observations[i].x;
		double o_y= observations[i].y;
		int min_id=-1;
		double min_dist=9999999.9;

		for(int j=0; j<predicted.size(); j++)
		{
			double p_x= predicted[j].x;
			double p_y= predicted[j].y;
			double distance= dist(p_x, p_y, o_x, o_y);
			if(distance<min_dist)
			{
				min_dist=distance;
				min_id= predicted[j].id;
			}

		}
		observations[i].id= min_id;
		// std::cout<<min_dist<<" "<<std::endl;
	}




}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// std::cout<<"num: "<<num_particles<<std::endl;

	double s_x = std_landmark[0];
	double s_y = std_landmark[1];
	double base= ( 1/(2*M_PI*s_x*s_y));
	weights.clear();
	for(int i=0; i<num_particles; i++)
	{
		
		double p_x= particles[i].x;
		double p_y= particles[i].y;
		double p_theta= particles[i].theta;

		vector<LandmarkObs> predictions;


		for(int j=0; j<map_landmarks.landmark_list.size(); j++)
		{
			auto landmark_list= map_landmarks.landmark_list[j];
			int l_id= landmark_list.id_i;
			double l_x= landmark_list.x_f;
			double l_y= landmark_list.y_f;

			double distance= dist(p_x, p_y, l_x, l_y);
			if(distance<=sensor_range)
				predictions.push_back(LandmarkObs{l_id, l_x,l_y});
		}

		vector<LandmarkObs> transformed_observations;

		for(int j=0; j<observations.size(); j++)
		{
			double cos_theta=cos(p_theta);
			double sin_theta=sin(p_theta);
			LandmarkObs obs= observations[j];
			double to_x= p_x+ cos_theta* obs.x- sin_theta* obs.y;
			double to_y= p_y+ sin_theta* obs.x+ cos_theta* obs.y;
			// obs.id=j+num_particles;
			// std::cout<<p_x<<" "<<p_y<<std::endl;
			transformed_observations.push_back(LandmarkObs{ obs.id, to_x, to_y });
		}

		dataAssociation(predictions,transformed_observations);
		// particles[i].weight = 1.0;
		// vector<double> sense_x;
		// vector<double> sense_y;
		// vector<int> associations;
		double particle_weight=1;

		for(int j=0; j<transformed_observations.size(); j++)
		{
			LandmarkObs t_obs=transformed_observations[j];
			double to_x1= t_obs.x;
			double to_y1= t_obs.y;
			int to_id= t_obs.id; // this also corresponds to the id in predicted landmarks.
			double pr_x, pr_y;
			// sense_x.push_back(to_x1);
			// sense_y.push_back(to_y1);
			// associations.push_back(to_id);

			for (int k = 0; k < predictions.size(); k++) {
				if (predictions[k].id == to_id) {
					pr_x = predictions[k].x;
					pr_y = predictions[k].y;
					break;
				}
			}

			double obs_w =  base* exp( -( pow(pr_x-to_x1,2)/(2*pow(s_x, 2)) + (pow(pr_y-to_y1,2)/(2*pow(s_y, 2))) ) );
			particle_weight *= obs_w;
			// weights[i]=particles[i].weight;
		}
		particles[i].weight= particle_weight;
		weights.push_back(particle_weight);
		// SetAssociations(particles[i],particles[i].associations,particles[i].sense_x,particles[i].sense_y);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<Particle> new_particles;
	default_random_engine gen;

	uniform_int_distribution<int> uniform_int_dist(0, num_particles-1);
	auto index = uniform_int_dist(gen);

	// average weight
	double avg_of_elems = accumulate(weights.begin(), weights.end(), 0.0)/num_particles;

	// uniform real distribution
	uniform_real_distribution<double> uniform_real_dist(0.0, avg_of_elems);

	double beta = 0.0;

	// resample wheel
	for (int i = 0; i < num_particles; i++) {
		beta += uniform_real_dist(gen) * 2.0;
		while (beta > weights[index]) {
		beta -= weights[index];
		index = (index + 1) % num_particles;
		}
		new_particles.push_back(particles[index]);
	}
	particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
