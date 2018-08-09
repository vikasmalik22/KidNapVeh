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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 12;
	double std_x, std_y, std_theta;

	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	default_random_engine gen;

	for (int i = 0; i < num_particles; i++)
	{
		Particle P;		
		// TODO: Sample  and from these normal distrubtions like this: 
		//	 sample_x = dist_x(gen);
		//	 where "gen" is the random engine initialized earlier.
		P.id = i;
		P.x = dist_x(gen);
		P.y = dist_y(gen);
		P.theta = dist_theta(gen);
		P.weight = 1.0;

		// add p to your particles vector
		particles.push_back(P);
		weights.push_back(P.weight);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	double std_x, std_y, std_theta;
	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];
	default_random_engine gen;

	for (int i = 0; i < num_particles; i++)
	{
		double particle_x = particles[i].x;
		double particle_y = particles[i].y;
		double particle_theta = particles[i].theta;

		double pred_x;
		double pred_y;
		double pred_theta;

		//avoid division by zero
		if (fabs(yaw_rate) > 0.0001) {
			pred_x = particle_x + (velocity / yaw_rate) * (sin(particle_theta + yaw_rate*delta_t) - sin(particle_theta));
			pred_y = particle_y + (velocity / yaw_rate) * (cos(particle_theta) - cos(particle_theta + yaw_rate*delta_t));
			pred_theta = particle_theta + (yaw_rate * delta_t);
		}
		else {
			pred_x = particle_x + (velocity*delta_t*cos(particle_theta));
			pred_y = particle_y + (velocity*delta_t*sin(particle_theta));
			pred_theta = particle_theta;
		}

		//add noise
		// This line creates a normal (Gaussian) distribution for x
		normal_distribution<double> dist_x(pred_x, std_x);
		normal_distribution<double> dist_y(pred_y, std_y);
		normal_distribution<double> dist_theta(pred_theta, std_theta);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i = 0; i < observations.size(); i++)
	{
		LandmarkObs obs = observations[i];

		double closest = numeric_limits<double>::max();

		// init id of landmark from map placeholder to be associated with the observation
		int map_id = -1;

		for (int j = 0; j < predicted.size(); j++)
		{
			LandmarkObs pred = predicted[j];

			double distance = dist(obs.x, obs.y, pred.x, pred.y);

			if (distance < closest)
			{
				closest = distance;
				map_id = pred.id;
			}
		}

		observations[i].id = map_id;
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

	for (int i = 0; i < num_particles; i++)
	{

		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		//// create a vector to hold the map landmark locations predicted to be within sensor range of the particle
		std::vector<LandmarkObs> LandmarkObs_InRange;

		for (int l = 0; l < map_landmarks.landmark_list.size(); l++)
		{
			float landmark_posx = map_landmarks.landmark_list[l].x_f;
			float landmark_posy = map_landmarks.landmark_list[l].y_f;
			int landmark_id = map_landmarks.landmark_list[l].id_i;

			double distance = dist(landmark_posx, landmark_posy, p_x, p_y);

			//if (fabs(landmark_posx - p_x) <= sensor_range && fabs(landmark_posy - p_y) <= sensor_range)
			if (distance <= sensor_range)
			{
				LandmarkObs_InRange.push_back(LandmarkObs{ landmark_id, landmark_posx, landmark_posy });
			}
		}

		std::vector<LandmarkObs> transf_obs;
		//Transform Observations into map coordinates from car coordinates
		for (int j = 0; j < observations.size(); j++)
		{
			double m_x = p_x + (cos(p_theta)*observations[j].x) - (sin(p_theta)*observations[j].y);
			double m_y = p_y + (sin(p_theta)*observations[j].x) + (cos(p_theta)*observations[j].y);
			transf_obs.push_back(LandmarkObs{ observations[j].id, m_x, m_y});
		}
		
		//DataAssociation
		dataAssociation(LandmarkObs_InRange, transf_obs);

		// reinit weight
		particles[i].weight = 1.0;
		
		for (int k = 0; k < transf_obs.size(); k++)
		{
			double map_x = transf_obs[k].x;
			double map_y = transf_obs[k].y;
			double mean_landmark_x;
			double mean_landmark_y;

			for (int m = 0; m < LandmarkObs_InRange.size(); m++)
			{
				if (transf_obs[k].id == LandmarkObs_InRange[m].id)
				{
					mean_landmark_x = LandmarkObs_InRange[m].x;
					mean_landmark_y = LandmarkObs_InRange[m].y;
					break;
				}
			}		

			double std_landmark_x = std_landmark[0];
			double std_landmark_y = std_landmark[1];
			double obs_w = (1 / (2 * M_PI*std_landmark_x*std_landmark_y)) * exp(-(pow(mean_landmark_x - map_x, 2) / (2 * pow(std_landmark_x, 2)) + (pow(mean_landmark_y - map_y, 2) / (2 * pow(std_landmark_y, 2)))));
			
			particles[i].weight *= obs_w;
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	vector<Particle> new_particles;

	// get all of the current weights
	vector<double> weights;
	for (int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
	}

	// generate random starting index for resampling wheel
	uniform_int_distribution<int> uni_particle_dist(0, num_particles - 1);
	int index = uni_particle_dist(gen);

	// get max weight
	double max_weight = *max_element(weights.begin(), weights.end());

	// uniform random distribution [0.0, max_weight)
	uniform_real_distribution<double> uni_weight_dist(0.0, max_weight);

	double beta = 0.0;

	// spin the resample wheel!
	for (int i = 0; i < num_particles; i++) {
		beta += uni_weight_dist(gen) * 2.0;
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
