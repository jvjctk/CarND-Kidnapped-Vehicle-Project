/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  if (!initialized()) {
    num_particles = 100;  // TODO: Set the number of particles
    particles.resize(num_particles);

    // Gausian normal distributions
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_theta(theta, std[2]);

    // TODO: Add random Gaussian noise to each particle
    std::default_random_engine gen; //initializing randon engine

    for (int i = 0; i < num_particles; i++) {
      particles[i].id = i;
      particles[i].x = dist_x(gen);
      particles[i].y = dist_y(gen);
      particles[i].theta = dist_theta(gen);
      particles[i].weight = 1.0;
    }
    is_initialized = true;
  }


}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  // 0-centered Gausian distributions
  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);
  
  std::default_random_engine gen; //initializing randon engine


  for (int i = 0; i < particles.size(); i++) {
    // Convert from car coordinates to map coordinates using equation xxxx
    particles[i].x += velocity * cos(particles[i].theta) * delta_t;
    particles[i].y += velocity * sin(particles[i].theta) * delta_t;
    particles[i].theta += yaw_rate * delta_t;

    // Adding noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }


}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */


  double dist2_value, min_value_dist2;
  for (int i = 0; i < observations.size(); i++) { 
	// processing each observation

    min_value_dist2 = std::numeric_limits<double>::infinity(); // setting as infinity

    for (int j = 0; j < predicted.size(); j++) {
	// processing each prediction
      dist2_value = dist2(observations[i].x, observations[i].y, predicted[j].x,
                     predicted[j].y);

      if (dist2_value < min_value_dist2) {
        min_value_dist2 = dist2_value;
        observations[i].id = predicted[j].id; // assigning id
      }
    }
  }



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  double diff_x, diff_y;
  
  double gauss_norm =  1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
  double varx = 2 * pow(std_landmark[0], 2);
  double vary = 2 * pow(std_landmark[1], 2);

  double weight_sum;
  

  for (int i = 0; i < num_particles; i++) {
    

    vector<LandmarkObs> observations_; // creating object for converted observations
    observations_.resize(observations.size());
    
    // Converting to corresponding map coordinates with formula
    for (int j = 0; j < observations.size(); j++) {
      observations_[j].id = j;
      observations_[j].x = particles[i].x + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta)* observations[j].y;
      observations_[j].y = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
    }

    
    // considering particles which are in the range
    vector<LandmarkObs> observable_landmarks;
    for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      
      if ((fabs((particles[i].x - map_landmarks.landmark_list[j].x_f)) <= sensor_range) && (fabs((particles[i].y - map_landmarks.landmark_list[j].y_f)) <= sensor_range))
	{
          observable_landmarks.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f});
        }
    }

    // data association
    dataAssociation(observable_landmarks, observations_);

    // Compute the weight for each particle
    particles[i].weight = 1.0;
    for (int j = 0; j < observations_.size(); j++) {
      for (int k = 0; k < observable_landmarks.size(); k++) {
        if (observations_[j].id == observable_landmarks[k].id) {
          diff_x = observable_landmarks[k].x - observations_[j].x;
          diff_y = observable_landmarks[k].y - observations_[j].y;

          particles[i].weight *= gauss_norm * exp(-(pow(diff_x, 2) / varx + pow(diff_y, 2) / vary));
        }
      }
    }
    weight_sum += particles[i].weight;
  }

 // Normalizing weights

  for (int i = 0; i < particles.size(); i++) {
    particles[i].weight /= weight_sum;
  }



}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  vector<Particle> particles_; // creating object
  particles_.resize(num_particles);
  

  //finding maximum weight 
  double weight_max = 0;
  
  for (int i = 0; i < particles.size(); i++) {
    if (particles[i].weight > weight_max) {
      weight_max = particles[i].weight;
    }
  }

  
  std::uniform_real_distribution<double> random_weight(0.0, 2 * weight_max);
  
  std::uniform_int_distribution<> dist_index(0, particles.size() - 1);
  std::default_random_engine gen;
  int idx = dist_index(gen);
  double beta = 0;
  
  // resampling
  for (int i = 0; i < num_particles; i++) {
    beta += random_weight(gen);

    while (beta > particles[idx].weight) {
      beta -= particles[idx].weight;
      idx = (idx + 1) % num_particles;
    }
    resampled_particles[i] = particles[idx];
  }
  particles = resampled_particles;


}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}