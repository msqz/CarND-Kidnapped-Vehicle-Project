/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <string>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  std::default_random_engine generator;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[3]);
  num_particles = 1000;
  for (int i = 0; i < num_particles; i++) {
    Particle particle = {
        .id = i + 1,
        .x = dist_x(generator),
        .y = dist_y(generator),
        .theta = dist_theta(generator),
        .weight = 1,
    };

    particles.push_back(particle);
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  std::default_random_engine generator;
  for (Particle particle : particles) {
    if (yaw_rate == 0) {
      double new_x = velocity * delta_t * cos(particle.theta);
      std::normal_distribution dist_x(new_x, std_pos[0]);
      particle.x += dist_x(generator);

      double new_y = velocity * delta_t * sin(particle.theta);
      std::normal_distribution dist_y(new_y, std_pos[1]);
      particle.y += dist_y(generator);
    } else {
      double new_x = (velocity / yaw_rate) * (sin(particle.theta) * yaw_rate * delta_t - sin(particle.theta));
      std::normal_distribution dist_x(new_x, std_pos[0]);
      particle.x += dist_x(generator);

      double new_y = (velocity / yaw_rate) * (cos(particle.theta) - cos(particle.theta) + yaw_rate * delta_t);
      std::normal_distribution dist_y(new_y, std_pos[1]);
      particle.y += dist_y(generator);

      double new_theta = yaw_rate * delta_t;
      std::normal_distribution dist_theta(new_theta, std_pos[2]);
      particle.theta += dist_theta(generator);
    }
  }
}

void ParticleFilter::dataAssociation(const Map &map, std::vector<LandmarkObs> &observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.
  for (LandmarkObs obs : observations) {
    double tmp_min = 0.0;
    for (Map::single_landmark_s landmark : map.landmark_list) {
      double distance = dist(obs.x, landmark.x_f, obs.y, landmark.y_f);
      if (distance < tmp_min || obs.id == 0) {
        tmp_min = distance;
        obs.id = landmark.id_i;
      }
    }
  }
}

std::vector<LandmarkObs> getObservationsInRange(const Particle &particle, double range, const std::vector<LandmarkObs> &observations) {
  std::vector<LandmarkObs> enclosed;

  for (LandmarkObs obs : observations) {
    if (dist(particle.x, obs.x, particle.y, obs.y) <= range) {
      enclosed.push_back(obs);
    }
  }

  return enclosed;
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

  for (Particle particle : particles) {
    std::vector<LandmarkObs> enclosed_obs = getObservationsInRange(particle, sensor_range, observations);
    // transform -- every particle sees in sensor_range, so it should have assigned observations only in that range
    // transformation is separate for each particle
    for (LandmarkObs obs : enclosed_obs) {
      obs.x = particle.x + cos(particle.theta) * obs.x - sin(particle.theta) * obs.y;
      obs.y = particle.y + sin(particle.theta) * obs.x + cos(particle.theta) * obs.y;
    }
    // associate -- between transformed observations and landmarks
    dataAssociation(map_landmarks, enclosed_obs);
    double weight = 1;

    // weight -- product of probabilties (exp distance between each observation and landmark)
    // those with many hight probabilities will survive resampling phase
    for (LandmarkObs obs: enclosed_obs) {
      double numerator = exp(-(()); ........
      double denominator = 0.0;
      double prob = numerator/denominator;
      weight *= prob;
    }
  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
}

Particle ParticleFilter::SetAssociations(Particle &particle, const std::vector<int> &associations,
                                         const std::vector<double> &sense_x, const std::vector<double> &sense_y) {
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseX(Particle best) {
  vector<double> v = best.sense_x;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best) {
  vector<double> v = best.sense_y;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
