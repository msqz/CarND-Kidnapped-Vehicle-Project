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
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>

#include "particle_filter.h"

using namespace std;

static default_random_engine generator;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[3]);
  num_particles = 10;
  for (int i_p = 0; i_p < num_particles; i_p++) {
    Particle particle = {
        .id = 0,
        .x = dist_x(generator),
        .y = dist_y(generator),
        .theta = dist_theta(generator),
        .weight = 1.0,
    };

    particles.push_back(particle);
    weights[i_p] = 1.0;
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  for (Particle &particle : particles) {
    if (yaw_rate == 0) {
      particle.x += velocity * delta_t * cos(particle.theta);
      particle.y += velocity * delta_t * sin(particle.theta);
    } else {
      particle.x += (velocity / yaw_rate) * (sin(particle.theta + (yaw_rate * delta_t)) - sin(particle.theta));
      particle.y += (velocity / yaw_rate) * (cos(particle.theta) - cos(particle.theta + (yaw_rate * delta_t)));
      particle.theta += yaw_rate * delta_t;
    }

    normal_distribution<double> dist_x(particle.x, std_pos[0]);
    normal_distribution<double> dist_y(particle.y, std_pos[1]);
    normal_distribution<double> dist_theta(particle.theta, std_pos[2]);

    particle.x = dist_x(generator);
    particle.y = dist_y(generator);
    particle.theta = dist_theta(generator);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.
}

std::vector<LandmarkObs> getObservationsInRange(const Particle &particle, double range, const std::vector<LandmarkObs> &observations) {
  std::vector<LandmarkObs> enclosed;

  for (const auto &observation : observations) {
    if (dist(particle.x, observation.x, particle.y, observation.y) <= range) {
      enclosed.push_back(observation);
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

  for (Particle &particle : particles) {
    std::vector<LandmarkObs> scope = observations;  //getObservationsInRange(particle, sensor_range, observations);
    // transform observations to map coordinates
    double sin_theta = sin(particle.theta);
    double cos_theta = cos(particle.theta);
    std::vector<LandmarkObs> transformed;
    for (auto &s : scope) {
      double x_t = particle.x + (cos_theta * s.x) - (sin_theta * s.y);
      double y_t = particle.y + (sin_theta * s.x) + (cos_theta * s.y);
      transformed.push_back(LandmarkObs{
          .id = 0,
          .x = x_t,
          .y = y_t,
      });
    }

    // convert landmarks list to map (for easier lookup)
    std::map<int, LandmarkObs> landmarks;
    for (auto l : map_landmarks.landmark_list) {
      landmarks.insert(std::pair<int, LandmarkObs>(l.id_i, LandmarkObs{
                                                               .id = l.id_i,
                                                               .x = l.x_f,
                                                               .y = l.y_f,
                                                           }));
    }

    // associate observations with nearest landmarks
    for (LandmarkObs &t : transformed) {
      std::map<int, double> distances;
      for (auto const &l : landmarks) {
        distances.insert(std::pair<int, double>(l.first, dist(t.x, t.y, l.second.x, l.second.y)));
      }
      auto nearest = std::min_element(distances.begin(), distances.end(),
                                      [](std::pair<int, double> a, std::pair<int, double> b) {
                                        return a.second < b.second;
                                      });
      t.id = nearest->first;
    }

    //debug tool
    // vector<int> associations;
    // vector<double> sense_x;
    // vector<double> sense_y;
    // set weights
    particle.weight = 1.0;
    for (LandmarkObs &t : transformed) {
      // associations.push_back(t.id);
      // sense_x.push_back(t.x);
      // sense_y.push_back(t.y);

      LandmarkObs landmark = landmarks[t.id];

      double std_x = std_landmark[0];
      double std_y = std_landmark[1];
      double mu_x = landmark.x;
      double mu_y = landmark.y;
      double x = t.x;
      double y = t.y;

      double gauss_norm = (2 * M_PI * std_x * std_y);
      double exponent = (pow(x - mu_x, 2) / (2 * pow(std_x, 2))) + (pow(y - mu_y, 2) / (2 * pow(std_y, 2)));
      double prob = exp(-exponent) / gauss_norm;
      particle.weight *= prob;
    }

    //SetAssociations(particle, associations,sense_x, sense_y);
  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::vector<double> weights_tmp;
  for (const Particle &p : particles) {
    weights_tmp.push_back(p.weight);
  }
  // TODO all zeros in weights_tmp
  std::discrete_distribution<int> distribution(weights_tmp.begin(), weights_tmp.end());
  std::vector<Particle> resampled;
  for (int i_p = 0; i_p < num_particles; i_p++) {
    int i_r = distribution(generator);
    resampled.push_back(particles[i_r]);
  }

  particles = resampled;
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

  return particle;
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
