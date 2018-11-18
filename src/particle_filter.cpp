/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <math.h>
#include <algorithm>
#include <chrono>
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
static map<int, LandmarkObs> landmarks;

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
    weights.push_back(particle.weight);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  auto start = chrono::high_resolution_clock::now();
  for (int i_p = 0; i_p < num_particles; i_p++) {
    if (yaw_rate == 0) {
      particles[i_p].x += velocity * delta_t * cos(particles[i_p].theta);
      particles[i_p].y += velocity * delta_t * sin(particles[i_p].theta);
    } else {
      particles[i_p].x += (velocity / yaw_rate) * (sin(particles[i_p].theta + (yaw_rate * delta_t)) - sin(particles[i_p].theta));
      particles[i_p].y += (velocity / yaw_rate) * (cos(particles[i_p].theta) - cos(particles[i_p].theta + (yaw_rate * delta_t)));
      particles[i_p].theta += yaw_rate * delta_t;
    }

    normal_distribution<double> dist_x(particles[i_p].x, std_pos[0]);
    normal_distribution<double> dist_y(particles[i_p].y, std_pos[1]);
    normal_distribution<double> dist_theta(particles[i_p].theta, std_pos[2]);

    particles[i_p].x = dist_x(generator);
    particles[i_p].y = dist_y(generator);
    particles[i_p].theta = dist_theta(generator);
  }

  auto end = chrono::high_resolution_clock::now();
  cout << "prediction(): " << chrono::duration_cast<chrono::microseconds>(end - start).count() << "us\n";
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.
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
  auto start = chrono::high_resolution_clock::now();

  for (int i_p = 0; i_p < num_particles; i_p++) {
    // convert landmarks list to map (for easier lookup)
    if (landmarks.empty()) {
      for (auto l : map_landmarks.landmark_list) {
        landmarks.insert(std::pair<int, LandmarkObs>(l.id_i, LandmarkObs{
                                                                 .id = l.id_i,
                                                                 .x = l.x_f,
                                                                 .y = l.y_f,
                                                             }));
      }
    }

    // transform observations to map coordinates
    double sin_theta = sin(particles[i_p].theta);
    double cos_theta = cos(particles[i_p].theta);
    std::vector<LandmarkObs> transformed;
    for (auto &o : observations) {
      double x_t = particles[i_p].x + (cos_theta * o.x) - (sin_theta * o.y);
      double y_t = particles[i_p].y + (sin_theta * o.x) + (cos_theta * o.y);
      LandmarkObs t = LandmarkObs{
          .id = 0,
          .x = x_t,
          .y = y_t,
      };

      // associate observation with nearest landmark
      std::map<int, double> distances;
      for (auto const &l : landmarks) {
        distances.insert(std::pair<int, double>(l.first, dist(t.x, t.y, l.second.x, l.second.y)));
      }

      auto nearest = std::min_element(distances.begin(), distances.end(),
                                      [](std::pair<int, double> a, std::pair<int, double> b) {
                                        return a.second < b.second;
                                      });
      t.id = nearest->first;

      // store transformed observation
      transformed.push_back(t);
    }

    //debug tool
    // vector<int> associations;
    // vector<double> sense_x;
    // vector<double> sense_y;

    // set weights
    particles[i_p].weight = 1.0;
    double std_x = std_landmark[0];
    double std_y = std_landmark[1];
    double gauss_norm = (2 * M_PI * std_x * std_y);
    double std_x_denom = (2 * pow(std_x, 2));
    double std_y_denom = (2 * pow(std_y, 2));

    for (LandmarkObs &t : transformed) {
      // associations.push_back(t.id);
      // sense_x.push_back(t.x);
      // sense_y.push_back(t.y);

      LandmarkObs landmark = landmarks[t.id];

      double x = t.x;
      double y = t.y;
      double mu_x = landmark.x;
      double mu_y = landmark.y;

      double exponent = (pow(x - mu_x, 2) / std_x_denom) + (pow(y - mu_y, 2) / std_y_denom);
      double prob = exp(-exponent) / gauss_norm;
      particles[i_p].weight *= prob;
    }
    weights[i_p] = particles[i_p].weight;

    // SetAssociations(particles[i_p], associations, sense_x, sense_y);
  }
  auto end = chrono::high_resolution_clock::now();
  cout << "updateWeights(): " << chrono::duration_cast<chrono::microseconds>(end - start).count() << "us\n";
}

void ParticleFilter::resample() {
  auto start = chrono::high_resolution_clock::now();
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::discrete_distribution<int> distribution(weights.begin(), weights.end());
  std::vector<Particle> resampled;
  for (int i_p = 0; i_p < num_particles; i_p++) {
    int i_r = distribution(generator);
    resampled.push_back(particles[i_r]);
  }

  particles = resampled;

  auto end = chrono::high_resolution_clock::now();
  cout << "resample(): " << chrono::duration_cast<chrono::microseconds>(end - start).count() << "us\n";
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
