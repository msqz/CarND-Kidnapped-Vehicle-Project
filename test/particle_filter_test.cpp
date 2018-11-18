#include "../src/particle_filter.h"
#include <iostream>
#include "../src/helper_functions.h"
#include "../src/map.h"

void testInit() {
  bool result = true;

  ParticleFilter pf;
  double x = 1.0;
  double y = 1.0;
  double theta = 0.1;
  double stdev[2] = {0.0, 0.0};
  pf.init(x, y, theta, stdev);

  for (const Particle &p : pf.particles) {
    if (p.x == 0) {
      std::cout << "uninitalized particle\n";
      result = false;
    }
    if (p.y == 0) {
      std::cout << "uninitalized particle\n";
      result = false;
    }
    if (p.theta == 0) {
      std::cout << "uninitalized particle\n";
      result = false;
    }
    if (p.weight != 1.0) {
      std::cout << "unitialized particle\n";
      result = false;
    }
  }

  if (result == false) {
    throw - 1;
  }
  std::cout << __FUNCTION__ << " passed\n";
}

void testPrediction() {
  bool result = true;

  ParticleFilter pf;
  double x = 1.0;
  double y = 1.0;
  double theta = 0.1;
  double stdev[2] = {0.0, 0.0};
  pf.init(x, y, theta, stdev);

  double delta_t = 1.0;
  double std_pos[2] = {0.0, 0.0};
  double velocity = 1.0;
  double yaw_rate = 0.1;
  pf.prediction(delta_t, std_pos, velocity, yaw_rate);

  double x_expected = 1.99;
  double y_expected = 1.15;
  double theta_expected = 0.2;

  for (const Particle &particle : pf.particles) {
    if (round(particle.x * 100) / 100 != x_expected) {
      std::cout << "x should be " << x_expected << " got " << particle.x << "\n";
      result = false;
    }
    if (round(particle.y * 100) / 100 != y_expected) {
      std::cout << "y should be " << y_expected << " got " << particle.y << "\n";
      result = false;
    }
    if (round(particle.theta * 100) / 100 != theta_expected) {
      std::cout << "theta should be " << theta_expected << " got " << particle.theta << "\n";
      result = false;
    }
  }

  if (result == false) {
    throw - 1;
  }
  std::cout << __FUNCTION__ << " passed\n";
}

void testUpdate() {
  bool result = true;

  ParticleFilter pf;
  double x = 1.0;
  double y = 1.0;
  double theta = 0.1;
  double stdev[2] = {0.0, 0.0};
  pf.init(x, y, theta, stdev);

  double sensor_range = 50;
  double std_landmark[2] = {0.0, 0.0};
  std::vector<LandmarkObs> observations{LandmarkObs{
      .id = 0,
      .x = 2.0,
      .y = 2.0,
  }};

  Map map_landmarks;
  map_landmarks.landmark_list.push_back(Map::single_landmark_s{
    .id_i = 1,
    .x_f = 2.0,
    .y_f = 2.0,
  });

  pf.updateWeights(sensor_range, std_landmark, observations, map_landmarks);
  double weight_expected = 0.0;

  for (const Particle &particle : pf.particles) {
    if (particle.weight != weight_expected) {
      std::cout << "weight should be " << weight_expected << " got " << particle.weight << "\n";
      result = false;
    }
  }

  if (result == false) {
    throw - 1;
  }
  std::cout << __FUNCTION__ << " passed\n";
}

int main() {
  testInit();
  testPrediction();
  testUpdate();
  return 0;
}