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
#include <limits>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
        default_random_engine gen;

        //setting standard deviation
        double std_x = std[0];
        double std_y = std[1];
        double std_theta = std[2];

        //Create normal distribution
        normal_distribution<double> dist_x(x, std_x);
        normal_distribution<double> dist_y(y, std_y);
        normal_distribution<double> dist_theta(theta, std_theta);
        
        num_particles = 500;
        for (int i = 0; i < num_particles; i++){
            Particle P;
            P.id = i;
            P.x = dist_x(gen);
            P.y = dist_y(gen);
            P.theta = dist_theta(gen);
            P.weight = 1;
            particles.push_back(P);
        }
        is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
        default_random_engine gen;

        //setting standard deviation
        double std_x = std_pos[0];
        double std_y = std_pos[1];
        double std_theta = std_pos[2];

        //Create normal distribution
        normal_distribution<double> dist_x(0, std_x);
        normal_distribution<double> dist_y(0, std_y);
        normal_distribution<double> dist_theta(0, std_theta);

        if (yaw_rate == 0){
            for(int i = 0; i<num_particles; i++){
                Particle P = particles[i];


                P.x = P.x + cos(P.theta) * velocity * delta_t + dist_x(gen);
                P.y = P.y + sin(P.theta) * velocity * delta_t + dist_y(gen);
                P.theta = P.theta +  dist_theta(gen);
                while (P.theta < 0){
                    P.theta += 2 * M_PI;
                }
                while(P.theta > 2 * M_PI){
                    P.theta -= 2 * M_PI;
                }
                particles[i] = P;
            }
        }
        else{
            for(int i = 0; i<num_particles; i++){
                Particle P = particles[i];

                P.x = P.x + (sin(P.theta + yaw_rate * delta_t) - sin(P.theta)) * velocity / yaw_rate + dist_x(gen);
                P.y = P.y + (cos(P.theta) - cos(P.theta + yaw_rate * delta_t)) * velocity / yaw_rate + dist_y(gen);
                P.theta = P.theta + yaw_rate * delta_t +  dist_theta(gen);
                while (P.theta < 0){
                    P.theta += 2 * M_PI;
                }
                while(P.theta > 2 * M_PI){
                    P.theta -= 2 * M_PI;
                }
                particles[i] = P;
            }
        }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
        /*
         * This part is implemented in updateWeights method
         */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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
        double mindist;
        double dist;
        float sig_x = std_landmark[0];
        float sig_y = std_landmark[1];

        for(int partind = 0; partind < num_particles; partind++) {

            Particle P = particles[partind];
            P.weight = 1;
            vector<int> associations;
            vector<double> sense_x;
            vector<double> sense_y;

            for(int i = 0; i < observations.size(); i++){
                LandmarkObs obs = observations[i];
                Map::single_landmark_s closest_landmark;

                //convert observation to map space
                double xpos = P.x + obs.x * cos(P.theta) - obs.y * sin(P.theta);
                double ypos = P.y + obs.x * sin(P.theta) + obs.y * cos(P.theta);
                mindist = numeric_limits<double>::max();
                //finding closest landmark
                for(int j = 0; j< map_landmarks.landmark_list.size(); j++){
                    Map::single_landmark_s landmark = map_landmarks.landmark_list[j];
                    float lm_x = landmark.x_f;
                    float lm_y = landmark.y_f;
                    int lm_id = landmark.id_i;
                    dist = sqrt((lm_x - xpos) * (lm_x - xpos) + (lm_y - ypos) * (lm_y - ypos));
                    if(dist < mindist){
                        mindist = dist;
                        closest_landmark = landmark;
                    }
                }
                double gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
                double exponent = ((closest_landmark.x_f - xpos)*(closest_landmark.x_f - xpos))/(2 * sig_x*sig_x)
                                + ((closest_landmark.y_f - ypos)*(closest_landmark.y_f - ypos))/(2 * sig_y*sig_y);

                P.weight *= gauss_norm * exp(-exponent);
                associations.push_back(closest_landmark.id_i);
                sense_x.push_back(xpos);
                sense_y.push_back(ypos);
            }
            P = SetAssociations(P, associations, sense_x, sense_y);
            particles[partind] = P;
        }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
        vector<Particle> newParticles;

        default_random_engine gen;
        uniform_int_distribution<int> indexgen(0, num_particles);

        int index = indexgen(gen);
        double beta = 0;
        double max_weight=0;
        for (int i = 0; i < num_particles; i++){
            Particle P = particles[i];
            if (P.weight>max_weight) {
                max_weight = P.weight;
            }
        }
        uniform_real_distribution<double> realgen(0, 2 * max_weight);
        for (int i = 0; i < num_particles; i++){
            beta += realgen(gen);
            while(beta > particles[index].weight){
                beta -= particles[index].weight;
                index = (index + 1) % num_particles;
            }
            newParticles.push_back(particles[index]);
        }
        particles = newParticles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

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
