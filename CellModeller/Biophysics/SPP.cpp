#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <math.h>
#include <random>
#include <tuple>
#include "SPP.h"


SPP::SPP(float dt_c, float gamma_c, float gamma_s_c, float radius_c, float W_s_c, float W_c_c, float f_pol_c, float D_r_c, float F_m_c, bool z_axis_c){
  this->dt = dt_c;
  this->gamma = gamma_c;
  this->gamma_s = gamma_s_c;
  this->radius = radius_c;
  this->W_s = W_s_c;
  this->W_c = W_c_c;
  this->f_pol = f_pol_c;
  this->D_r = D_r_c;
  this->F_m = F_m_c;
  this->z_axis = z_axis_c;
  this->dimensions = 2;
  if (z_axis_c){this->dimensions = 3;}
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  this->generator = std::default_random_engine(seed);
  this->normal = std::normal_distribution<float>(0.0,1.0);
}

void SPP::addCell(float pos_x, float pos_y, float pos_z, float polar_1, float polar_2){
  this->cell_centers.push_back(pos_x);
  this->cell_centers.push_back(pos_y);
  this->cell_centers.push_back(pos_z);
  this->cell_polarization.push_back(polar_1);
  this->cell_polarization.push_back(polar_2);
}

void SPP::step(){
  // find contacting pairs
  std::vector<std::tuple<int, int, float>> contacts = this->find_contacts();
  this->move_cells(contacts);
  // get new polarizations
  this->repolarize();
}

std::vector<std::tuple<int, int, float>> SPP::find_contacts(){
  // using brute force, could be optimized
  int cant_cells = this->cell_polarization.size()/2;
  std::vector<std::tuple<int, int, float>> contacts;
  float distance;
  for (int i=0; i<cant_cells; i++){
    for (int j=i+1; j<cant_cells; j++){
      printf("checking contacts\n");
      distance = 0;
      distance += pow((this->cell_centers[3*i]-this->cell_centers[3*j]), 2);
      distance += pow((this->cell_centers[3*i+1]-this->cell_centers[3*j+1]), 2);
      if (this->z_axis){distance += pow((this->cell_centers[3*i+2]-this->cell_centers[3*j+2]), 2);}
      distance = sqrt(distance);
      printf("The distance is %f\n", distance);
      if (distance <= 2*this->radius){
        printf("collision!\n");
        contacts.push_back(std::make_tuple(i, j, distance));
      }
    }
  }
  return contacts;
}

void SPP::move_cells(std::vector<std::tuple<int, int, float>> contacts){
  // get forces for each cell (eq 1, 2, 3)
  int cant_cells = this->cell_polarization.size()/2;
  int cant_contacts = contacts.size();
  int a, b, c;
  float force = 0;
  float contact_force;
  std::tuple<float, float, float> normal;
  for (int i=0; i<cant_cells; i++){
    normal = std::make_tuple(0, 0, 0);
    // calculate p_i, where cell_polarization[2*i] is the angle w.r.t the x axis in the xy plane
    // and cell_polarization[2*i+1] is the angle w.r.t the z axis (angles in radians)
    std::get<0>(normal) = this->F_m * cos(this->cell_polarization[2*i]);
    std::get<1>(normal) = this->F_m * sin(this->cell_polarization[2*i]);
    if (this->z_axis){
      std::get<0>(normal) *= sin(this->cell_polarization[2*i+1]);
      std::get<1>(normal) *= sin(this->cell_polarization[2*i+1]);
      std::get<2>(normal) = this->F_m * cos(this->cell_polarization[2*i+1]);
    }
    
    // contacts
    for (int j=0; j<cant_contacts; j++){
      if (std::get<0>(contacts[j])==i || std::get<1>(contacts[j])==i){
        // normal force
        a = std::get<0>(contacts[j]);
        b = std::get<1>(contacts[j]);
        if (b == i){c = b; b = a; a = c;}
        contact_force = 0;
        contact_force = std::get<2>(contacts[j]) - this->radius;
        contact_force *= (this->W_s + this->W_c) / this->radius;
        contact_force = 2 / this->radius * (this->W_s - contact_force);
        std::get<0>(normal) -= contact_force * (this->cell_centers[3*a]-this->cell_centers[3*b]) / std::get<2>(contacts[j]);
        std::get<1>(normal) -= contact_force * (this->cell_centers[3*a+1]-this->cell_centers[3*b+1]) / std::get<2>(contacts[j]);
        std::get<2>(normal) -= contact_force * (this->cell_centers[3*a+2]-this->cell_centers[3*b+2]) / std::get<2>(contacts[j]);
        // velocities * gamma
        // std::get<0>(normal) += this->gamma * (this->cell_velocities[3*a]-this->cell_velocities[3*b]);
        // std::get<1>(normal) += this->gamma * (this->cell_velocities[3*a+1]-this->cell_velocities[3*b+1]);
        // std::get<2>(normal) += this->gamma * (this->cell_velocities[3*a+2]-this->cell_velocities[3*b+2]);
      }
    }
    // divide by gamma_s
    std::get<0>(normal) /= this->gamma_s;
    std::get<1>(normal) /= this->gamma_s;
    std::get<2>(normal) /= this->gamma_s;

    // upgrade positions (pos + const*force)
    this->cell_centers[3*i] += std::get<0>(normal) * this->dt;
    this->cell_centers[3*i+1] += std::get<1>(normal) * this->dt;
    this->cell_centers[3*i+2] += std::get<2>(normal) * this->dt;
  }
}


void SPP::repolarize(){
  // assumes there´s no desired polarization
  int cant_cells = this->cell_polarization.size();

  // add noise
  std::vector<float> noise;
  for (int i=0; i<cant_cells; i++){
    noise.push_back(this->normal(this->generator));
  }
  // update polarization
  float delta_i;
  float C = sqrt(2 * this->D_r);
  for (int i=0; i<cant_cells/2; i++){
    delta_i = -this->f_pol * (this->cell_polarization[2*i]);
    delta_i += C * noise[2*i];
    this->cell_polarization[2*i] += delta_i;
  }
  if (z_axis){
    for (int i=0; i<cant_cells; i++){
      delta_i = -this->f_pol * (this->cell_polarization[2*i+1]);
      delta_i += C * noise[2*i+1];
      this->cell_polarization[2*i+1] += delta_i;
    }
  }

}
