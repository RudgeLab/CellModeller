#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <math.h>
#include <random>
#include <tuple>
#include "SPP.h"


std::tuple<float, float, float> normalize(float x, float y, float z){
  float normalize = pow(x, 2) + pow(y, 2) + pow(z, 2);
  normalize = sqrt(normalize);
  x = x/normalize;
  y = y/normalize;
  z = z/normalize;
  return std::make_tuple(x, y, z);
}

std::tuple<float, float> to_spherical(std::vector<float> dir){
  float theta, phi, arg;
  if (dir[0]==0){
    theta = PI/2.0;
    if (dir[1] < 0){
      theta *= -1.0;
    }
  }
  else{
    theta = atan(dir[1] / dir[0]);
    if (dir[0] > 0 && dir[1] < 0){
      theta += 2.0 * PI;
    }
    if (dir[0] < 0){
      theta += PI;
    }
  }
  if (dir[2]==0){phi = PI/2.0;}
  else{
    arg = sqrt(pow(dir[0], 2) + pow(dir[1], 2)) / dir[2];
    phi = atan(arg);
    if (dir[2]<0){phi += PI;}
  }
  std::tuple<float, float> theta_phi;
  theta_phi = std::make_tuple(theta, phi);
  return theta_phi;
}

float to_polar(std::vector<float> dir){
  float theta;
  if (dir[0]==0){
    theta = PI / 2.0;
    if (dir[1] < 0){theta += PI;}
  }
  else{
    theta = atan(dir[1] / dir[0]);
    if (dir[0] < 0){theta += PI;}
    if (dir[0] > 0 && dir[1] < 0){theta += 2.0 * PI;}
  }
  return theta;
}

std::tuple<float, float, float> to_2Dcartesian(float theta){
  float x, y, z;
  x = cos(theta);
  y = sin(theta);
  z = 0;
  std::tuple<float, float, float> abc;
  abc = std::make_tuple(x, y, z);
  return abc;
}

std::tuple<float, float, float> to_3Dcartesian(float theta, float phi){
  float x, y, z;
  x = cos(theta) * sin(phi);
  y = sin(theta) * sin(phi);
  z = cos(phi);
  std::tuple<float, float, float> abc;
  abc = std::make_tuple(x, y, z);
  return abc;
}

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

void SPP::addCell(float pos_x, float pos_y, float pos_z, float dir_1, float dir_2, float dir_3){
  std::vector<float> center;
  center.push_back(pos_x);
  center.push_back(pos_y);
  center.push_back(pos_z);
  this->cell_centers.push_back(center);
  std::vector<float> direction;
  direction.push_back(dir_1);
  direction.push_back(dir_2);
  direction.push_back(dir_3);
  this->cell_directions.push_back(direction);
}

void SPP::step(){
  // find contacting pairs
  std::vector<std::tuple<int, int, float>> contacts = this->find_contacts();
  this->move_cells(contacts);
  // get new directions
  printf("directions_x before replarizing are %f, %f\n", this->cell_directions[0][0], this->cell_directions[1][0]);
  this->repolarize();
  printf("directions_x after replarizing are %f, %f\n", this->cell_directions[0][0], this->cell_directions[1][0]);
}

std::vector<std::tuple<int, int, float>> SPP::find_contacts(){
  // using brute force, could be optimized
  int cant_cells = this->cell_directions.size();
  std::vector<std::tuple<int, int, float>> contacts;
  float distance;
  std::vector<float> center_i;
  std::vector<float> center_j;
  for (int i=0; i<cant_cells; i++){
    for (int j=i+1; j<cant_cells; j++){
      center_i = this->cell_centers[i];
      center_j = this->cell_centers[j];
      distance = 0;
      distance += pow((center_i[0]-center_j[0]), 2);
      distance += pow((center_i[1]-center_j[1]), 2);
      if (this->z_axis){distance += pow((center_i[2]-center_j[2]), 2);}
      distance = sqrt(distance);
      if (distance <= 2*this->radius){
        contacts.push_back(std::make_tuple(i, j, distance));
      }
    }
  }
  return contacts;
}

void SPP::move_cells(std::vector<std::tuple<int, int, float>> contacts){
  // get forces for each cell (eq 1, 2, 3)
  int cant_cells = this->cell_directions.size();
  int cant_contacts = contacts.size();
  int a, b, c;
  float force = 0;
  float contact_force;
  std::tuple<float, float, float> normal;
  std::vector<float>  direction;
  for (int i=0; i<cant_cells; i++){
    direction = this->cell_directions[i];
    normal = std::make_tuple(0, 0, 0);
    std::get<0>(normal) = this->F_m * direction[0];
    std::get<1>(normal) = this->F_m * direction[1];
    std::get<2>(normal) = this->F_m * direction[2];
    
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
        std::get<0>(normal) -= contact_force * (this->cell_centers[b][0]-this->cell_centers[a][0]) / std::get<2>(contacts[j]);
        std::get<1>(normal) -= contact_force * (this->cell_centers[b][1]-this->cell_centers[a][1]) / std::get<2>(contacts[j]);
        std::get<2>(normal) -= contact_force * (this->cell_centers[b][2]-this->cell_centers[a][2]) / std::get<2>(contacts[j]);
      }
    }
    // divide by gamma_s
    std::get<0>(normal) /= this->gamma_s;
    std::get<1>(normal) /= this->gamma_s;
    std::get<2>(normal) /= this->gamma_s;

    // upgrade positions
    this->cell_centers[i][0] += std::get<0>(normal) * this->dt;
    this->cell_centers[i][1] += std::get<1>(normal) * this->dt;
    this->cell_centers[i][2] += std::get<2>(normal) * this->dt;
  }
}


void SPP::repolarize(){
  // assumes there´s no desired direction
  int cant_cells = this->cell_directions.size();

  // add noise
  std::vector<float> noise;
  for (int i=0; i<2*cant_cells; i++){
    noise.push_back(this->normal(this->generator));
    //noise.push_back(0);
  }
  // update direction
  float delta_i;
  float C = sqrt(2 * this->D_r);
  float x, y, z, new_theta, new_phi;
  std::tuple<float, float, float> abc;
  std::tuple<float,float> theta_phi;
  std::vector<float> dir_i;
  for (int i=0; i<cant_cells; i++){
    dir_i = this->cell_directions[i];
    if (this->z_axis){
      theta_phi = to_spherical(dir_i);
      //delta_i = -this->f_pol * std::get<0>(theta_phi);
      //delta_i += C * noise[2*i];
      new_theta = std::get<0>(theta_phi) + C * noise[2*i];
      //delta_i = -this->f_pol * std::get<1>(theta_phi);
      //delta_i += C * noise[2*i+1];
      new_phi = std::get<1>(theta_phi) + C * noise[2*i+1];
      abc = to_3Dcartesian(new_theta, new_phi);
      dir_i[0] = std::get<0>(abc);
      dir_i[1] = std::get<1>(abc);
      dir_i[2] = std::get<2>(abc);
    }
    else {
      new_theta = to_polar(dir_i);
      //delta_i = -this->f_pol * new_theta;
      //delta_i += C * noise[i];
      new_theta += C * noise[i];
      abc = to_2Dcartesian(new_theta);
      dir_i[0] = std::get<0>(abc);
      dir_i[1] = std::get<1>(abc);
    }
    this->cell_directions[i][0] = dir_i[0];
    this->cell_directions[i][1] = dir_i[1];
    this->cell_directions[i][2] = dir_i[2];
  }
}