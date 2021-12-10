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

float point_plane_distance(std::vector<float> point, std::tuple<float, float, float> plane_norm, std::tuple<float, float, float> plane_point){
  // Calculates the distance between a point and a plane Ax+By+Cz+D=0
  float distance, A, B, C, D;
  A = std::get<0>(plane_norm);
  B = std::get<1>(plane_norm);
  C = std::get<2>(plane_norm);
  D = -(A * std::get<0>(plane_point) + B * std::get<1>(plane_point) + C * std::get<2>(plane_point));
  distance = abs(A * point[0] + B * point[1] + C * point[2] + D);
  distance /= sqrt(pow(A, 2) + pow(B, 2) + pow(C, 2));
  return distance;
}

std::tuple<float, float> to_spherical(std::vector<float> dir){
  //converts a normalized point (x,y,z) into its spherical coordinates
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
  //converts a normalized point (x,y,z) into its polar coordinates
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
  //converts an angle into the normalized point (x,y,0)
  float x, y, z;
  x = cos(theta);
  y = sin(theta);
  z = 0;
  std::tuple<float, float, float> abc;
  abc = std::make_tuple(x, y, z);
  return abc;
}

std::tuple<float, float, float> to_3Dcartesian(float theta, float phi){
  //converts a point from its spherical coordinates into the normalized point (x,y,z)
  float x, y, z;
  x = cos(theta) * sin(phi);
  y = sin(theta) * sin(phi);
  z = cos(phi);
  std::tuple<float, float, float> abc;
  abc = std::make_tuple(x, y, z);
  return abc;
}

SPP::SPP(float dt_c, float gamma_c, float gamma_s_c, float radius_c, float W_s_c, float W_c_c, float f_pol_c, float D_r_c, float F_m_c, bool z_axis_c){
  // class initializer
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
  // initialize the noise generator
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  this->generator = std::default_random_engine(seed);
  this->normal = std::normal_distribution<float>(0.0,1.0);
}

void SPP::addCell(float pos_x, float pos_y, float pos_z, float dir_1, float dir_2, float dir_3){
  // adds a new cell into the simulation
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

void SPP::addPlane(float p_x, float p_y, float p_z, float n_x, float n_y, float n_z){
  // adds a new plane to the simulation
  // the plane it's defined by a normal vector and a point in the plane
  std::tuple<float, float, float> point;
  point = std::make_tuple(p_x, p_y, p_z);
  this->plane_pts.push_back(point);
  std::tuple<float, float, float> norm;
  norm = std::make_tuple(n_x, n_y, n_z);
  this->plane_norms.push_back(norm);
}

void SPP::step(){
  // find contacting pairs
  std::vector<std::tuple<int, int, float>> cell_contacts = this->find_cell_contacts();
  // find plane collisions
  std::vector<std::tuple<int, int, float>> plane_contacts = this->find_plane_contacts();
  // move cells according to their polarization and collisions
  this->move_cells(cell_contacts, plane_contacts);
  // get new directions
  this->repolarize();
}

std::vector<std::tuple<int, int, float>> SPP::find_cell_contacts(){
  // using brute force, could be optimized
  int cant_cells = this->cell_directions.size();
  std::vector<std::tuple<int, int, float>> contacts;
  float distance;
  std::vector<float> center_i;
  std::vector<float> center_j;
  for (int i=0; i<cant_cells; i++){
    for (int j=i+1; j<cant_cells; j++){
      // we calculate the distance between the centers of each pair of cells
      center_i = this->cell_centers[i];
      center_j = this->cell_centers[j];
      distance = 0;
      distance += pow((center_i[0]-center_j[0]), 2);
      distance += pow((center_i[1]-center_j[1]), 2);
      if (this->z_axis){distance += pow((center_i[2]-center_j[2]), 2);}
      distance = sqrt(distance);
      if (distance <= 2*this->radius){
        // if they are close enough, we save the contact as a tuple (i,j,d),
        // where i and j are the indexes of the contacting cells and d is the distance between their centers
        contacts.push_back(std::make_tuple(i, j, distance));
      }
    }
  }
  return contacts;
}

std::vector<std::tuple<int, int, float>> SPP::find_plane_contacts(){
  int cant_cells = this->cell_directions.size();
  int cant_planes = this->plane_norms.size();
  std::vector<std::tuple<int, int, float>> contacts;
  float distance;
  std::vector<float> center_i;
  for (int i=0; i<cant_cells; i++){
    for (int j=0; j<cant_planes; j++){
      // we compute the distance between each cell and plane
      distance = point_plane_distance(this->cell_centers[i], this->plane_norms[j], this->plane_pts[j]);
      if (distance <= this->radius){
        // if they are close enough, we save the contact as a tuple (i, j, d), where:
        // i is the index of the cell, j is the index of the plane, and d is the perpendicular distance between them
        contacts.push_back(std::make_tuple(i, j, distance));
      }
    }
  }
  return contacts;
}

void SPP::move_cells(std::vector<std::tuple<int, int, float>> cell_contacts, std::vector<std::tuple<int, int, float>> plane_contacts){
  // get forces for each cell (eq 1, 2, 3)
  int cant_cells = this->cell_directions.size();
  int cant_contacts = cell_contacts.size();
  int cant_plane_contacts = plane_contacts.size();
  int a, b, c;
  float force = 0;
  float contact_force;
  std::tuple<float, float, float> movement;
  std::vector<float>  direction;
  for (int i=0; i<cant_cells; i++){
    direction = this->cell_directions[i];
    // the movement is given by the polarization of the cell...
    movement = std::make_tuple(0, 0, 0);
    std::get<0>(movement) = this->F_m * direction[0];
    std::get<1>(movement) = this->F_m * direction[1];
    std::get<2>(movement) = this->F_m * direction[2];
    
    // and its contacts/collisions
    // contacts
    for (int j=0; j<cant_contacts; j++){
      if (std::get<0>(cell_contacts[j])==i || std::get<1>(cell_contacts[j])==i){
        // normal force
        a = std::get<0>(cell_contacts[j]);
        b = std::get<1>(cell_contacts[j]);
        // we want a to be the index of the cell we are checking right now
        if (b == i){c = b; b = a; a = c;}
        // calculate the force of this contact
        contact_force = std::get<2>(cell_contacts[j]) - this->radius;
        contact_force *= (this->W_s + this->W_c) / this->radius;
        contact_force = 2 / this->radius * (this->W_s - contact_force);

        // this force is given in the direction of the normal vector between both cells
        std::get<0>(movement) -= contact_force * (this->cell_centers[b][0]-this->cell_centers[a][0]) / std::get<2>(cell_contacts[j]);
        std::get<1>(movement) -= contact_force * (this->cell_centers[b][1]-this->cell_centers[a][1]) / std::get<2>(cell_contacts[j]);
        std::get<2>(movement) -= contact_force * (this->cell_centers[b][2]-this->cell_centers[a][2]) / std::get<2>(cell_contacts[j]);
      }
    }

    // plane contacts
    for (int j=0; j<cant_plane_contacts; j++){
      if (std::get<0>(plane_contacts[j]) == i){
        contact_force = std::get<2>(plane_contacts[j]); - this->radius;
        contact_force *= (this->W_s + this->W_c) / this->radius;
        contact_force = 2 / this->radius * (this->W_s - contact_force);
        // this force is given in the direction of the normal vector between the cell and the plane
        std::get<0>(movement) -= contact_force * std::get<0>(this->plane_norms[std::get<1>(plane_contacts[j])]);
        std::get<1>(movement) -= contact_force * std::get<1>(this->plane_norms[std::get<1>(plane_contacts[j])]);
        std::get<2>(movement) -= contact_force * std::get<2>(this->plane_norms[std::get<1>(plane_contacts[j])]);
      }
    }


    // divide by gamma_s
    std::get<0>(movement) /= this->gamma_s;
    std::get<1>(movement) /= this->gamma_s;
    std::get<2>(movement) /= this->gamma_s;

    // upgrade positions
    this->cell_centers[i][0] += std::get<0>(movement) * this->dt;
    this->cell_centers[i][1] += std::get<1>(movement) * this->dt;
    this->cell_centers[i][2] += std::get<2>(movement) * this->dt;
  }
}


void SPP::repolarize(){
  // assumes there´s no desired direction


  int cant_cells = this->cell_directions.size();

  // create noise
  std::vector<float> noise;
  for (int i=0; i<2*cant_cells; i++){
    noise.push_back(this->normal(this->generator));
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
    // if we are in three dimensions, we work in spherical coordinates, if not we work in polar coordinates
    if (this->z_axis){
      // convert to spherical coordinates
      theta_phi = to_spherical(dir_i);
      // add noise in both angles
      new_theta = std::get<0>(theta_phi) + C * noise[2*i];
      new_phi = std::get<1>(theta_phi) + C * noise[2*i+1];
      // go back to cartesian coordinates
      abc = to_3Dcartesian(new_theta, new_phi);
      dir_i[0] = std::get<0>(abc);
      dir_i[1] = std::get<1>(abc);
      dir_i[2] = std::get<2>(abc);
    }
    else {
      // convert to polar coordinates
      new_theta = to_polar(dir_i);
      // add_noise to the angle
      new_theta += C * noise[i];
      // go back to cartesian coordinates
      abc = to_2Dcartesian(new_theta);
      dir_i[0] = std::get<0>(abc);
      dir_i[1] = std::get<1>(abc);
    }
    // update the direction
    this->cell_directions[i][0] = dir_i[0];
    this->cell_directions[i][1] = dir_i[1];
    this->cell_directions[i][2] = dir_i[2];
  }
}