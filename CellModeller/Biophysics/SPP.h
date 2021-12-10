#include <vector>
#include <random>
#include <tuple>
#define PI 3.1415


std::tuple<float, float, float> normalize(float x, float y, float z);
float point_plane_distance(std::vector<float> point, std::tuple<float, float, float> plane_norm, std::tuple<float, float, float> plane_point);
std::tuple<float, float> to_spherical(std::vector<float> dir);
float to_polar(std::vector<float> dir);
std::tuple<float, float, float> to_2Dcartesian(float theta);
std::tuple<float, float, float> to_3Dcartesian(float theta, float phi);

class SPP{
  private:
    bool z_axis;
    float gamma, gamma_s, W_s, W_c, f_pol, D_r, F_m;
    int dimensions;
    std::default_random_engine generator;
    std::normal_distribution<float> normal;
  public:
    std::vector<std::vector<float>> cell_centers;
    std::vector<std::vector<float>> cell_directions;
    std::vector<std::tuple<float, float, float>> plane_pts;
    std::vector<std::tuple<float, float, float>> plane_norms;
    float dt, radius;
    SPP(float dt_c=0.01, float gamma_c=1, float gamma_s_c=1, float radius_c=1, float W_s_c=1, float W_c_c=1, float f_pol_c=1, float D_r_c=0.001, float F_m_c=1, bool z_axis_c=true);
    void addCell(float pos_x=0, float pos_y=0, float pos_z=0, float dir_1=0, float dir_2=0, float dir_3=0);
    void addPlane(float p_x, float p_y, float p_z, float n_x, float n_y, float n_z);
    void step();
    std::vector<std::tuple<int, int, float>> find_cell_contacts();
    std::vector<std::tuple<int, int, float>> find_plane_contacts();
    void move_cells(std::vector<std::tuple<int, int, float>> cell_contacts, std::vector<std::tuple<int, int, float>> plane_contacts);
    void repolarize();
};