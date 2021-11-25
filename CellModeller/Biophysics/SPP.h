#include <vector>
#include <random>
#include <tuple>
#define PI 3.1415


std::tuple<float, float, float> normalize(float x, float y, float z=0);

class SPP{
  private:
    bool z_axis;
    float gamma, gamma_s, W_s, W_c, f_pol, D_r, F_m;
    int dimensions;
    std::default_random_engine generator;
    std::normal_distribution<float> normal;
  public:
    std::vector<float> cell_centers;
    std::vector<float> cell_directions;
    float dt, radius;
    SPP(float dt_c=0.01, float gamma_c=1, float gamma_s_c=1, float radius_c=1, float W_s_c=1, float W_c_c=0, float f_pol_c=1, float D_r_c=1, float F_m_c=1, bool z_axis_c=true);
    void addCell(float pos_x=0, float pos_y=0, float pos_z=0, float dir_1=0, float dir_2=0, float dir_3=PI);
    void step();
    std::vector<std::tuple<int, int, float>> find_contacts();
    void move_cells(std::vector<std::tuple<int, int, float>> contacts);
    void repolarize();
};