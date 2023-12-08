/* -------------------------------------------------------
AppPotts_eng class header for abnormal grain growth
--
Read-Shockley implementation developed by Efrain Hernandez-Rivera (2017--2018)
US Army Research Laboratory
--
THIS SOFTWARE IS MADE AVAILABLE ON AN "AS IS" BASIS
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, NEITHER
EXPRESSED OR IMPLIED
------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(potts/eng,AppPotts_eng)

#else

#ifndef SPK_APP_POTTS_ENG_H
#define SPK_APP_POTTS_ENG_H

#include <map>
#include <vector>
#include <unordered_map>
#include <string>
#include <set>
#include "app_potts.h"
//#include <Eigen/Geometry>

namespace SPPARKS_NS {

class AppPotts_eng : public AppPotts {
  public:

    AppPotts_eng(class SPPARKS *, int, char **);
    ~AppPotts_eng();
    void input_app(char *, int, char **);
    void init_app();
    void grow_app();
    void site_event_rejection(int, class RandomPark *);
    void init_vector(int, int, double, double, double, double);

    void inclination_calculation(int i, int si, std::vector<double> & vector_lab_frame, int center);
    void compute_normal_vector_2DLinear_matrix(int i, int si, std::vector<double> & vector_lab_frame, int center);
    int int_mod(int base, int divide);
    void vector_average(std::vector<double> vector_lab_frame_1, std::vector<double> vector_lab_frame_2, std::vector<double> & average_normal_lab_frame);
    void compute_normal_vector_3DLinear_matrix(int i, int si, std::vector<double> & vector_lab_frame, int center);

    void convert_inclination_refercrystal(int i, int i0, std::vector<double> & average_normal_lab_frame, std::vector<double> & inclination);
    void vector_to_quternion(std::vector<double> & normals, std::vector<double> & reference_normals, std::vector<double> & q);

    void euler2quat(std::vector<double> eulerAngle, std::vector<double> & quternion_result);
    void mat2quat(const double O[3][3], double q[4]);
    void symmat(double ***sym);
    void quat_mult(const double qi[4], const double qj[4], double q[4]);
    void quaternions(const std::vector<double> qi, const std::vector<double> qj, std::vector<double> & q_result);

    // new neighbor function
    void get_nearest_neighbor(int i, int neighbor_size);
    void get_boarder_place(int center, int & i_dis, int & j_dis, int & k_dis);

    void read_misorientation();
    void init_energy();
    void read_energy();
    void compute_inclination_bins();
    void init_inclination_table();
    void read_inclination_table();
    double compute_inclination_energy(std::vector<double> & input_i_frame);
    double energy_ave_split(int i, int nei);
    void read_linear_vector_matrix();
    double site_energy(int);





  protected:
    double *phi1,*phi2,*Phi; //pointer to 3 rotation angles
    int *spin;
    int dimension;  //simulation dimension
    double thetam; //High-low angle divider
    double Jij; //Interaction energy
    int Osym; //Symmetry Operator flag
    double **symquat;
    //smooth algorithm parameters
    int interval, neighbor_length;
    std::string smthAlgo;
    int clip;
    int cleanstep=0;
    float storage_flag;
    int num_read=0;
    int num_calc=0;
    float test_error_incl=0;
    // double **symquat; //Symmetry Operator in quaternion space
    int nx,ny,nz;
    // The calculated normal vector of site i or nei based on i's grain in lab frame
    std::vector<std::vector<double>> vector_i_lab_frame;
    // The calculated normal vector of site i or nei based on nei's grain in lab frame
    std::vector<std::vector<double>> vector_nei_lab_frame;
    // The averaged normal vector between site i and nei based on i's grain in lab frame
    std::vector<std::vector<double>> average_normal_i_lab_frame;
    // The averaged normal vector between site i and nei based on nei's grain in lab frame
    std::vector<std::vector<double>> average_normal_nei_lab_frame;
    std::vector<std::vector<double>> average_normal_lab_frame_tri;
    // The inclination between site i and nei based on i's grain in i's frame (after convert)
    std::vector<std::vector<double>> inclination_i_frame;
    // The inclination between site i and nei based on nei's grain in nei's frame (after convert)
    std::vector<std::vector<double>> inclination_nei_frame;

    std::vector<double> reference_axis = {0, 0, 0};
    std::vector<double> reference_quaternion = {1, 0, 0, 0};
    std::vector<std::vector<double>> inclination_quaternion_i_frame;
    std::vector<std::vector<double>> inclination_quaternion_nei_frame;
    std::vector<std::vector<double>> inclination_quaternion_symmetry_i_frame;
    std::vector<std::vector<double>> inclination_quaternion_symmetry_nei_frame;

    // inclination energy
    double maxIncEnergy;
    int numBins;
    double E_delta;
    double E_m;
    double current_gamma = 0.0;

    // triple energy
    std::set<int> nei_set;
    std::string triple_energy;

    //map to store misorientations
    std::map<std::pair<int,int>, double> misos;
    std::map<std::pair<int,int>, double> energy;
    std::vector<std::vector<int>> structure_matrix_global; //2d site ID matrix for L2
    std::vector<std::vector<std::vector<int>>> structure_matrix_global_3d;//3d site ID matrix for L2
    // vector matrix for L2
    std::vector<std::vector<double>> linear_vector_matrix_i;
    std::vector<std::vector<double>> linear_vector_matrix_j;
    std::vector<std::vector<std::vector<double>>> linear_vector_matrix_i_3d;
    std::vector<std::vector<std::vector<double>>> linear_vector_matrix_j_3d;
    std::vector<std::vector<std::vector<double>>> linear_vector_matrix_k_3d;
    // inclination energy
    std::vector<std::vector<double>> eulers;
    std::vector<std::vector<double>> incTable;
    std::vector<double> incBins;
    // store inclination vector
    std::vector< std::unordered_map<int, std::vector<double> > > inc_storage;

    // New neighbor functions
    int numneigh_near;
    std::vector<int> neighbor_near;


  private:
   double pfraction;
   int multi,nthresh;
};

}

#endif
#endif
