#include <amino.h>
// #include <amino/mat.h>
// #include "amino/rx/rxtype.h"
#include <amino/rx/scenegraph.h>
#include <amino/rx/scene_plugin.h>
#include <amino/rx/scene_win.h>
#include <amino/rx/scene_fk.h>
// #include <amino/rx/scene_sdl.h>
// #include "amino/rx/scene_gl.h"
// #include "amino/rx/scene_geom.h"
#include <cmath>
#include <unistd.h> 
#include <random>
#include <vector>
#include <array>
#include <chrono>
#include <float.h>

static inline void
aafeq(const char *name, double a, double b, double tol)
{
    if (!aa_feq(a, b, tol)) {
        fprintf(stderr, "FAILED: %s\n", name);
        fprintf(stderr, "a: %f, b: %f\n", a, b);
        abort();
    }
}

static inline bool
aafeq2(const char *name, double a, double b, double tol)
{
    if (!aa_feq(a, b, tol)) {
        fprintf(stderr, "FAILED: %s\n", name);
        fprintf(stderr, "a: %f, b: %f\n", a, b);
        return false;
    }
    return true;
}

static inline void
admeq(const char *name, const double *A, const double *B, double tol, size_t size)
{
    aafeq( name, aa_la_ssd(size,A,B), 0, tol );
}

static inline bool
admeq2(const char *name, const double *A, const double *B, double tol, size_t size)
{
    return aafeq2( name, aa_la_ssd(size,A,B), 0, tol );
}

static inline void
mat_print_raw(const double *data, const int rows, const int cols)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            std::cout << data[rows*j+i] << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

static inline void
mat_print_pretty(const double *data, const int rows, const int cols)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double x = data[rows*j+i];
            if (fabs(x) < 1e-4) {
                std::cout << 0 << "  ";
            } else {
                std::cout << x << "  ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

static inline void
array_print(const double *data, const int n)
{
    for (int i = 0; i < n; i++) {
        std::cout << data[i];
        std::cout << "\n";
    }
    std::cout << "\n";
}

static bool
check_joint_limits(std::array<double, 7> &sol, 
                   const double joint_limits[14])
{
    for (int i = 0; i < sol.size(); i++) {
        double lower_limit = joint_limits[i*2];
        double upper_limit = joint_limits[i*2+1];
        int attempts = 0;
        while (attempts < 10 && sol[i] < lower_limit) {
            sol[i] += 2*M_PI;
            ++attempts;
        }
        attempts = 0;
        while (attempts < 10 && sol[i] > upper_limit) {
            sol[i] -= 2*M_PI;
            ++attempts;
        }
        if (sol[i] < lower_limit || sol[i] > upper_limit) {
            return false; // if any angle in the solution fails
        } 
    }

    return true; // if all angles in the solution are in limits
}

static void
randomize_config(double config_data[7],
                 const double joint_limits[14])
{
    std::random_device rd;
    std::mt19937 gen(rd());  // here you could set the seed, but std::random_device already does that
    for (int i = 0; i < 7; i++) {
        double lower_lim = joint_limits[i*2];
        double upper_lim = joint_limits[i*2+1];
        std::uniform_real_distribution<double> dis1(lower_lim, upper_lim);
        config_data[i] = dis1(gen);
        // std::cout << "Joint " << i+1 << " = " << config_data[i] << "\n";
    }
}

static inline void
calc_offsets(const struct aa_rx_sg *sg,
             const char *fr_name_base,
             const char *fr_name_2,
             const char *fr_name_4,
             const char *fr_name_6,
             const char *fr_name_ee,
             double &d1, double &d3, double &d5, double &d7)
{
    struct aa_rx_fk *fk = aa_rx_fk_malloc(sg);
    double config_zero[7] = {0}; // 7 instead of 13, not count fixed frames. Set all to 0
    struct aa_dvec config_vec = AA_DVEC_INIT(7, config_zero, 1);
    aa_rx_fk_all(fk, &config_vec);

    /* From base to joint 2 */
    aa_rx_frame_id frame_base = aa_rx_sg_frame_id(sg, fr_name_base);
    aa_rx_frame_id frame_2 = aa_rx_sg_frame_id(sg, fr_name_2);
    double qutr_zero_data[7];
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_2, qutr_zero_data);
    d1 = qutr_zero_data[6];
    // std::cout << "Offset from base to joint 2 = " << d1 << "\n\n";

    /* From joint 2 to joint 4 */
    aa_rx_frame_id frame_4 = aa_rx_sg_frame_id(sg, fr_name_4);
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_4, qutr_zero_data);
    d3 = qutr_zero_data[6] - d1;
    // std::cout << "Offset from joint 2 to 4 = " << d3 << "\n\n";

    /* From joint 4 to joint 6 */
    aa_rx_frame_id frame_6 = aa_rx_sg_frame_id(sg, fr_name_6);
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_6, qutr_zero_data);
    d5 = qutr_zero_data[6] - d1 - d3;
    // std::cout << "Offset from joint 4 to 6 = " << d5 << "\n\n";

    /* From joint 6 to ee */
    aa_rx_frame_id frame_ee = aa_rx_sg_frame_id(sg, fr_name_ee);
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_ee, qutr_zero_data);
    d7 = qutr_zero_data[6] - d1 - d3 - d5;
    // std::cout << "Offset from joint 6 to end effector = " << d7 << "\n\n";

    aa_rx_fk_destroy(fk);
}

static inline bool
calc_SW(const double &d_SE,
        const double &d_EW,
        const double &d_WT,
        const struct amino::Vec3 &p_S,
        const struct amino::QuatTran &qutr_T,
        double &d_SW,
        double v_SW_data[3]) 
{
    /* Equation 2 of M. Gong et al. */
    double W_vec_data[] = {0, 0, -d_WT};
    double p_W_data[3];
    qutr_T.transform(W_vec_data, p_W_data);
    struct amino::Vec3 p_W{p_W_data};

    /* Check limit of d_SW */
    v_SW_data[3];
    aa_la_vsub(3, p_W.data, p_S.data, v_SW_data);
    d_SW = aa_la_norm(3, v_SW_data);
    double difference = d_SW - (d_SE + d_EW);
    if (difference >= -0.02 && difference <= 0) {
        // std::cout << "WARNING: d_SW approximately equal d_SE + d_EW --> pose is near workspace limit!\n";
    } else if (difference > 0) {
        std::cout << "ERROR: d_SW > d_SE + d_EW\n";
        // assert(d_SW - (d_SE + d_EW) <= 0);
        return false;
    }
    return true;
}

static inline void
calc_sigma(const double v_SW_data[3],
           struct amino::Quat &qu_sigma)
{
    double X_sigma_data[3],
           Y_sigma_data[3],
           Z_sigma_data[3];
    
    /* X-column */
    memcpy(X_sigma_data, v_SW_data, 3*sizeof(double));
    aa_la_normalize(3, X_sigma_data);

    /* Y-column */
    double Z_0_data[3] = {0, 0, 1};
    aa_la_cross(Z_0_data, X_sigma_data, Y_sigma_data);
    aa_la_normalize(3, Y_sigma_data);
    double temp = aa_la_norm(3, Y_sigma_data);
    if (temp <= AA_EPSILON) { // special case when X_sigma and Z_0 are parallel
        double Y_0_data[3] = {0, 1, 0};
        memcpy(Y_sigma_data, Y_0_data, 3*sizeof(double));
    }

    /* Z-column */
    aa_la_cross(X_sigma_data, Y_sigma_data, Z_sigma_data);
    aa_la_normalize(3, Z_sigma_data);

    /* Whole sigma-D coordinate system when psi = 0 */
    struct amino::RotMat R_sigma{X_sigma_data[0], Y_sigma_data[0], Z_sigma_data[0],
                                 X_sigma_data[1], Y_sigma_data[1], Z_sigma_data[1],
                                 X_sigma_data[2], Y_sigma_data[2], Z_sigma_data[2]}; // aa_tf_rotmat
    qu_sigma = qu_sigma.from_rotmat(R_sigma.data); // aa_tf_quat
    aa_tf_qminimize(qu_sigma.data);
}

static double
calc_elbow_param(const struct aa_rx_sg *sg,
                 const double config_data[7],
                 const double p_T_data[3],
                 const struct amino::Quat &qu_T,
                 const char *fr_name_base,
                 const char *fr_name_2,
                 const char *fr_name_4,
                 const char *fr_name_6,
                 const char *fr_name_ee)
{
    double elbow_ang_param = 0;

    /* Correct signs per axis-conventions of URDF file */
    double q1 = -config_data[0];
    double q2 = -config_data[1];
    double q3 = -config_data[2];
    double q4 = -config_data[3];

    /* Find offsets between frames/joints */
    double d1 = 0, 
           d3 = 0, 
           d5 = 0, 
           d7 = 0;
    calc_offsets(sg,
                 fr_name_base,
                 fr_name_2,
                 fr_name_4,
                 fr_name_6,
                 fr_name_ee,
                 d1, d3, d5, d7);

    /* Other names for the offsets */
    double d_BS = d1; // B=base(0), S=shoulder(q2), E=elbow(q4), W=wrist(q6), T=tool(q7)
    double d_SE = d3;
    double d_EW = d5;
    double d_WT = d7;

    /* Calculate vector and distance from S to W */
    struct amino::Vec3 p_S{0, 0, d_BS}; // position of shoulder 
    struct amino::Vec3 p_T{p_T_data}; // user-assigned position of tool 
    // We have orientation of tool as qu_T
    struct amino::QuatTran qutr_T{qu_T, p_T}; // user-assigned pose/transformation of tool
    aa_tf_qminimize(qutr_T.r.data);
    double v_SW_data[3] = {0};
    double d_SW = 0;
    calc_SW(d_SE,
            d_EW,
            d_WT,
            p_S,
            qutr_T,
            d_SW,
            v_SW_data);

    /* Columns for coordinate system sigma-D when psi (elbow self-motion angle) = 0 */
    struct amino::Quat qu_sigma_D_psi_0; // aa_tf_quat
    calc_sigma(v_SW_data,
               qu_sigma_D_psi_0);
    aa_tf_qminimize(qu_sigma_D_psi_0.data);

    /* Find coordinate system sigma-0 */
    double p_W_data_2[3] = {d_EW*sin(q4), 0, d_BS+d_SE+d_EW*cos(q4)};
    double v_SW_data_2[3] = {0};
    aa_la_vsub(3, p_W_data_2, p_S.data, v_SW_data_2);
    struct amino::Quat qu_sigma_0;
    calc_sigma(v_SW_data_2,
               qu_sigma_0);
    aa_tf_qminimize(qu_sigma_0.data);

    /* Calculate qu_S or R_S */
    double qu_S_given[4] = {0};
    aa_tf_eulerzyz2quat(q1, q2, q3, qu_S_given);
    aa_tf_qminimize(qu_S_given);

    /* Calculate qu_x_psi or R_x_psi */
    double qu_x_psi[4] = {0};
    double qu_temp[4] = {0};
    aa_tf_qmulnorm(qu_S_given, qu_sigma_0.data, qu_temp);
    aa_tf_qcmul(qu_sigma_D_psi_0.data, qu_temp, qu_x_psi);
    aa_tf_qnormalize(qu_x_psi);
    aa_tf_qminimize(qu_x_psi);
    
    /* Calculate psi (self-motion angle) from R_x_psi */
    double R_x_psi[9] = {0};
    aa_tf_quat2rotmat(qu_x_psi, R_x_psi);
    
    // std::cout << "\nDebug 1 - R_x_psi =\n";
    // mat_print_pretty(R_x_psi, 3, 3);

    elbow_ang_param = atan2(R_x_psi[5], R_x_psi[4]);

    return elbow_ang_param;
}

static int SCREEN_WIDTH = 800;
static int SCREEN_HEIGHT = 600;


/* ik7dof takes the DH parameter link lengths, desired position of tool,
   desired orientation of tool, and returns the joint angles needed */
std::vector<std::array<double, 7>> ik7dof(const struct aa_rx_sg *sg,
                                          const double p_T_data[3],
                                          const struct amino::Quat &qu_T,
                                          const double &elbow_ang_param,
                                          const double joint_limits[14],
                                          const char *fr_name_base,
                                          const char *fr_name_2,
                                          const char *fr_name_4,
                                          const char *fr_name_6,
                                          const char *fr_name_ee) 
{
    // std::cout << "Got into ik7dof() successfully!\n\n";

    /* Following  "Analytical Inverse Kinematics and Self-Motion 
    Application for 7-DOF Redundant Manipulator" - M. Gong et al.*/

    /* Math arm is solved as a ZYZ-Y-ZYZ arm */

    /* Initialize arrays of solutions */
    std::vector<std::array<double, 7>> sols;
    std::array<double, 7> dummy_sol;
    for (int i = 0; i < 7; i++) dummy_sol[i] = 99;
    for (int i = 0; i < 8; i++) sols.push_back(dummy_sol);

    /* Verifying input parameter of elbow self motion angle*/
    if (elbow_ang_param < -M_PI || elbow_ang_param > M_PI) {
        std::cout << "ERROR: Elbow self-motion angle parameter needs to be from -pi to pi, inclusive.\n\n";
        assert(elbow_ang_param >= -M_PI && elbow_ang_param <= M_PI);
    }

    /* Find offsets between frames/joints */
    double d1 = 0, 
           d3 = 0, 
           d5 = 0, 
           d7 = 0;
    calc_offsets(sg,
                 fr_name_base,
                 fr_name_2,
                 fr_name_4,
                 fr_name_6,
                 fr_name_ee,
                 d1, d3, d5, d7);

    /* Other names for the offsets */
    double d_BS = d1; // B=base(0), S=shoulder(q2), E=elbow(q4), W=wrist(q6), T=tool(q7)
    double d_SE = d3;
    double d_EW = d5;
    double d_WT = d7;

    /* Calculate vector and distance from S to W */
    struct amino::Vec3 p_S{0, 0, d_BS}; // position of shoulder 
    struct amino::Vec3 p_T{p_T_data}; // user-assigned position of tool 
    // We have orientation of tool as qu_T
    struct amino::QuatTran qutr_T{qu_T, p_T}; // user-assigned pose/transformation of tool
    aa_tf_qminimize(qutr_T.r.data);
    double v_SW_data[3] = {0};
    double d_SW = 0;
    bool d_SW_constraint_sat = calc_SW(d_SE,
                                       d_EW,
                                       d_WT,
                                       p_S,
                                       qutr_T,
                                       d_SW,
                                       v_SW_data);
    
    if (!d_SW_constraint_sat) {
        std::vector<std::array<double, 7>> empty_vec; // newly declared vectors are 0 size
        return empty_vec;
    }

    /* Find joint 4 */
    double cos_SEW  = (pow(d_SE,2)+pow(d_EW,2)-pow(d_SW,2)) / (2*d_SE*d_EW);
    double q4_raw = (M_PI - acos(cos_SEW));
    for (int i = 0; i < 4; i++) sols[i][3] = q4_raw;
    for (int i = 4; i < 8; i++) sols[i][3] = -q4_raw;

    /* Columns for coordinate system sigma-D when psi (elbow self-motion angle) = 0 */
    struct amino::Quat qu_sigma_D_psi_0; // aa_tf_quat
    calc_sigma(v_SW_data,
               qu_sigma_D_psi_0);
    aa_tf_qminimize(qu_sigma_D_psi_0.data);

    /* Find coordinate system sigma-D with given psi (elbow self-motion angle) */
    struct amino::XAngle x_angle_psi{elbow_ang_param}; 
    struct amino::Quat qu_x_psi{x_angle_psi}; // aa_tf_quat
    double qu_sigma_D_data[4];
    aa_tf_qmulnorm(qu_sigma_D_psi_0.data, qu_x_psi.data, qu_sigma_D_data);
    aa_tf_qminimize(qu_sigma_D_data);


    /* Find quaternion of shoulder, quaternion of wrist, and extract euler angles */
    for (int i = 0; i < 8; i++)
    {
        /* Quaternion of shoulder spherical joint*/
        // Find coordinate system sigma-0
        double q4 = sols[i][3];
        double p_W_data_2[3] = {d_EW*sin(q4), 0, d_BS+d_SE+d_EW*cos(q4)};
        double v_SW_data_2[3];
        aa_la_vsub(3, p_W_data_2, p_S.data, v_SW_data_2);
        struct amino::Quat qu_sigma_0;
        calc_sigma(v_SW_data_2,
                   qu_sigma_0);
        aa_tf_qminimize(qu_sigma_0.data);

        // Calculate quat/rotation of shoulder spherical joint
        double qu_S_data[4];
        aa_tf_qmulc(qu_sigma_D_data, qu_sigma_0.data, qu_S_data);
        aa_tf_qnormalize(qu_S_data);
        aa_tf_qminimize(qu_S_data);

        // std::cout << "qu_S_data (" << i << ") =\n";
        // array_print(qu_S_data, 4);


        /* Quaternion of wrist spherical joint*/
        struct amino::YAngle y_angle_q4{q4};
        struct amino::Quat qu_y_q4{y_angle_q4}; // aa_tf_quat
        double qu_0_4_data[4];
        aa_tf_qmulnorm(qu_S_data, qu_y_q4.data, qu_0_4_data);
        aa_tf_qminimize(qu_0_4_data);
        double qu_W_data[4];
        aa_tf_qcmul(qu_0_4_data, qu_T.data, qu_W_data); // qu_T means qu_0_7, rotation from 0 to 7
        aa_tf_qnormalize(qu_W_data);
        aa_tf_qminimize(qu_W_data);


        /* Extract euler angles */
        double euler_S[3]; // shoulder
        switch (i) {
            case 0:
            case 1:
            case 4:
            case 5:
                aa_tf_quat2eulerzyz(qu_S_data, euler_S);
                break;

            default:
                aa_tf_quat2eulerzyz_negative(qu_S_data, euler_S);
                break;
        } 
        sols[i][0] = euler_S[0]; // q1. NEED TO ADD THE SIGN CONVENTIONS HERE
        sols[i][1] = euler_S[1]; // q2
        sols[i][2] = euler_S[2]; // q3 

        double euler_W[3]; // wrist
        switch (i) {
            case 0:
            case 2:
            case 4:
            case 6:
                aa_tf_quat2eulerzyz(qu_W_data, euler_W);
                break;

            default:
                aa_tf_quat2eulerzyz_negative(qu_W_data, euler_W);
                break;
        } 
        sols[i][4] = euler_W[0]; // q5. NEED TO ADD THE SIGN CONVENTIONS HERE
        sols[i][5] = euler_W[1]; // q6
        sols[i][6] = euler_W[2]; // q7 

        /* Correct signs per axis-conventions of URDF file */
        sols[i][0] = -sols[i][0]; // q1
        sols[i][1] = -sols[i][1]; // q2
        sols[i][2] = -sols[i][2]; // q3 
        sols[i][3] = -sols[i][3]; // q4
        sols[i][4] = -sols[i][4]; // q5
        sols[i][5] = sols[i][5]; // q6
        sols[i][6] = sols[i][6]; // q7 
    }


    // std::cout << "\nAll potential solutions:\n\n";
    // for (int i = 0; i < 8; i++) {
    //     std::cout << "Solution index " << i << ":\n";
    //     for (int j = 0; j < sols[i].size(); j++) {
    //         std::cout << "q" << j+1 << " = " << sols[i][j] << "\n";
    //     } 
    //     std::cout << "\n";
    // }
    
    /* Joint limit checks */
    std::vector<int> indices_in_limit;
    for (int i = 0; i < 8; i++) {
        bool is_in_limits = true;
        is_in_limits = check_joint_limits(sols[i], joint_limits);
        if (is_in_limits) {
            indices_in_limit.push_back(i);
        }
    }

    // std::cout << "\nDebug - Number of solutions within joint limits: "\
    // << indices_in_limit.size() << "\n";

    /* Test with forward kinematics */
    struct aa_rx_fk *fk = aa_rx_fk_malloc(sg);
    aa_rx_frame_id frame_ee = aa_rx_sg_frame_id(sg, fr_name_ee);
    std::vector<int> indices_final;
    int valid_count = 0;
    for (const int& i : indices_in_limit) {
        // Given: qutr_T calculated above

        // Calculated
        double config_data_res[7] = {sols[i][0], 
                                     sols[i][1], 
                                     sols[i][2], 
                                     sols[i][3], 
                                     sols[i][4], 
                                     sols[i][5], 
                                     sols[i][6]}; 
        struct aa_dvec config_vec_res = AA_DVEC_INIT(7, config_data_res, 1);
        aa_rx_fk_all(fk, &config_vec_res);
        double qutr_T_res[7];
        aa_rx_fk_get_abs_qutr(fk, frame_ee, qutr_T_res);
        // Minimize quaternion part of qutr
        aa_tf_qminimize(qutr_T_res); // only affects first 4 elements
        aa_tf_qminimize(qutr_T.data); 
        
        // Test equal
        if ( admeq2( "1 solution configuration invalid after limit checking...",\
         qutr_T.data, qutr_T_res, AA_EPSILON, 7 ) ) {
            indices_final.push_back(i);
            valid_count++;
        }
    }

    /* Push the solutions that passed */
    std::vector<std::array<double, 7>> sols_final;
    for (const int& i : indices_final) {
        sols_final.push_back(sols[i]);
    }
    // /* Print the solutions that passed */
    // std::cout << "\nDebug 1 - Number of valid solutions (matching ee): "\
    // << sols_final.size() << "\n\n";
    
    /* Clean up allocated structures */
    aa_rx_fk_destroy(fk);
    // aa_rx_sg_destroy(sg);
    // aa_rx_win_destroy(win);

    // std::cout << "Got to the end of ik7dof()!\n\n";
    return sols_final;
}



/* -------------------------------------------------------------------------------------- */
int main(int argc, char ** argv) 
{
    std::cout << "Got into main() successfully!\n";

    /* Following  "Analytical Inverse Kinematics and Self-Motion 
    Application for 7-DOF Redundant Manipulator" - M. Gong et al.*/

    /* Create scene graph for Schunk arm */
    // Create scene graph
    struct aa_rx_sg *sg = aa_rx_sg_create();
    const char* scene_name = "schunk";
    // const char* path = "/home/billhuynh-dyalab/git/domains/build/schunk/libschunk.so";
    // const char* path = "/home/billhuynh/git/domains/build/schunk/libschunk.so";
    const char* path = "/home/billhuynh-dyalab/git/constrained-infeasibility/utils/schunk-cylinder/libschunk.so";
    // const char* path = "/home/billhuynh/git/constrained-infeasibility/utils/schunk-cylinder/libschunk.so";
    aa_rx_dl_sg_at(path, scene_name, sg, "");
    assert(sg);

    // Initialize scene graph internal structures
    int r = aa_rx_sg_init(sg); 
    assert(0 == r);

    struct aa_rx_fk *fk = aa_rx_fk_malloc(sg);

    double config_data[7] = {0}; // 7 instead of 13, not count fixed frames. Set all to 0
    struct aa_dvec config_vec = AA_DVEC_INIT(7, config_data, 1);

    /* Joint limits, here specifically for Schunk arm */
    double joint_limits[14] = {-3.12159, 3.12159, // {lower_limit_1, upper_limit_1, lower_limit_2, ...}
                               -2.12, 2.12,
                               -3.12159, 3.12159,
                               -2.16, 2.16,
                               -3.12159, 3.12159,
                               -2.07, 2.07,
                               -2.94, 2.94};

    /* User-defined frame names */
    const char *fr_name_base = "robot_podest_joint";
    const char *fr_name_2 = "robot_2_joint";
    const char *fr_name_4 = "robot_4_joint";
    const char *fr_name_6 = "robot_6_joint";
    const char *fr_name_ee = "robot_ee_joint";

    /* Randomize joints, get Td from FK, run IK, compare, repeat */    
    double total_time = 0;
    double total_percent_error = 0;
    int solved_num = 0;
    int num_runs = 1000000;
    for (int i = 0; i < num_runs; i++) {
        // std::cout << "------------ Run #" << i+1 << " of IK ------------\n\n"; 

        /* Sample 7 random joint angles */
        randomize_config(config_data, joint_limits);
        // std::cout << "\nDebug - scene's joint angles:\n";

        // config_data[0] = -1.10819; // Set to something known, for debug
        // config_data[1] = -0.370763;
        // config_data[2] = 1.85778;
        // config_data[3] = -1.45939;
        // config_data[4] = 1.40562;
        // config_data[5] = -1.19261;
        // config_data[6] = 0.973711;

        // config_data[0] = -1.5; // Set to something known, for debug
        // config_data[1] = 0;
        // config_data[2] = 0;
        // config_data[3] = M_PI/3;
        // config_data[4] = -2;
        // config_data[5] = M_PI/3;
        // config_data[6] = 2.3;

        // array_print(config_data, 7);

        aa_dvec_view(&config_vec, 7, config_data, 1);
        aa_rx_fk_all(fk, &config_vec);
        aa_rx_frame_id frame_ee = aa_rx_sg_frame_id(sg, fr_name_ee);
        double qutr_data[7];
        aa_rx_fk_get_abs_qutr(fk, frame_ee, qutr_data); // get rot and trans of ee

        // /* Visualize */ 
        // struct aa_rx_win *win = 
        //     aa_rx_win_default_create("GIVEN CONFIGURATION", SCREEN_WIDTH, SCREEN_HEIGHT);
        // aa_rx_win_set_sg(win, sg); // set the scenegraph for the window
        // aa_rx_win_set_config(win, 7, config_data);
        // aa_rx_win_run();
        // /* Clean up allocated structures */
        // aa_rx_win_destroy(win);

        /* Given position of tool (randomized) */
        double p_T_data[3];
        memcpy(p_T_data, qutr_data+4, 3*sizeof(double));

        /* Given orientation/rotation of tool (randomized) */
        double qu_T_data[4];
        memcpy(qu_T_data, qutr_data, 4*sizeof(double));
        struct amino::Quat qu_T{qu_T.from_quat(qu_T_data)};

        /* Additional parameters */
        double elbow_ang_param = calc_elbow_param(sg,
                                                  config_data,
                                                  p_T_data, 
                                                  qu_T,
                                                  fr_name_base,
                                                  fr_name_2,
                                                  fr_name_4,
                                                  fr_name_6,
                                                  fr_name_ee);
        // std::cout << "Debug - Calculated elbow_ang_param = " << elbow_ang_param << "\n\n";
    
        // elbow_ang_param = -0.5721; // from -pi to pi -- Debug

        /* Call ik function */
        std::vector<std::array<double, 7>> sols_final;
        auto beg = std::chrono::high_resolution_clock::now();
        sols_final = ik7dof(sg,
                            p_T_data, 
                            qu_T,
                            elbow_ang_param, 
                            joint_limits,
                            fr_name_base,
                            fr_name_2,
                            fr_name_4,
                            fr_name_6,
                            fr_name_ee);
        auto end = std::chrono::high_resolution_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg).count();
        total_time += duration;

        /* Print the solutions that passed */
        // std::cout << "\nDebug 2 - Number of valid solutions (matching ee): "\
        // << sols_final.size() << "\n\n";
        // for (int i = 0; i < sols_final.size(); i++) {
        //     std::cout << "Solution #" << i+1 << ":\n";
        //     for (int j = 0; j < 7; j++) {
        //         std::cout << "q" << j+1 << " = " << sols_final[i][j] << "\n";
        //     }
        //     std::cout << "\n";
        // }

        /* Display time elapsed */
        // std::cout << "Elapsed time: " << duration << " microseconds\n\n";

        /* Calculate errors compared to initial 7 joint angles started with */
        double percent_error = 0;
        int sols_count = sols_final.size();
        if (sols_count == 0) {
            percent_error = 100;
        } else {
            // int closest_sol_index = -1;
            solved_num += 1;
            double min_ssd = DBL_MAX;
            for (int i = 0; i < sols_count; i++) {
                double config_temp[7] = {0};
                for (int k = 0; k < 7; k++) { // from std_array<double> to double *
                    config_temp[k] = sols_final[i][k];
                }
                double ssd = aa_la_ssd(7, config_temp, config_data);
                // std::cout << "Debug - ssd = " << ssd << "\n";
                if (ssd < min_ssd) {
                    min_ssd = ssd;
                    // closest_sol_index = i;
                }
            }
            // std::cout << "\nDebug - min_ssd = " << min_ssd << "\n";
            // std::cout << "\nDebug - closest_sol_index = " << closest_sol_index << "\n\n";
            percent_error = sqrt(min_ssd) / aa_la_norm(7, config_data) * 100;
            // std::cout << "\nDebug - percent_error = " << percent_error << "\n";
        }

        total_percent_error += percent_error;
    }

    /* Display averages */
    std::cout << "\nAcross " << num_runs << " ik runs:\n";

    double solved_percentage = double(solved_num) / double(num_runs) * 100;
    std::cout << "Solved percentage: " << solved_percentage << " %\n";

    double average_time = total_time / num_runs;
    std::cout << "Average elapsed time: " << average_time << " microseconds\n";

    double average_error = total_percent_error / num_runs;
    std::cout << "Average percent error (vector of 7 angles): " << average_error << " %\n\n";

    /* Clean up allocated structures */
    aa_rx_fk_destroy(fk);
    aa_rx_sg_destroy(sg);
    // aa_rx_win_destroy(win);
    // // SDL_Quit();

    std::cout << "Got to the end of main()! Returning 0 next...\n\n";

    return 0;
}