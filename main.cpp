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

static inline bool
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

static int SCREEN_WIDTH = 800;
static int SCREEN_HEIGHT = 600;


/* ik7dof takes the DH parameter link lengths, desired position of tool,
   desired orientation of tool, and returns the joint angles needed */
void ik7dof(const struct aa_rx_sg *sg,
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
    struct aa_rx_fk *fk = aa_rx_fk_malloc(sg);
    double config_data[7] = {0}; // 7 instead of 13, not count fixed frames. Set all to 0
    struct aa_dvec config_vec = AA_DVEC_INIT(7, config_data, 1);
    aa_rx_fk_all(fk, &config_vec);

    // From base to joint 2
    aa_rx_frame_id frame_base = aa_rx_sg_frame_id(sg, fr_name_base);
    aa_rx_frame_id frame_2 = aa_rx_sg_frame_id(sg, fr_name_2);
    double qutr_zero_data[7];
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_2, qutr_zero_data);
    double d1 = qutr_zero_data[6];
    // std::cout << "Offset from base to joint 2 = " << d1 << "\n\n";

    // From joint 2 to joint 4
    aa_rx_frame_id frame_4 = aa_rx_sg_frame_id(sg, fr_name_4);
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_4, qutr_zero_data);
    double d3 = qutr_zero_data[6] - d1; // NOTE: offset from joint 2 to joint 4 is along joint 2's x-axis
    // std::cout << "Offset from joint 2 to 4 = " << d3 << "\n\n";

    // From joint 4 to joint 6
    aa_rx_frame_id frame_6 = aa_rx_sg_frame_id(sg, fr_name_6);
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_6, qutr_zero_data);
    double d5 = qutr_zero_data[6] - d1 - d3;
    // std::cout << "Offset from joint 4 to 6 = " << d5 << "\n\n";

    // From joint 6 to ee
    aa_rx_frame_id frame_ee = aa_rx_sg_frame_id(sg, fr_name_ee);
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_ee, qutr_zero_data);
    double d7 = qutr_zero_data[6] - d1 - d3 - d5;
    // std::cout << "Offset from joint 6 to end effector = " << d7 << "\n\n";


    /* Other names for the offsets */
    double d_BS = d1; // B=base(0), S=shoulder(q2), E=elbow(q4), W=wrist(q6), T=tool(q7)
    double d_SE = d3;
    double d_EW = d5;
    double d_WT = d7;

    /* Position of shoulder */
    struct amino::Vec3 p_S{0, 0, d_BS};
    // std::cout << "Position of S:\n";
    // array_print(p_S.data, 3);

    /* User-assigned position of tool */
    struct amino::Vec3 p_T{p_T_data};
    // std::cout << "Position of T:\n";
    // array_print(p_T.data, 3);

    // We have orientation of tool as qu_T

    /* User-assigned pose/transformation of tool */
    struct amino::QuatTran qutr_T{qu_T, p_T};

    /* Equation 2 of M. Gong et al. */
    double W_vec_data[] = {0, 0, -d_WT};
    // double W_vec_data[] = {0, -d_WT, 0.0023}; // found out that END_EFFECTOR_GRASP Z-axis is not aligned with world Z
    double p_W_data[3];
    qutr_T.transform(W_vec_data, p_W_data);
    struct amino::Vec3 p_W{p_W_data};
    // std::cout << "Position of W:\n";
    // array_print(p_W.data, 3);

    /* Check limit of d_SW */
    double v_SW_data[3];
    aa_la_vsub(3, p_W.data, p_S.data, v_SW_data);
    // std::cout << "Displacement from S to W:\n";
    // array_print(v_SW_data, 3);
    double d_SW = aa_la_norm(3, v_SW_data);
    double difference = d_SW - (d_SE + d_EW);
    // std::cout << "Distance from S to W: " << d_SW << "\n";
    // std::cout << "d_SE + d_EW = " << d_SE + d_EW << "\n";
    if (difference >= -0.02 && difference <= 0) {
        std::cout << "WARNING: d_SW approximately equal d_SE + d_EW --> pose is near workspace limit!\n";
    } else if (difference > 0) {
        std::cout << "ERROR: d_SW > d_SE + d_EW\n";
        assert(d_SW - (d_SE + d_EW) <= 0);
    }

    /* Find joint 4 */
    double cos_SEW  = (pow(d_SE,2)+pow(d_EW,2)-pow(d_SW,2)) / (2*d_SE*d_EW);
    double q4_raw = (M_PI - acos(cos_SEW));
    for (int i = 0; i < 4; i++) sols[i][3] = q4_raw;
    for (int i = 4; i < 8; i++) sols[i][3] = -q4_raw;


    /* Columns for coordinate system sigma-D when psi (elbow self-motion angle) = 0 */
    double X_sigma_data[3];
    memcpy(X_sigma_data, v_SW_data, 3*sizeof(double));
    aa_la_normalize(3, X_sigma_data);

    double Z_0_data[3] = {0, 0, 1};
    double Y_sigma_data[3];
    aa_la_cross(Z_0_data, X_sigma_data, Y_sigma_data);
    aa_la_normalize(3, Y_sigma_data);
    double temp = aa_la_norm(3, Y_sigma_data);
    if (temp <= AA_EPSILON) { // special case when X_sigma and Z_0 are parallel
        double Y_0_data[3] = {0, 1, 0};
        memcpy(Y_sigma_data, Y_0_data, 3*sizeof(double));
    }

    double Z_sigma_data[3];
    aa_la_cross(X_sigma_data, Y_sigma_data, Z_sigma_data);
    aa_la_normalize(3, Z_sigma_data);

    struct amino::RotMat R_sigma_D_psi_0{X_sigma_data[0], Y_sigma_data[0], Z_sigma_data[0],
                                         X_sigma_data[1], Y_sigma_data[1], Z_sigma_data[1],
                                         X_sigma_data[2], Y_sigma_data[2], Z_sigma_data[2]}; // aa_tf_rotmat
    struct amino::Quat qu_sigma_D_psi_0{R_sigma_D_psi_0}; // aa_tf_quat
    aa_tf_qminimize(qu_sigma_D_psi_0.data);


    /* Find coordinate system sigma-D with given psi (elbow self-motion angle) */
    struct amino::XAngle x_angle_psi{elbow_ang_param}; 
    struct amino::Quat qu_x_psi{x_angle_psi}; // aa_tf_quat
    double qu_sigma_D_data[4];
    aa_tf_qmulnorm(qu_sigma_D_psi_0.data, qu_x_psi.data, qu_sigma_D_data);
    aa_tf_qminimize(qu_sigma_D_data);

    /* Find quaternion of should, quaternion of wrist, and extract euler angles */
    for (int i = 0; i < 8; i++)
    {
        /* Quaternion of shoulder spherical joint*/

        // Find coordinate system sigma-0
        double q4 = sols[i][3];
        double p_W_data_2[3] = {d_EW*sin(q4), 0, d_BS+d_SE+d_EW*cos(q4)};
        double v_SW_data_2[3];
        aa_la_vsub(3, p_W_data_2, p_S.data, v_SW_data_2);  
        double X_sigma_0_data[3];
        memcpy(X_sigma_0_data, v_SW_data_2, 3*sizeof(double));
        aa_la_normalize(3, X_sigma_0_data);

        double Y_sigma_0_data[3];
        aa_la_cross(Z_0_data, X_sigma_0_data, Y_sigma_0_data);
        aa_la_normalize(3, Y_sigma_0_data);
        double temp_2 = aa_la_norm(3, Y_sigma_0_data);
        if (temp_2 <= AA_EPSILON) { // special case when X_sigma and Z_0 are parallel
            double Y_0_data[3] = {0, 1, 0};
            memcpy(Y_sigma_0_data, Y_0_data, 3*sizeof(double));
        }

        double Z_sigma_0_data[3];
        aa_la_cross(X_sigma_0_data, Y_sigma_0_data, Z_sigma_0_data);
        aa_la_normalize(3, Z_sigma_0_data);

        struct amino::RotMat R_sigma_0{X_sigma_0_data[0], Y_sigma_0_data[0], Z_sigma_0_data[0],
                                       X_sigma_0_data[1], Y_sigma_0_data[1], Z_sigma_0_data[1],
                                       X_sigma_0_data[2], Y_sigma_0_data[2], Z_sigma_0_data[2]}; // aa_tf_rotmat
        struct amino::Quat qu_sigma_0{R_sigma_0}; // aa_tf_quat
        aa_tf_qminimize(qu_sigma_0.data);

        // Calculate quat/rotation of shoulder spherical joint
        double qu_S_data[4];
        aa_tf_qmulc(qu_sigma_D_data, qu_sigma_0.data, qu_S_data);
        aa_tf_qnormalize(qu_S_data);
        aa_tf_qminimize(qu_S_data);


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


    std::cout << "\nAll potential solutions:\n\n";
    for (int i = 0; i < 8; i++) {
        std::cout << "Solution index " << i << ":\n";
        for (int j = 0; j < sols[i].size(); j++) {
            std::cout << "q" << j+1 << " = " << sols[i][j] << "\n";
        } 
        std::cout << "\n";
    }
    
    /* Joint limit checks */
    std::vector<int> indices_in_limit;
    for (int i = 0; i < 8; i++) {
        bool is_in_limits = true;
        is_in_limits = check_joint_limits(sols[i], joint_limits);
        if (is_in_limits) {
            indices_in_limit.push_back(i);
        }
    }

    std::cout << "\nNumber of solutions within joint limits: " << indices_in_limit.size() << "\n";

    /* Test with forward kinematics */
    // TODO: TEST WITH QUAT-TRANS INSTEAD OF CONVERTING TO T_MATRIX[12]
    std::vector<int> indices_final;
    for (const int& i : indices_in_limit) {
        // std::cout << "\nTesting solution set with index " << i << "...\n";
        
        // Given
        double T_0_7_given[12];
        aa_tf_qv2tfmat(qu_T.data, p_T_data, T_0_7_given);
        // std::cout << "\nEnd effector transf matrix given:\n";
        // mat_print_pretty(T_0_7_given, 3, 4);
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
        double p_0_7_res[3];
        memcpy(p_0_7_res, qutr_T_res+4, 3*sizeof(double));
        double qu_0_7_res[4];
        memcpy(qu_0_7_res, qutr_T_res, 4*sizeof(double));
        double T_0_7_res[12];
        aa_tf_qv2tfmat(qu_0_7_res, p_0_7_res, T_0_7_res);
        // std::cout << "\nEnd effector transf matrix from result:\n";
        // mat_print_pretty(T_0_7_res, 3, 4);
        // Test equal
        if ( admeq2( "result T == given T", T_0_7_res, T_0_7_given, AA_EPSILON, 12 ) ) {
            indices_final.push_back(i);
        }
    }


    /* Print the solutions that passed */
    std::cout << "\nValid solutions:\n\n";
    for (const int& i : indices_final) {
        std::cout << "Solution index " << i << ":\n";
        for (int j = 0; j < sols[i].size(); j++) {
            std::cout << "q" << j+1 << " = " << sols[i][j] << "\n";
        } 
        std::cout << "\n";
    }

    
    /* Clean up allocated structures */
    aa_rx_fk_destroy(fk);
    // aa_rx_sg_destroy(sg);
    // aa_rx_win_destroy(win);

    // std::cout << "Got to the end of ik7dof()!\n\n";
    return;
}



/* -------------------------------------------------------------------------------------- */
int main(int argc, char ** argv) 
{
    std::cout << "Got into main() successfully!\n\n";

    /* Following  "Analytical Inverse Kinematics and Self-Motion 
    Application for 7-DOF Redundant Manipulator" - M. Gong et al.*/

    /* Create scene graph for Schunk arm */
    // Create scene graph
    struct aa_rx_sg *sg = aa_rx_sg_create();
    const char* scene_name = "schunk";
    // const char* path = "/home/billhuynh-dyalab/git/domains/build/schunk/libschunk.so";
    // const char* path = "/home/billhuynh/git/domains/build/schunk/libschunk.so";
    // const char* path = "/home/billhuynh-dyalab/git/constrained-infeasibility/utils/schunk-cylinder/libschunk.so";
    const char* path = "/home/billhuynh/git/constrained-infeasibility/utils/schunk-cylinder/libschunk.so";
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
    for (int i = 0; i < 1; i++) {
        // std::cout << "------------ Run #" << i+1 << " of IK ------------\n\n"; 

        /* Sample 7 random joint angles */
        randomize_config(config_data, joint_limits);
        std::cout << "Scene's joint angles:\n";

        config_data[0] = -1.10819; // Set to something known, for debug
        config_data[1] = -0.370763;
        config_data[2] = 1.85778;
        config_data[3] = -1.45939;
        config_data[4] = 1.40562;
        config_data[5] = -1.19261;
        config_data[6] = 0.973711;

        // config_data[0] = -1.5; // Set to something known, for debug
        // config_data[1] = 0;
        // config_data[2] = 0;
        // config_data[3] = M_PI/3;
        // config_data[4] = -2;
        // config_data[5] = M_PI/3;
        // config_data[6] = 2.3;

        array_print(config_data, 7);
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
        // double elbow_ang_param = 0; // from -pi to pi
        double elbow_ang_param = -0.5721; // from -pi to pi

        /* Call ik function */
        ik7dof(sg,
               p_T_data, 
               qu_T,
               elbow_ang_param, 
               joint_limits,
               fr_name_base,
               fr_name_2,
               fr_name_4,
               fr_name_6,
               fr_name_ee);

        // Debug
        std::cout << "Scene's joint angles, shown again:\n";
        array_print(config_data, 7);
    }


    /* Clean up allocated structures */
    aa_rx_fk_destroy(fk);
    aa_rx_sg_destroy(sg);
    // aa_rx_win_destroy(win);
    // // SDL_Quit();

    std::cout << "Got to the end of main()! Returning 0 next...\n\n";

    return 0;
}