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

void
aafeq(const char *name, double a, double b, double tol)
{
    if (!aa_feq(a, b, tol)) {
        fprintf(stderr, "FAILED: %s\n", name);
        fprintf(stderr, "a: %f, b: %f\n", a, b);
        abort();
    }
}

static inline void
admeq(const char *name, const double *A, const double *B, double tol, size_t size)
{
    aafeq( name, aa_la_ssd(size,A,B), 0, tol );
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
            if (fabs(x) < 1e-6) {
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

static inline void
check_joint_limits(double &q, 
                   const double joint_limits[14], 
                   const int &joint_num)
{
    int i = 0;
    while (i < 10 && q < joint_limits[(joint_num-1)*2]) {
        q += 2*M_PI;
        ++i;
    }
    i = 0;
    while (i < 10 && q > joint_limits[(joint_num-1)*2+1]) {
        q -= 2*M_PI;
        ++i;
    }
    if (q < joint_limits[(joint_num-1)*2] || q > joint_limits[(joint_num-1)*2+1]) {
        std::cout << "ERROR: Joint " << joint_num << " cannot be within its limits"; 
        std::cout << ", for the given pose and parameters.\n\n";
        assert(q >= joint_limits[(joint_num-1)*2] && q <= joint_limits[(joint_num-1)*2+1]);
    }
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
            const double &elbow_sign_param,
            const double &elbow_ang_param,
            const double &arm_sign_param,
            const double &wrist_sign_param,
            const double joint_limits[14],
            const char *fr_name_base,
            const char *fr_name_2,
            const char *fr_name_4,
            const char *fr_name_6,
            const char *fr_name_ee) 
{
    std::cout << "\nGot into ik7dof() successfully!\n\n";

    /* Following  "Analytical Inverse Kinematics and Self-Motion 
    Application for 7-DOF Redundant Manipulator" - M. Gong et al.*/

    /* Verifying input parameter of elbow joint sign */
    double abs_val = fabs(elbow_sign_param);
    if (abs_val - 1 != 0) {
        std::cout << "ERROR: Elbow joint sign parameter needs to be either -1 or 1.\n\n";
        assert(abs_val - 1 == 0);
    }

    /* Verifying input parameter of elbow self motion angle*/
    if (elbow_ang_param < -M_PI || elbow_ang_param > M_PI) {
        std::cout << "ERROR: Elbow self-motion angle parameter needs to be from -pi to pi, inclusive.\n\n";
        assert(elbow_ang_param >= -M_PI && elbow_ang_param <= M_PI);
    }

    /* Verifying input parameter of arm sign */
    double abs_val_2 = fabs(arm_sign_param);
    if (abs_val_2 - 1 != 0) {
        std::cout << "ERROR: Arm sign parameter needs to be either -1 or 1.\n\n";
        assert(abs_val_2 - 1 == 0);
    }

    /* Verifying input parameter of wrist sign */
    double abs_val_3 = fabs(wrist_sign_param);
    if (abs_val_3 - 1 != 0) {
        std::cout << "ERROR: Wrist sign parameter needs to be either -1 or 1.\n\n";
        assert(abs_val_3 - 1 == 0);
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
    std::cout << "Offset from base to joint 2 = " << d1 << "\n\n";

    // From joint 2 to joint 4
    aa_rx_frame_id frame_4 = aa_rx_sg_frame_id(sg, fr_name_4);
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_4, qutr_zero_data);
    double d3 = qutr_zero_data[6] - d1; // NOTE: offset from joint 2 to joint 4 is along joint 2's x-axis
    std::cout << "Offset from joint 2 to 4 = " << d3 << "\n\n";

    // From joint 4 to joint 6
    aa_rx_frame_id frame_6 = aa_rx_sg_frame_id(sg, fr_name_6);
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_6, qutr_zero_data);
    double d5 = qutr_zero_data[6] - d1 - d3;
    std::cout << "Offset from joint 4 to 6 = " << d5 << "\n\n";

    // From joint 6 to ee
    aa_rx_frame_id frame_ee = aa_rx_sg_frame_id(sg, fr_name_ee);
    aa_rx_fk_get_rel_qutr(fk, frame_base, frame_ee, qutr_zero_data);
    double d7 = qutr_zero_data[6] - d1 - d3 - d5;
    std::cout << "Offset from joint 6 to end effector = " << d7 << "\n\n";


    /* Other names for the offsets */
    double d_BS = d1; // B=base(0), S=shoulder(q2), E=elbow(q4), W=wrist(q6), T=tool(q7)
    double d_SE = d3;
    double d_EW = d5;
    double d_WT = d7;

    /* Position of shoulder */
    struct amino::Vec3 p_S{0, 0, d_BS};
    std::cout << "Position of S:\n";
    array_print(p_S.data, 3);

    /* User-assigned position of tool */
    struct amino::Vec3 p_T{p_T_data};
    std::cout << "Position of T:\n";
    array_print(p_T.data, 3);

    // We have orientation of tool as qu_T

    /* User-assigned pose/transformation of tool */
    struct amino::QuatTran qutr_T{qu_T, p_T};

    /* Equation 2 of M. Gong et al. */
    double W_vec_data[] = {0, 0, -d_WT};
    // double W_vec_data[] = {0, -d_WT, 0.0023}; // found out that END_EFFECTOR_GRASP Z-axis is not aligned with world Z
    double p_W_data[3];
    qutr_T.transform(W_vec_data, p_W_data);
    struct amino::Vec3 p_W{p_W_data};
    std::cout << "Position of W:\n";
    array_print(p_W.data, 3);

    // // /* DEBUG: Second way of finding vector TW */
    // // aa_rx_frame_id frame_ee = aa_rx_sg_frame_id(sg, "END_EFFECTOR_GRASP");
    // // aa_rx_frame_id W_frame = aa_rx_sg_frame_id(sg, fr_name_6);
    // // double rel_TW_qv_data[7];
    // // aa_rx_fk_get_rel_qutr(fk, frame_ee, W_frame, rel_TW_qv_data);
    // // std::cout << "Quat-trans of frame W wrt frame T:\n";
    // // array_print(rel_TW_qv_data, 7);

    // /* Check limit of d_SW */
    // double v_SW_data[3];
    // aa_la_vsub(3, p_W.data, p_S.data, v_SW_data);
    // struct amino::Vec3 v_SW{v_SW_data};
    // std::cout << "Displacement from S to W:\n";
    // array_print(v_SW.data, 3);
    // double d_SW = aa_la_norm(3, v_SW.data);
    // double difference = d_SW - (d_SE + d_EW);
    // std::cout << "Distance from S to W: " << d_SW << "\n";
    // std::cout << "d_SE + d_EW = " << d_SE + d_EW << "\n";
    // if (difference >= -0.02 && difference <= 0) {
    //     std::cout << "WARNING: d_SW approximately equal d_SE + d_EW --> pose is near workspace limit!\n";
    // } else if (difference > 0) {
    //     std::cout << "ERROR: d_SW > d_SE + d_EW\n";
    //     assert(d_SW - (d_SE + d_EW) <= 0);
    // }

    // /* Find joint 4 */
    // double cos_SEW  = (pow(d_SE,2)+pow(d_EW,2)-pow(d_SW,2)) / (2*d_SE*d_EW);
    // double q4 = elbow_sign_param * (M_PI - acos(cos_SEW));


    // /* Columns for coordinate system sigma-D when psi (elbow self-motion angle) = 0 */
    // double X_sigma_data[3];
    // memcpy(X_sigma_data, v_SW_data, 3*sizeof(double));
    // aa_la_normalize(3, X_sigma_data);
    // // std::cout << "X_sigma_data:\n";
    // // array_print(X_sigma_data, 3);

    // double Z_0_data[3] = {0, 0, 1};
    // double Y_sigma_data[3];
    // aa_la_cross(Z_0_data, X_sigma_data, Y_sigma_data);
    // aa_la_normalize(3, Y_sigma_data);
    // double temp = aa_la_norm(3, Y_sigma_data);
    // if (temp <= AA_EPSILON) { // special case when X_sigma and Z_0 are parallel
    //     double Y_0_data[3] = {0, 1, 0};
    //     memcpy(Y_sigma_data, Y_0_data, 3*sizeof(double));
    // }
    // // std::cout << "Y_sigma_data:\n";
    // // array_print(Y_sigma_data, 3);

    // double Z_sigma_data[3];
    // aa_la_cross(X_sigma_data, Y_sigma_data, Z_sigma_data);
    // aa_la_normalize(3, Z_sigma_data);
    // // std::cout << "Z_sigma_data:\n";
    // // array_print(Z_sigma_data, 3);

    // struct amino::RotMat R_sigma_D_psi_0{X_sigma_data[0], Y_sigma_data[0], Z_sigma_data[0],
    //                                      X_sigma_data[1], Y_sigma_data[1], Z_sigma_data[1],
    //                                      X_sigma_data[2], Y_sigma_data[2], Z_sigma_data[2]}; // aa_tf_rotmat
    // struct amino::Quat qu_sigma_D_psi_0{R_sigma_D_psi_0}; // aa_tf_quat


    // /* Find coordinate system sigma-D with given psi (elbow self-motion angle) */
    // /* One way */
    // // double v_SW_data_normed[3];
    // // memcpy(v_SW_data_normed, v_SW_data, 3*sizeof(double));
    // // aa_la_normalize(3, v_SW_data_normed);
    // // struct amino::AxisAngle axang_SW_psi{v_SW_data_normed, elbow_ang_param}; // aa_tf_axang
    // // struct amino::Quat qu_SW_psi{axang_SW_psi}; // aa_tf_quat
    // // double qu_sigma_D_data_1[4];
    // // aa_tf_qmulnorm(qu_SW_psi.data, qu_sigma_D_psi_0.data, qu_sigma_D_data_1); // with angle-axis, order of matmul not important
    // // std::cout << "\nChecking 1st way result:\n";
    // // array_print(qu_sigma_D_data_1, 4);
    // /* Another way */
    // struct amino::XAngle x_angle_psi{elbow_ang_param}; 
    // struct amino::Quat qu_x_psi{x_angle_psi}; // aa_tf_quat
    // double qu_sigma_D_data_2[4];
    // aa_tf_qmulnorm(qu_sigma_D_psi_0.data, qu_x_psi.data, qu_sigma_D_data_2);
    // std::cout << "\nChecking 2nd way result for sigma_D quaternion:\n";
    // array_print(qu_sigma_D_data_2, 4);


    // /* Find coordinate system sigma-0 */
    // double p_W_data_2[3] = {d_EW*sin(q4), 0, d_BS+d_SE+d_EW*cos(q4)};
    // std::cout << "Position of W when calculating sigma_0:\n";
    // array_print(p_W_data_2, 3);

    // double v_SW_data_2[3];
    // aa_la_vsub(3, p_W_data_2, p_S.data, v_SW_data_2);
    // std::cout << "Displacement from S to W when calculating sigma_0:\n";
    // array_print(v_SW_data_2, 3);

    // double X_sigma_data_2[3];
    // memcpy(X_sigma_data_2, v_SW_data_2, 3*sizeof(double));
    // aa_la_normalize(3, X_sigma_data_2);
    // // std::cout << "X_sigma_data_2:\n";
    // // array_print(X_sigma_data_2, 3);

    // // double Z_0_data[3] = {0, 0, 1};
    // double Y_sigma_data_2[3];
    // aa_la_cross(Z_0_data, X_sigma_data_2, Y_sigma_data_2);
    // aa_la_normalize(3, Y_sigma_data_2);
    // double temp_2 = aa_la_norm(3, Y_sigma_data_2);
    // if (temp_2 <= AA_EPSILON) { // special case when X_sigma and Z_0 are parallel
    //     double Y_0_data[3] = {0, 1, 0};
    //     memcpy(Y_sigma_data_2, Y_0_data, 3*sizeof(double));
    // }
    // // std::cout << "Y_sigma_data_2:\n";
    // // array_print(Y_sigma_data_2, 3);

    // double Z_sigma_data_2[3];
    // aa_la_cross(X_sigma_data_2, Y_sigma_data_2, Z_sigma_data_2);
    // aa_la_normalize(3, Z_sigma_data_2);
    // // std::cout << "Z_sigma_data_2:\n";
    // // array_print(Z_sigma_data_2, 3);

    // struct amino::RotMat R_sigma_0{X_sigma_data_2[0], Y_sigma_data_2[0], Z_sigma_data_2[0],
    //                                X_sigma_data_2[1], Y_sigma_data_2[1], Z_sigma_data_2[1],
    //                                X_sigma_data_2[2], Y_sigma_data_2[2], Z_sigma_data_2[2]}; // aa_tf_rotmat
    // struct amino::Quat qu_sigma_0{R_sigma_0}; // aa_tf_quat
    // // array_print(qu_sigma_0.data, 4);


    // /* Find rotation of shoulder spherical joint */
    // // double qu_sigma_0_inv_data[4];
    // // aa_tf_qinv(qu_sigma_0.data, qu_sigma_0_inv_data);
    // // double qu_S[4];
    // // aa_tf_qmulnorm(qu_sigma_D_data_2, qu_sigma_0_inv_data, qu_S);
    // // array_print(qu_S, 4);
    // double qu_S_data[4];
    // aa_tf_qmulc(qu_sigma_D_data_2, qu_sigma_0.data, qu_S_data);
    // // array_print(qu_S_data, 4);
    // double R_S_data[9];
    // aa_tf_quat2rotmat(qu_S_data, R_S_data);
    // std::cout << "Matrix R_S = R_0_3:\n";
    // mat_print_raw(R_S_data, 3, 3);

    // /* Calculate joint 2 */
    // double q2 = arm_sign_param * acos(R_S_data[8]);

    // /* Calculate joint 1*/
    // double q1 = atan2( arm_sign_param*R_S_data[7] , arm_sign_param*R_S_data[6] );

    // /* Calculate joint 3*/
    // double q3 = atan2( arm_sign_param*R_S_data[5] , -arm_sign_param*R_S_data[2] );


    // /* Find rotation of last three joints (5, 6, 7) combined, qu_W */
    // struct amino::YAngle y_angle_q4{q4};
    // struct amino::Quat qu_y_q4{y_angle_q4}; // aa_tf_quat
    // double qu_0_4_data[4];
    // aa_tf_qmulnorm(qu_S_data, qu_y_q4.data, qu_0_4_data);
    // /* One way*/
    // // double qu_W_data[4];
    // // aa_tf_qmulc(qu_T.data, qu_0_4_data, qu_W_data); // qu_T means qu_0_7, rotation from 0 to 7
    // // std::cout << "Quaternion of rot from 4 to 7:\n";
    // // array_print(qu_W_data, 4);
    // /* Another way, matches frame transitions better on paper */
    // double qu_W_data[4];
    // aa_tf_qcmul(qu_0_4_data, qu_T.data, qu_W_data); // qu_T means qu_0_7, rotation from 0 to 7
    // std::cout << "Quaternion of R_4_7, 2nd way:\n";
    // array_print(qu_W_data, 4);

    // double R_W_data[9];
    // aa_tf_quat2rotmat(qu_W_data, R_W_data);
    // std::cout << "Matrix R_W = R_4_7:\n";
    // mat_print_raw(R_W_data, 3, 3);

    // /* Calculate joint 6 */
    // double q6 = wrist_sign_param * acos(R_W_data[8]);

    // /* Calculate joint 5 */
    // double q5 = atan2( wrist_sign_param*R_W_data[7] , wrist_sign_param*R_W_data[6] );

    // /* Calculate joint 7 */
    // double q7 = atan2( wrist_sign_param*R_W_data[5] , -wrist_sign_param*R_W_data[2] );


    // /* Joint limit checks */
    // std::cout << "Joint q4 before limit check = " << q4 << "\n";
    // check_joint_limits(q4, joint_limits, 4);
    // std::cout << "Joint q4 after limit check = " << q4 << "\n\n";

    // std::cout << "Joint q2 before limit check = " << q2 << "\n";
    // check_joint_limits(q2, joint_limits, 2);
    // std::cout << "Joint q2 after limit check = " << q2 << "\n\n";

    // std::cout << "Joint q1 before limit check = " << q1 << "\n";
    // check_joint_limits(q1, joint_limits, 1);
    // std::cout << "Joint q1 after limit check = " << q1 << "\n\n";

    // std::cout << "Joint q3 before limit check = " << q3 << "\n";
    // check_joint_limits(q3, joint_limits, 3);
    // std::cout << "Joint q3 after limit check = " << q3 << "\n\n";

    // std::cout << "Joint q6 before limit check = " << q6 << "\n";
    // check_joint_limits(q6, joint_limits, 6);
    // std::cout << "Joint q6 after limit check = " << q6 << "\n\n";

    // std::cout << "Joint q5 before limit check = " << q5 << "\n";
    // check_joint_limits(q5, joint_limits, 5);
    // std::cout << "Joint q5 after limit check = " << q5 << "\n\n";

    // std::cout << "Joint q7 before limit check = " << q7 << "\n";
    // check_joint_limits(q7, joint_limits, 7);
    // std::cout << "Joint q7 after limit check = " << q7 << "\n\n";
    

    // /* Check with forward kinematics, using modified (proximal) DH per M. Gong et al. */
    // // Given
    // double T_0_7_given[12];
    // aa_tf_qv2tfmat(qu_T.data, p_T_data, T_0_7_given);
    // std::cout << "End effector transf matrix given:\n";
    // mat_print_raw(T_0_7_given, 3, 4);
    // // Calculated
    // double config_data_res[7] = {q1, q2, q3, q4, q5, q6, q7}; // 7 instead of 13, not count fixed frames
    // struct aa_dvec config_vec_res = AA_DVEC_INIT(7, config_data_res, 1);
    // aa_rx_fk_all(fk, &config_vec_res);
    // double qv_0_7_res[7];
    // aa_rx_fk_get_abs_qutr(fk, frame_ee, qv_0_7_res);
    // double p_0_7_res[3];
    // memcpy(p_0_7_res, qv_0_7_res+4, 3*sizeof(double));
    // double qu_0_7_res[4];
    // memcpy(qu_0_7_res, qv_0_7_res, 4*sizeof(double));
    // double T_0_7_res[12];
    // aa_tf_qv2tfmat(qu_0_7_res, p_0_7_res, T_0_7_res);
    // std::cout << "End effector transf matrix from result:\n";
    // mat_print_pretty(T_0_7_res, 3, 4);

    // // struct aa_dmat T_result = AA_DMAT_INIT(3, 4, T_0_7_res, 3);
    // // struct aa_dmat T_given = AA_DMAT_INIT(3, 4, T_0_7_given, 3);
    // admeq( "result T == given T", T_0_7_res, T_0_7_given, AA_EPSILON, 12 );


    /* Clean up allocated structures */
    aa_rx_fk_destroy(fk);
    // aa_rx_sg_destroy(sg);
    // aa_rx_win_destroy(win);

    std::cout << "Got to the end of ik7dof()!\n\n";
    return;
}






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
    for (int i = 0; i < 1; i++) {
        /* Sample 7 random joint angles */
        randomize_config(config_data, joint_limits);
        std::cout << "Scene's joint angles:\n";

        config_data[0] = 0; // Set to something known, for debug
        config_data[1] = M_PI_4;
        config_data[2] = 0;
        config_data[3] = M_PI_4;
        config_data[4] = 0;
        config_data[5] = M_PI_4;
        config_data[6] = 0;

        array_print(config_data, 7);
        aa_dvec_view(&config_vec, 7, config_data, 1);
        aa_rx_fk_all(fk, &config_vec);
        aa_rx_frame_id frame_ee = aa_rx_sg_frame_id(sg, fr_name_ee);
        double qutr_data[7];
        aa_rx_fk_get_abs_qutr(fk, frame_ee, qutr_data); // get rot and trans of ee

        /* Visualize */ 
        struct aa_rx_win *win = 
            aa_rx_win_default_create("Scenegraph win test", SCREEN_WIDTH, SCREEN_HEIGHT);
        aa_rx_win_set_sg(win, sg); // set the scenegraph for the window
        aa_rx_win_set_config(win, 7, config_data);
        aa_rx_win_run();

        /* Given position of tool (randomized) */
        double p_T_data[3];
        memcpy(p_T_data, qutr_data+4, 3*sizeof(double));

        /* Given orientation/rotation of tool (randomized) */
        double qu_T_data[4];
        memcpy(qu_T_data, qutr_data, 4*sizeof(double));
        struct amino::Quat qu_T{qu_T.from_quat(qu_T_data)};

        /* Additional parameters */
        double elbow_sign_param = 1; // either 1 or -1
        double elbow_ang_param = 0; // from -pi to pi
        double arm_sign_param = 1; // either 1 or -1
        double wrist_sign_param = 1; // either 1 or -1

        /* Call ik function */
        ik7dof(sg,
               p_T_data, 
               qu_T, 
               elbow_sign_param, 
               elbow_ang_param, 
               arm_sign_param, 
               wrist_sign_param,
               joint_limits,
               fr_name_base,
               fr_name_2,
               fr_name_4,
               fr_name_6,
               fr_name_ee);

        /* Clean up allocated structures */
        aa_rx_win_destroy(win);
    }


    /* Clean up allocated structures */
    aa_rx_fk_destroy(fk);
    aa_rx_sg_destroy(sg);
    // aa_rx_win_destroy(win);
    // // SDL_Quit();

    std::cout << "Got to the end of main()! Returning 0 next...\n\n";

    return 0;
}