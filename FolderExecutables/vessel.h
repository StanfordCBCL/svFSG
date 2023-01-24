// gnr_vessel.h
#ifndef VESSEL
#define VESSEL

using std::string;
using std::vector;
using std::cout;

class vessel {
public:
    string vessel_name;
    string file_name;
    string gnr_name;
    string equil_gnr_name;
    string exp_name;
    string hns_name;

    //Time variables
    int nts; //total number of time steps
    double dt; //time increment
    int sn; //current time step index
    double s; //actual current time

    //Geometric quantities
    double A_h, B_h, H_h; //Traction-free reference
    double A_mid_h; //Tf Midpoint reference
    double a_h, b_h, h_h; //In-vivo reference
    double a_mid_h; //In vivo midpoint
    double lambda_z_h; //Reference in vivo stretch

    vector<double> a, a_mid, b, h; //Evolving in vivo reference
    vector<double> a_pas, a_mid_pas, h_pas; //Evolving passive configuration
    vector<double> A, A_mid, B, H, lambda_z_pre; //Evolving traction free

    //Number of constituents
    int n_alpha;
    int n_pol_alpha;
    int n_native_alpha;

    //Constituent inflammatory status
    vector<int> alpha_infl;

    //Material properties
    vector<double> c_alpha_h, eta_alpha_h, g_alpha_h, G_alpha_h;
    vector<double> epsilon_pol_min;

    //Mass fractions, referential apparent densities, kinetic quantities
    vector<double> phi_alpha_h, rhoR_alpha_h, mR_alpha_h, k_alpha_h;
    vector<double> K_sigma_p_alpha_h, K_sigma_d_alpha_h, K_tauw_p_alpha_h, K_tauw_d_alpha_h;

    //True mass densities, initial volume fractions
    vector<double> rho_hat_alpha_h, epsilonR_alpha_0;

    //Reference true mass density
    double rhoR_h;

    //Histories
    vector<double> rhoR, rho, rhoR_alpha, mR_alpha, k_alpha;
    vector<double> epsilonR_alpha, epsilon_alpha;
    vector<double> ups_infl_p, ups_infl_d;

    //Reference loading quantities
    double P_h, f_h, bar_tauw_h, Q_h, P_prev;
    vector<double> sigma_h;

    //Current loading quantities
    double lambda_th_curr, lambda_z_curr;
    double P, f, bar_tauw, Q;
    vector<double> sigma, Cbar, CC, lambda_alpha_tau, lambda_z_tau;
    double mb_equil; //Current mechanobiological equil. state


    //HANSHAKE NEW VARIABLES
    double sigma_inv, sigma_inv_h;

    vector<double> F_s;
    double J_di;
    vector<double> F_curr;
    vector<double> Fi_curr;
    vector<double> lambda_act; //active radius history

    //Active stress quantities
    vector<int> alpha_active; //boolean for active constituents
    vector<double> a_act; //active radius history
    double T_act, T_act_h; //Max active stress, homeostatic max active stress
    double k_act; //Active remodeling parameters
    double lambda_0; //Min contractile stretch
    double lambda_m; //Max contractile stretch
    double CB; //Basal VC to VD ratio
    double CS; //Scaling factor for VC to VD ratio

    //Mechanobiologically equilibrated quantities
    double a_e; //equilibrated radius
    double h_e; //equilibrated thickness
    double rho_c_e; //equilibrated collagen density
    double rho_m_e; //equilibrated smc density
    double f_z_e; //equilibrated axial force
    double mb_equil_e; //Current mechanobiological equil. state

    //Phenomenologic immune parameters
    double Ki_p_h, Ki_d_h;
    double tevg_val;

    //Flags
    int num_exp_flag; //indicates whether doing reg G&R step or a numerical experiment
    int pol_only_flag; //indicates whether other constituents are produced
    int wss_calc_flag; //indicates if GnR should update its own current WSS
    int app_visc_flag; //indicates whether to use the empirical correction for viscosity from Secomb 2017

    //Initialization parameters
    //Initializing constituents
    double c1_e, c2_e, c1_m, c2_m, c1_ct, c2_ct, c1_cz, c2_cz, c1_cd1, c2_cd1, c1_cd2, c2_cd2;
    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents
    double eta_e_h, eta_m_h, eta_ct_h, eta_cz_h, eta_cd1_h, eta_cd2_h;
    //Pre-stretch parameters
    double G_e_rh, G_e_th, G_e_zh;
    double g_e_h, g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h;
    //Mass density parameters
    double rho_hat_h;
    //Homeostatic mass fractions of constituents
    double phi_e_h, phi_m_h, phi_ct_h, phi_cz_h, phi_cd1_h, phi_cd2_h;
    //Homeostatic mass densities
    double rhoR_e_h, rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h;
    //Degradation parameters
    double k_e_h, k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h;
    //Stress mediated production
    double K_sigma_p_e_h, K_sigma_p_m_h, K_sigma_p_ct_h, K_sigma_p_cz_h, K_sigma_p_cd1_h, K_sigma_p_cd2_h;
    //Stress mediated degradation
    double K_sigma_d_e_h, K_sigma_d_m_h, K_sigma_d_ct_h, K_sigma_d_cz_h, K_sigma_d_cd1_h, K_sigma_d_cd2_h;
    //Wall Shear Stress mediated production
    double K_tauw_p_e_h, K_tauw_p_m_h, K_tauw_p_ct_h, K_tauw_p_cz_h, K_tauw_p_cd1_h, K_tauw_p_cd2_h;
    //Wall Shear Stress mediated degradation
    double K_tauw_d_e_h, K_tauw_d_m_h, K_tauw_d_ct_h, K_tauw_d_cz_h, K_tauw_d_cd1_h, K_tauw_d_cd2_h;
    //Other
    //double rhoR_0 = 0.0;
    //double eta = 0;
    //double g_alpha = 0;

    //Conversions
    double um_to_m = pow(10, -6);
    double mm_to_m = pow(10, -3);
    double kPa_to_Pa = pow(10, 3);

    std::ofstream GnR_out, Equil_GnR_out, Exp_out, Mat_out, HnS_out;

    vessel(); //Default constructor
    //Vessel(string file_name); //File name constructor ***ELS USE DELEGATING CONSTRUCTOR***
    //Vessel(string file_name, string vessel_type); //File name and type constructor ***ELS USE DELEGATING CONSTRUCTOR***
    void load();
    void save();
    void print();
    void printTEVGOutputs();
    void printNativeOutputs();
    void printOutputsHandshake();
    void printNativeEquilibratedOutputs();
    void initializeNative(string native_name, double n_days_inp = 10, double dt_inp = 1);
    void initializeNativeExplicit(string native_name, double n_days_inp = 10, double dt_inp = 1);
    void initializeNativeHandshake(string native_name, double n_days_inp = 10, double dt_inp = 1);
    void initializeTEVG(string scaffold_name, string immune_name,vessel const &native_vessel, double infl_scale_trans, double n_days_inp = 10, double dt_inp = 1);
    void initializeTEVGExplicit(string scaffold_name, string immune_name,vessel const &native_vessel, double infl_scale_trans, double n_days_inp = 10, double dt_inp = 1);
    void initializeTEVGHandshake(string scaffold_name, string immune_name, double tevg_arg, double n_days_inp = 10, double dt_inp = 1);

    //TODO
    //~Vessel() destructor

};



#endif /* GNR_VESSEL */