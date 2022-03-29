// Library includes
//#include <string>

//Define global variables and physical constants
//Include all needed libraries
#include "resources.hh"

//Inclusion of function declaration batches
#include "functions.hh"

//Function for cleaning up file levels of 'type', -1 for all.
void clean_up(std::vector<std::vector<std::string> > matrix_created_files, int type);
template <typename T> std::string Num2Str ( T Number );
template <typename T> bool EmptyOption(T x);
template <typename T> double calc_RAM( T data ) {
    std::size_t S;
    S = sizeof data;
    return ((double)S)*1E-6; //MB
}

int main (int argc, char *argv[]) {

	//Option struct
	OPTIONS opt;
	std::string input_delim;

	//Simulation data struct
	SIMULATION sim;

	//Dedicated data structure itterators
	unsigned int i,j;

	//Look at the time!
	char *time_char;
    time_t rawtime;
    struct tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    time_char = asctime(timeinfo);

    sim.start_t = clock();
    
    clock_t EXEC_TIMES_tick;

	//Execution path
	std::string selfpath = get_selfpath();
	selfpath.erase(selfpath.end()-6,selfpath.end());

	//###########
	//CREATE LIST OF ALL CREATED FILES
	// - FOR OPENING, CLOSING, COPY,
	// - DELETE, CREATE, and CLEANUP
	//###########

	//FOLDERS
	std::string output = "output/";
	std::string debug = "debug/";
	std::string input = "input/";

	std::vector<std::string> list_created_files;
	std::vector<std::vector<std::string> > matrix_created_files;

	//On new program launch clear levels 0-2 to prepare clean simulation
	//On resume of old simulation, use status file in output data and clean level 0
	//On simulation completion, collect levels 1-2 into new unique folder and clear level 0

	//Create row 0: Temporary files
	list_created_files.push_back(selfpath + output + "TEMPORARY_DATA.tmp"); //0
	matrix_created_files.push_back(list_created_files);
	list_created_files.clear();

	//Create row 1: Debug data
	list_created_files.push_back(selfpath + output + debug + "system_energy_ang_momentum.dat"); //0
	list_created_files.push_back(selfpath + output + debug + "execution_times.dat"); //1
	matrix_created_files.push_back(list_created_files);
	list_created_files.clear();

	//Create row 2: Output data
	list_created_files.push_back(selfpath + output + "close_encounters.dat"); //0
	list_created_files.push_back(selfpath + output + "removals.dat");//1
	list_created_files.push_back(selfpath + output + "msg_log.txt");//2
	list_created_files.push_back(selfpath + output + "summary.txt");//3
	list_created_files.push_back(selfpath + output + "simulation_resume.cfg");//4
	list_created_files.push_back(selfpath + output + "state_data.dat");//5
    list_created_files.push_back(selfpath + output + "snapshot_data.dat");//6
	matrix_created_files.push_back(list_created_files);
	list_created_files.clear();

	//Initial state column format:
	//x,y,z,vx,vy,vz (position in m, speeds in m/s i.e. SI units)

	//Properties column format:
    //Object type (1 = massive, 0 = test particle w/ mass, 2 = test particle w/o mass)
    //Mass

	//Create row 3: Input files
	list_created_files.push_back(selfpath + input + "initial_state.dat");//0
	list_created_files.push_back(selfpath + input + "properties.dat");//1
	list_created_files.push_back(selfpath + input + "GHOPE.cfg");//2
    list_created_files.push_back(selfpath + input + "snapshots.dat");//3
	matrix_created_files.push_back(list_created_files);
	list_created_files.clear();

	// ######### LOAD EXTERNAL CONFIG ##########
    std::cout << "Loading config from " << matrix_created_files[3][2] << std::endl;
    load_config(matrix_created_files[3][2],&opt);

    std::cout << std::setprecision(10);

    //REMOVE LEFT OVER FILES
	clean_up(matrix_created_files,0);
	clean_up(matrix_created_files,1);
	clean_up(matrix_created_files,2);

	//Logfile out
    std::ofstream  log_out;
	// INIT OUTPUT
	log_out.open(matrix_created_files[2][2].c_str(), std::ios::app);

	output_stream out(std::cout, log_out);
	out.type = 1;

	out << "Cleaned up files from levels 0 trough 2 to prepare for new simulation";

	out << "############################################";
	out << "Running G-HOPE Version: " + Num2Str(GHOPE_version);
	out << "Copyright @ Daniel Kastinen [The Swedish Institute of Space Physics]";
	out << "Credit is to be supplied according to GitLab wiki instructions";
	out << "############################################";


	// ######### Create data structures ##########

	//Temporary data 
    tensor<double> temp_tensor;
    tensor<double> temp_tensor1;
    tensor<double> temp_tensor2;
    std::vector<double> temp_vec;
    std::vector<double> temp_vec2;
    std::vector<std::vector<double> > temp_mat;
    std::vector<double>::const_iterator d_vec_itter_start;
	std::vector<double>::const_iterator d_vec_itter_stop;

    //support_data_structure
    tensor<double> interaction_data;

    std::vector<double> snapshot_list_time;

    //Object properties, column format:
    //0 body id
    //1 Object type (-1 = removed, 1 = massive, 0 = test particle w/ mass, 2 = test particle w/o mass)
    //2 Mass (blank or 0 for massless)
    //3 Density (optional)
    std::vector<std::vector<double> > property_matrix;
    tensor<double> m_vector;
    tensor<double> rho_vector;
    tensor<double> radius_vector;
    tensor<int> type_vector;
    tensor<int> body_id_vector;

    //Vectors for inter-step storing before file i/o
    std::vector<std::vector<double> > save_mat_time_diag;
    std::vector<std::vector<double> > save_mat_energy_diag;


    tensor<double> save_mat_state;
    tensor<double> snapshot_data;

    //Phase space data
    int DIMS[3];
    tensor<double> q_vec;
    tensor<double> p_vec;
    tensor<double> p_vec_temp;
    tensor<double> q_vec_temp;

    tensor<double> Q_vec;
    tensor<double> P_vec;
    

    //Time and itteration tracking
    std::vector<double> T_vec;
    std::vector<double> EXEC_TIMES;


    // ######### LOAD INITIAL CONDITION DATA ##########
	
    input_delim = (char)opt.delim_int;

    //get IC data
    if(file_exists(matrix_created_files[3][0])) {
        temp_mat.clear();
        out << "Loading data file from: " << matrix_created_files[3][0];
        load_data_matrix(&temp_mat,matrix_created_files[3][0],input_delim);
        

        DIMS[0] = temp_mat.size();
        DIMS[1] = 3;
        q_vec.resize(DIMS,2);
        p_vec.resize(DIMS,2);

        for(i=0; i<temp_mat.size();i++) {
            //Re-formating of state data goes here
            //E.g if initial state is not in SI units
            for(j=0; j < 3; j++) {
                q_vec(i,j) = temp_mat[i][j]*AU;
                p_vec(i,j) = temp_mat[i][j+3];
            }

        }
    }
    else {
    	throw("No phase space initial conditions file");
    }

    //Get prop data
    if(file_exists(matrix_created_files[3][1])) {
        property_matrix.clear();
        out << "Loading data file from: " << matrix_created_files[3][1];
        load_data_matrix(&property_matrix,matrix_created_files[3][1],input_delim);
    }
    else {
    	throw("No object properties file");
    }

    //get snapshot list file if enabled
    if(opt.snapshot_list) {
        if(file_exists(matrix_created_files[3][3])) {
            temp_mat.clear();
            out << "Loading data file from: " << matrix_created_files[3][3];
            load_data_matrix(&temp_mat,matrix_created_files[3][3],input_delim);
        }
        else {
            throw("Snapshot list enabled but no file found");
        }

        for(i=0; i < temp_mat.size(); i++) {
            snapshot_list_time.push_back(temp_mat[i][0]*day2s);
            if(opt.snapshot_list && (temp_mat[i][0] < opt.start_time || temp_mat[i][0] > opt.end_time)) {
                throw("Snapshots cannot be taken outside integration interval");
            }
        }
        std::sort(snapshot_list_time.begin(),snapshot_list_time.end()); //sort accoring to order of calculation
    }

    out << "" << "Extracting data into smaller structures.";

    DIMS[0] = property_matrix.size();
    body_id_vector.resize(DIMS,1);
    type_vector.resize(DIMS,1);
    m_vector.assign(DIMS,1,0.0);
    rho_vector.assign(DIMS,1,0.0);
    radius_vector.assign(DIMS,1,0.0);

    for(i=0; i < property_matrix.size(); i++) {
        body_id_vector(i) = (int)property_matrix[i][0];
    	type_vector(i) = (int)property_matrix[i][1];

        //If user has input mass and density add them, else set to zero
        if(property_matrix[i].size() == 3) {
            m_vector(i) = property_matrix[i][2];
            if(property_matrix[i].size() >= 4) {
                rho_vector(i) = property_matrix[i][3]; 
                radius_vector(i) = pow(3.0*property_matrix[i][2]/(4*PI*property_matrix[i][3]),2.0/3.0);
            }
        }
    }
	
    if(property_matrix.size() != q_vec.size(0)) {
    	throw("Amount of bodies in initial conditions file does not match amount in properties file");
    }
    else {
    	sim.body_n = q_vec.size(0);
    }

    out << "" << "Mapping velocity to momentum.";
    //Input is in velocity, phase space uses momentum
    //Map V -> P
    //Also count body distribution
    sim.interacting_body_n = 0;
    sim.noninteracting_body_n = 0;
    for(i=0; i < sim.body_n; i++) {
    	if(type_vector(i) == 1) {
    		sim.interacting_body_n++;
    	}
    	else {
    		sim.noninteracting_body_n++;
    	}
    	p_vec(i,0) *= m_vector(i);
    	p_vec(i,1) *= m_vector(i);
    	p_vec(i,2) *= m_vector(i);
    }

    out << "" << "Containing data sizes: ";
    out       << "  - Position   " + Num2Str(q_vec.size(0)) + " x " + Num2Str(q_vec.size(1)) + " = " + Num2Str(q_vec.size());
    out       << "  - Momentum   " + Num2Str(p_vec.size(0)) + " x " + Num2Str(p_vec.size(1)) + " = " + Num2Str(p_vec.size());
    out       << "  - Masses     " + Num2Str(m_vector.size(0)) + " = " + Num2Str(m_vector.size());

    if(opt.verbose_level == 2) {
        out << "" << "INITIAL STATE OF INTERACTING BODIES: ";
        for(i=0; i < sim.body_n; i++) {
            if(type_vector[i] == 1) {
                temp_tensor1 = q_vec.col(i);
                temp_tensor2 = p_vec.col(i);
                out << "ID: " + Num2Str(body_id_vector[i]);
                out << "Position: " + Num2Str(q_vec(i,0)/AU) + " AU, " + Num2Str(q_vec(i,1)/AU)  + " AU, " + Num2Str(q_vec(i,1)/AU) + " AU. Norm = " + Num2Str(abs_v(temp_tensor1)/AU) + " AU";
                out << "Velocity: " + Num2Str(p_vec(i,0)/m_vector(i)*0.001) + " km/s, " + Num2Str(p_vec(i,1)/m_vector(i)*0.001) + " km/s, " + Num2Str(p_vec(i,1)/m_vector(i)*0.001) + " km/s. Norm = " + Num2Str(abs_v(temp_tensor2)/m_vector(i)*0.001) + " km/s";
            }
        }
        out << "";
    }

    //Find central body index in data struct
    if(EmptyOption(opt.central_body)) {
        if(opt.verbose_level > 0) {
            out << "No central body selected, finding the most massive";
        }
        sim.central_body_index = (unsigned int)(max_element(m_vector.begin(),m_vector.end()) - m_vector.begin());
    }
    else {
        for(i=0; i < sim.body_n; i++) {
            if(body_id_vector(i) == opt.central_body) {
                sim.central_body_index = i;
            }
        }
    }

    out << "";
    if(opt.snapshot_list && snapshot_list_time.size() > 0) {
        out << "Non empty snapshot list found with option enabled, printing list: ";
        for(i=0; i < snapshot_list_time.size(); i++) {
            out << "Snapshot " + Num2Str(i+1) + ": " + Num2Str(snapshot_list_time[i]*s2year) + " y";
        }
    }
    out << "";
    

    if(opt.verbose_level > 0) {
        out << "Bodies counted. Total number          : " + Num2Str(sim.body_n);
        out << "Number of interacting bodies          : " + Num2Str(sim.interacting_body_n);
        out << "Number of non-interacting bodies      : " + Num2Str(sim.noninteracting_body_n);
        out << "Central body located in data as index : " + Num2Str(sim.central_body_index);
        out << "                              with ID : " + Num2Str(body_id_vector[sim.central_body_index]) << "";
    }

    //Additional Re-formating of property data goes here
    //E.g if mass is not in SI units


    // ########################################
    // # Integrator preparation
    // ########################################
    // ###########
    // # Calculating all simulation parameters from input
    // ###########
    if(EmptyOption(opt.steps)) {

        if(EmptyOption(opt.dt))
            {throw("Time step cannot be empty if loops option is empty");}
        else if(EmptyOption(opt.start_time))
            {throw("Start time cannot be empty if loops option is empty");}
        else if(EmptyOption(opt.end_time))
            {throw("End time cannot be empty if loops option is empty");}

        //dt, start and end time in days
        opt.end_time = opt.end_time*day2s;
        opt.start_time = opt.start_time*day2s;
        opt.dt = opt.dt*day2s;
        opt.steps = (unsigned int)floor((opt.end_time - opt.start_time)/opt.dt);
    }
    else if(EmptyOption(opt.dt)) {

        if(EmptyOption(opt.steps))
            {throw("Loops cannot be empty if time step option is empty");}
        else if(EmptyOption(opt.start_time))
            {throw("Start time cannot be empty if time step option is empty");}
        else if(EmptyOption(opt.end_time))
            {throw("End time cannot be empty if time step option is empty");}

        //start and end time in days
        opt.end_time = opt.end_time*day2s;
        opt.start_time = opt.start_time*day2s;

        opt.dt = (opt.end_time - opt.start_time)/((double)opt.steps);

    }
    else if(EmptyOption(opt.start_time)) {

        if(EmptyOption(opt.steps))
            {throw("Loops cannot be empty if start time option is empty");}
        else if(EmptyOption(opt.dt))
            {throw("Time step cannot be empty if start time option is empty");}
        else if(EmptyOption(opt.end_time))
            {throw("End time cannot be empty if start time option is empty");}

        //start and end time in days
        opt.end_time = opt.end_time*day2s;
        opt.dt = opt.dt*day2s;
        opt.start_time = opt.end_time - opt.dt*((double)opt.steps);
    }
    else if(EmptyOption(opt.end_time)) {

        if(EmptyOption(opt.steps))
            {throw("Loops cannot be empty if end time option is empty");}
        else if(EmptyOption(opt.dt))
            {throw("Time step cannot be empty if end time option is empty");}
        else if(EmptyOption(opt.start_time))
            {throw("Start time cannot be empty if end time option is empty");}

        //start and end time in days
        opt.start_time = opt.start_time*day2s;
        opt.dt = opt.dt*day2s;
        opt.end_time = opt.start_time + opt.dt*((double)opt.steps);
    }

    out << "" << "### Stepper settings configured to ";
    out <<       "    Time step               : " + Num2Str(opt.dt) + " s";
    out <<       "                            : " + Num2Str(opt.dt*s2day) + " d";
    out <<       "    Step count              : " + Num2Str(opt.steps);
    out <<       "    Start time              : " + Num2Str(opt.start_time*s2year) + " y";
    out <<       "    End time                : " + Num2Str(opt.end_time*s2year) + " y" << "";

    if(EmptyOption(opt.save_interval_itter) && EmptyOption(opt.save_inteval_time)) {
        sim.save_number = 1;
        sim.save_interval = opt.steps;
        opt.save_inteval_time = opt.end_time - opt.start_time;
    }
    else if(EmptyOption(opt.save_interval_itter)) {
        opt.save_inteval_time = opt.save_inteval_time*day2s; //Input in days
        sim.save_interval = opt.save_inteval_time/opt.dt; // auto floor from d -> uint
        sim.save_number = opt.steps/sim.save_interval; // auto floor uint/uint
    }
    else if(EmptyOption(opt.save_inteval_time)) {
        sim.save_number = (unsigned int)(floor(opt.steps/opt.save_interval_itter));
        sim.save_interval = opt.save_interval_itter;
    }
    else {
        out << "Both save interval for itterations and time were given, using the most dense option.";
        opt.save_inteval_time = opt.save_inteval_time*day2s; //Input in days

        //less itterations per save = more dense
        if(floor(opt.save_inteval_time/opt.dt) > opt.save_interval_itter) {
            sim.save_number = (unsigned int)(floor(opt.steps/opt.save_interval_itter));
            sim.save_interval = opt.save_interval_itter;
        }
        else {
            sim.save_interval = opt.save_inteval_time/opt.dt; // auto floor from d -> uint
            sim.save_number = opt.steps/sim.save_interval; // auto floor uint/uint 
        }
    }

    //if we are trying to save at a larger interval then the integration
    //then just save at the end
    //OR
    //if we are trying to save at an interval less then the step
    //then save at each step
    if(sim.save_interval > opt.steps) {sim.save_interval = opt.steps;}
    else if(sim.save_interval == 0) { sim.save_interval = 1; }

    //As we must save according to itterative steps, change the output time to correspond to true save intervals
    opt.save_inteval_time = (double)sim.save_interval*(opt.dt);

    double align_T;
    if(opt.dt*((double)sim.save_interval)*((double)sim.save_number) < (opt.end_time - opt.start_time) ) {
        sim.last_align = 1;
        align_T =  (opt.end_time - opt.start_time) - opt.dt*((double)sim.save_interval)*((double)sim.save_number);
    }
    else {
        sim.last_align = 0;
        align_T = 0;
    }

    out << "" << "### Save settings configured to ";
    out <<       "    Save interval itterations : " + Num2Str(sim.save_interval);
    out <<       "    Save interval time        : " + Num2Str(opt.save_inteval_time*s2year) + " y, OR";
    out <<       "                              : " + Num2Str(opt.save_inteval_time*s2day) + " d";
    out <<       "    Save number               : " + Num2Str(sim.save_number + sim.last_align);
    out <<       "    Last time aligning step   : " + Num2Str(sim.last_align);
    out <<       " Regular steps " + Num2Str(sim.save_number) + "*" + Num2Str(sim.save_interval) + "=" + Num2Str(sim.save_number*sim.save_interval) + " at step: " + Num2Str(opt.dt*s2day) + " d";
    out <<       " Aligning steps " + Num2Str(sim.last_align) + "*" + Num2Str(sim.save_interval) + "=" + Num2Str(sim.last_align*sim.save_interval) + " at step: " + Num2Str(align_T*s2day/((double)sim.save_interval)) + " d";
    out <<       " Snapshots NOT INCLUDED, these are counted separatly.";
    out << "";

    // ###########
    // # Declaring variables and structures
    // ###########
    //Substep-list used by StepperBSfixed
	std::vector<unsigned int> N_substep;
	N_substep.push_back(2);
	N_substep.push_back(4);
	N_substep.push_back(6);
    N_substep.push_back(8);
    N_substep.push_back(10);
    N_substep.push_back(12);
    N_substep.push_back(14);
    N_substep.push_back(16);

    //Vars used by StepperBS
    const Int nvar = sim.body_n*6;
    //Need to find correct atol and rtol
    const Doub Atol = opt.ATOL, rtol = opt.RTOL, h1 = opt.dt, hmin = 0.0, x1 = opt.start_time, x2 = opt.start_time + opt.dt*((double)sim.save_interval);
    VecDoub ystart(nvar);
    VecDoub Mdoub(sim.body_n);
    VecDoub rdoub(sim.body_n);
    VecDoub Typedoub(sim.body_n);

    //Optimized 2 Struct used by StepperBS
    //Avrage step execution             : 0.000319562 s
    //at rtol = atol = 1e-3
    //~13% speedup
    //ADD INLINE FUNCTIONS
    //Avrage step execution             : 0.000294936 s
    //~20% speedup (+7%)
    //Changed vec<vec> > to MatDoub
    //Avrage step execution             : 0.000281218 s
    //~24% speedup (+4%)
    struct rhs_grav {
        VecDoub M,Type;
        MatDoub R;

        rhs_grav(VecDoub mM, VecDoub TypeP) : M(mM), Type(TypeP), R(mM.size(),mM.size()) {}
        void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
            int N,Nt,D,i,j,I,J,O,I_part;
            Doub R3;
            D = 3;
            Nt = y.size();
            N = Nt/2;
            O = N/D;
            //Set the dq/dt
            for(i=0;i < N; i++) {
                I = i/D;// auto floor
                dydx[i] = y[i+N]/M[I];
            }

            //calc help matrix of inter-object distances for speedup
            for(I=0; I < O-1; I++) {
                i = I*3; 
                for(J=I+1;J < O; J++) {
                    j = J*3;
                    R[I][J] = sqrt( inline_SQR(y[i] - y[j]) + inline_SQR(y[i+1] - y[j+1]) + inline_SQR(y[i+2] - y[j+2]));
                }
            }

            //Set the dp/dt
            for(i=N;i < Nt; i++) {
                I = (i-N)/D;// auto floor
                I_part = (i-N) - I*D;
                dydx[i] = 0.0;
                for(J=0;J < O; J++) {
                    j = J*D + I_part;
                    if(Type[J] == 1 && I != J) {
                        if(I < J) {
                            R3 = inline_CUB(R[I][J]);
                        }
                        else {
                            R3 = inline_CUB(R[J][I]);
                        }
                        dydx[i] = dydx[i] + M[J]*(y[j] - y[i-N])/R3;
                    }
                }
                dydx[i] = M[I]*G*dydx[i];
            }
        }
    }; 


    //Optimized Struct used by StepperBS with PR-effect
    struct rhs_grav_PR {
        VecDoub M,Type,r;
        int SOL;
        bool PR_on;
        rhs_grav_PR(VecDoub mM, VecDoub TypeP, VecDoub rr, int ssol, bool PR) : M(mM), Type(TypeP), r(rr), SOL(ssol),PR_on(PR) {}
        void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
            int N,Nt,D,i,j,I,J,O,I_part,sol;
            Doub R3,Qpr,Rsol,vdotx;
            Doub vx,vy,vz,xx,yy,zz,EMcoef;
            Qpr = 1.0;
            D = 3;
            Nt = y.size();
            N = Nt/2;
            O = N/D;
            //Set the dq/dt
            for(i=0;i < N; i++) {
                I = i/D;
                dydx[i] = y[i+N]/M[I];
            }

            //calc help matrix of inter-object distances for speedup
            std::vector<std::vector<double> > R;
            R.resize(O);
            for(I=0; I < O-1; I++) {
                R[I].resize(O);
                i = I*3;
                for(J=I+1;J < O; J++) {
                    j = J*3;
                    R[I][J] = sqrt( (y[i] - y[j])*(y[i] - y[j]) + (y[i+1] - y[j+1])*(y[i+1] - y[j+1]) + (y[i+2] - y[j+2])*(y[i+2] - y[j+2]));
                }
            }

            //Set the dp/dt
            for(i=N;i < Nt; i++) {
                I = (i-N)/D;
                I_part = (i-N) - I*D;
                dydx[i] = 0.0;
                for(J=0;J < O; J++) {
                    j = J*D + I_part;
                    if(Type[J] == 1 && I != J) {
                        if(I < J) {
                            R3 = R[I][J]*R[I][J]*R[I][J];
                        }
                        else {
                            R3 = R[J][I]*R[J][I]*R[J][I];
                        }
                        dydx[i] = dydx[i] + M[J]*(y[j] - y[i-N])/R3;
                    }
                }
                dydx[i] = M[I]*G*dydx[i];
            }

            if(PR_on) {
                //PR calculation, additive to dp/dt
                for(J=0;J < O; J++) {
                    j = J*D;
                    sol = SOL*D;
                    if(Type[J] != 1) {
                        if(SOL < J) {
                            Rsol = R[SOL][J];
                        }
                        else {
                            Rsol = R[J][SOL];
                        }
                        vx = (y[j+N  ] - y[sol+N  ]); xx = (y[j  ] - y[sol  ]);
                        vy = (y[j+N+1] - y[sol+N+1]); yy = (y[j+1] - y[sol+1]);
                        vz = (y[j+N+2] - y[sol+N+2]); zz = (y[j+2] - y[sol+2]);

                        EMcoef = (r[J]*r[J]*L0*Qpr)/(4.0*Rsol*Rsol*c0*M[J]);
                        vdotx = vx*xx + vy*yy + vz*zz;

                        dydx[j  ] = dydx[j  ] + EMcoef*((1.0 - vdotx/(Rsol*c0))*xx/Rsol - vx/c0 );
                        dydx[j+1] = dydx[j+1] + EMcoef*((1.0 - vdotx/(Rsol*c0))*yy/Rsol - vy/c0 );
                        dydx[j+2] = dydx[j+2] + EMcoef*((1.0 - vdotx/(Rsol*c0))*zz/Rsol - vz/c0 );
                    }
                }
            }

        }
    }; 

    // ###########
    // # Setup of the data structures for the integrator
    // # Applying coordinate transformations for the integrator
    // ###########

    std::vector<double> integrator_coef;

    switch(opt.integrator) {
    	case 0:
    		//Only fill the data structures if we are going to use the StepperBS algorithm
			for(i=0; i< sim.body_n; i++) {
				ystart[3*i  ]              = q_vec(i,0);
				ystart[3*i+1]              = q_vec(i,1);
				ystart[3*i+2]              = q_vec(i,1);
				ystart[3*(i+sim.body_n)  ] = p_vec(i,0);
				ystart[3*(i+sim.body_n)+1] = p_vec(i,1);
				ystart[3*(i+sim.body_n)+2] = p_vec(i,1);
				Mdoub[i] = m_vector(i);
				Typedoub[i] = type_vector(i);
                rdoub[i] = radius_vector(i);
			}
    		out << "Integrator setup for adaptive Bulirsch–Stoer stepper algorithm (Numerical recipes 3) complete.";
		   	break;
        case 1:

            DIMS[0] = sim.body_n;
            DIMS[1] = sim.body_n;
            DIMS[2] = 3;
            interaction_data.resize(DIMS,3);

            integrator_coef.push_back(0.5);
            integrator_coef.push_back(1.0);
            integrator_coef.push_back(0.5);

            out << "Integrator setup for Kinetic/Potential drift kick drift Hamiltonian split complete.";
            break;
    }

    // ###########
    // # Setup all needed structs to use the avalible integrators
    // ###########

    //Setup the system and stepper for StepperBS
    Output out_ODE;
    //rhs_grav_PR SYSTEM(Mdoub,Typedoub,rdoub,sim.central_body_index,opt.PR_effect);
    rhs_grav SYSTEM(Mdoub,Typedoub);
    Odeint<StepperBS<rhs_grav> > SYSTEM_ode(ystart,x1,x2,Atol,rtol,h1,hmin,out_ODE,SYSTEM);

    // ######### Preparation for integration ##########
    double T_current;
    double t_step,next_T_target;
    bool next_is_snapshot;

    T_current = opt.start_time;
    T_vec.push_back(T_current);

    // ###########
    // # Synching test-particles (Adaptive-BS)
    // ###########


    // #### Save inital energy and total angular momentum if option on
    if(opt.save_energy_diag) {
        temp_tensor = Total_angular_momentum(q_vec,p_vec,m_vector,type_vector);
        temp_vec.clear();
        temp_vec.push_back(0.0);
        temp_vec.push_back(T_vec.back());
        temp_vec.push_back(temp_tensor(0));
        temp_vec.push_back(temp_tensor(1));
        temp_vec.push_back(temp_tensor(2));
        temp_vec.push_back(abs_v(temp_tensor));
        temp_vec.push_back(Total_energy(q_vec,p_vec,m_vector,type_vector));
        save_mat_energy_diag.push_back(temp_vec);
    }


    // ######### CHECKING RAM REQ ##########
    double state_RAM_usage = 0;
    state_RAM_usage += calc_RAM(q_vec);
    state_RAM_usage += calc_RAM(p_vec);
    state_RAM_usage += calc_RAM(type_vector);
    state_RAM_usage += calc_RAM(m_vector);
    state_RAM_usage += calc_RAM(body_id_vector);
    state_RAM_usage += calc_RAM(save_mat_energy_diag);
    sim.save_buildup = (unsigned int)floor( opt.max_RAM/state_RAM_usage );

    out << "" << "One state data set calculated to use " + Num2Str(state_RAM_usage) + " Mb RAM.";
    out << "Buffering max " + Num2Str(sim.save_buildup) + " states before HDD write to fully use RAM" << "";

	// ######### ---------------- ##########
    // ######### INTEGRATION LOOP ##########
    // ######### ---------------- ##########

    sim.t_elapsed = 0;
    sim.t_estimate = 0;
    sim.save_counter = 0;
    sim.console_print_c = 0;
    sim.nok = 0;
    sim.nbad = 0;
    sim.snapshot_counter = 0;

    //saving initial state
    DIMS[0] = 10; //10 columns (row size)
    DIMS[1] = sim.body_n; //column size (number of rows)
    DIMS[2] = 1; //Depth (number of matrix slices)
    snapshot_data.resize(DIMS,3);

    DIMS[2] = sim.save_buildup >= (sim.save_number + sim.last_align+1) ? (sim.save_number + sim.last_align+1) : sim.save_buildup;
    out << "Reiszing save state RAM storage data structure to: " + Num2Str(DIMS[0]) + "x" + Num2Str(DIMS[1]) + "x" + Num2Str(DIMS[2]);
    save_mat_state.resize(DIMS,3);

    switch(opt.coord) {
        case 0:
            out << "No coordinate transforms to be performed";
            break;
        case 1:
            out << "Heliocentric Earth Ecliptic Vernal Equinox chosen, performing transformation every state save.";
            out << "Assuming Earth has index 3 in data structure. THIS NEEDS FIXING";
            map2centric(q_vec,p_vec,m_vector,sim.central_body_index);
            map2ecliptic(q_vec,p_vec,3);

            break;
    }

    merge_data(&save_mat_state,0,&q_vec,&p_vec,&m_vector,&type_vector,&body_id_vector,T_vec.back());

    out << "" << "Starting integration loop";

    double LOAD_OVERHEAD_TIME = (double)(clock() - sim.start_t)/CLOCKS_PER_SEC;
    sim.steps_performed = 0;
    sim.loops_performed = 0;
	for(i=0; i < sim.save_number + sim.last_align; i++) {
        EXEC_TIMES.clear();
        EXEC_TIMES_tick = clock();
        sim.loops_performed++;

		if(sim.steps_performed == 0 || ((T_current - opt.start_time)/(opt.end_time - opt.start_time) > ((double)sim.console_print_c + 1)/100.0) && opt.verbose_level > 0) {
			sim.console_print_c += floor((T_current - opt.start_time)/(opt.end_time - opt.start_time)*100.0) - sim.console_print_c;
            out << "Progress " + Num2Str(sim.console_print_c+1) + 
             "%| Starting step " +  Num2Str(sim.steps_performed+1) + "/" + Num2Str(opt.steps) +
             ", " + Num2Str(T_current*s2year) + " y out of" +
             " " + Num2Str((double)opt.steps*opt.dt*s2year) + " y" +
             " | Time elapsed " + Num2Str(sim.t_elapsed) + 
             " s, Estimated time left " + Num2Str(sim.t_estimate) + " s";
		}

        //See if next itteration should be a snapshot itteration
        if(opt.snapshot_list && sim.snapshot_counter < snapshot_list_time.size() ) {
            next_T_target = T_current + opt.dt*((double)sim.save_interval);
            if(next_T_target > snapshot_list_time[sim.snapshot_counter]) {
                next_is_snapshot = true;
                out << "Performing snapshot " + Num2Str(sim.snapshot_counter) + ": ";
                out << "Propagating from " + Num2Str(T_current*s2year) + " y to " + Num2Str(snapshot_list_time[sim.snapshot_counter]*s2year) + " y";
            }
            else {
                next_is_snapshot = false;
            }
        }
        else {
            next_is_snapshot = false;
        }
        
        //TIME DIAG: 0 CONSOLE MSG
        EXEC_TIMES.push_back((double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC);
        EXEC_TIMES_tick = clock();
        switch(opt.integrator) {
        	case 0: //Adaptive BS-alg
                if(i==sim.save_number) { // Only happends if sim.last_align == 1 and lopp on that step
                    SYSTEM_ode.x2 = opt.end_time;
                }
                else if(next_is_snapshot) { // Do a sub-step to the snapshot and save it, then continue as usual
                    next_T_target = SYSTEM_ode.x2; // Save next target
                    SYSTEM_ode.x2 = snapshot_list_time[sim.snapshot_counter]; // change next target to snapshot
                    SYSTEM_ode.integrate(); //propagate to snapshot
                    //transfer data
                    for(j=0; j< sim.body_n; j++) {
                        q_vec(j,0) = ystart[3*j  ];
                        q_vec(j,1) = ystart[3*j+1];
                        q_vec(j,2) = ystart[3*j+2];
                        p_vec(j,0) = ystart[3*(j+sim.body_n)  ];
                        p_vec(j,1) = ystart[3*(j+sim.body_n)+1];
                        p_vec(j,2) = ystart[3*(j+sim.body_n)+2];
                    }
                    // create snapshot and save
                    out << "Saving snapshot to file";
                    switch(opt.coord) {
                        case 1:
                            map2centric(q_vec,p_vec,m_vector,sim.central_body_index);
                            map2ecliptic(q_vec,p_vec,3);
                            break;
                    }
                    merge_data(&snapshot_data,0,&q_vec,&p_vec,&m_vector,&type_vector,&body_id_vector,snapshot_list_time[sim.snapshot_counter]);
                    save_tensor(matrix_created_files[2][6], &snapshot_data, opt.precision, input_delim); 
                    
                    //set next integration to range from snapshot to next target
                    SYSTEM_ode.x1 = SYSTEM_ode.x2;
                    SYSTEM_ode.x2 = next_T_target;

                    //snapshot done, save progress
                    out << "Resuming normal operations, propagating to " + Num2Str(next_T_target*s2year) + " y";
                    sim.snapshot_counter++;
                }
                
                SYSTEM_ode.integrate();
                //TIME DIAG: 1 INTEGRATION STEP
                EXEC_TIMES.push_back((double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC);
                EXEC_TIMES_tick = clock();

                SYSTEM_ode.x1 = SYSTEM_ode.x2;
                SYSTEM_ode.x2 = SYSTEM_ode.x2 + opt.dt*((double)sim.save_interval);

                T_current = SYSTEM_ode.x;
                T_vec.push_back(T_current);


                //Transfer data into standard structs
                for(j=0; j< sim.body_n; j++) {
                    q_vec(j,0) = ystart[3*j  ];
                    q_vec(j,1) = ystart[3*j+1];
                    q_vec(j,2) = ystart[3*j+2];
                    p_vec(j,0) = ystart[3*(j+sim.body_n)  ];
                    p_vec(j,1) = ystart[3*(j+sim.body_n)+1];
                    p_vec(j,2) = ystart[3*(j+sim.body_n)+2];
                }

                sim.nok = SYSTEM_ode.nok;
                sim.nbad = SYSTEM_ode.nbad;
                sim.steps_performed = SYSTEM_ode.nok;

                //TIME DIAG: 2 DATA TRANSFER
                EXEC_TIMES.push_back((double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC);
                EXEC_TIMES_tick = clock();
	        	break;
            case 1: //Kin-Pot DKD split
                if(i==sim.save_number) { // Only happends if sim.last_align == 1 and lopp on that step
                    t_step = (opt.end_time - T_current)/((double)sim.save_interval);
                    out << "Last aligning step: switching to step " + Num2Str(t_step*s2day) + " d from " + Num2Str(opt.dt*s2day) + " d";
                }
                else if(next_is_snapshot) { // Do a sub-step to the snapshot and save it, then continue as usual
                    //lower time step (to ensure accuracy)
                    //propagate to snapshot
                    t_step = (snapshot_list_time[sim.snapshot_counter] - T_current)/((double)sim.save_interval); // change next target to snapshot
                    for(j=0; j < sim.save_interval; j++) {
                        H_kin(q_vec,p_vec,m_vector,type_vector,t_step*integrator_coef[0]);
                        H_pot(q_vec,p_vec,m_vector,type_vector,t_step*integrator_coef[1],interaction_data);
                        H_kin(q_vec,p_vec,m_vector,type_vector,t_step*integrator_coef[2]);
                    }
                    sim.steps_performed += sim.save_interval;
                    T_current += t_step*(double)sim.save_interval;
   
                    // create snapshot and save
                    out << "Saving snapshot to file";
                    switch(opt.coord) {
                        case 1:
                            map2centric(q_vec,p_vec,m_vector,sim.central_body_index);
                            map2ecliptic(q_vec,p_vec,3);
                            break;
                    }
                    merge_data( &snapshot_data,
                        0,
                        &q_vec,
                        &p_vec,
                        &m_vector,
                        &type_vector,
                        &body_id_vector,
                        snapshot_list_time[sim.snapshot_counter]);
                    save_tensor(matrix_created_files[2][6], &snapshot_data, opt.precision, input_delim); 
                    
                    //set next integration to range from snapshot to next target
                    
                    t_step = ( next_T_target - snapshot_list_time[sim.snapshot_counter])/((double)sim.save_interval);
                    //snapshot done, save progress
                    out << "Resuming normal operations, propagating to " + Num2Str(next_T_target*s2year) + " y";
                    sim.snapshot_counter++;
                }
                else {
                    t_step = opt.dt;
                }

                for(j=0; j < sim.save_interval; j++) {
                        H_kin(q_vec,p_vec,m_vector,type_vector,t_step*integrator_coef[0]);
                        H_pot(q_vec,p_vec,m_vector,type_vector,t_step*integrator_coef[1],interaction_data);
                        H_kin(q_vec,p_vec,m_vector,type_vector,t_step*integrator_coef[2]);
                }
                //TIME DIAG: 1 INTEGRATION STEP
                EXEC_TIMES.push_back((double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC);
                EXEC_TIMES_tick = clock();
                
                sim.steps_performed += sim.save_interval;
                
                sim.nok += sim.save_interval;
                sim.nbad += 0;

                T_current += t_step*(double)sim.save_interval;
                T_vec.push_back(T_current);

                //TIME DIAG: 2 DATA TRANSFER
                EXEC_TIMES.push_back((double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC);
                EXEC_TIMES_tick = clock();
                break;
            case 2: //non-adaptive BS alg
                
                break;

	    }

        sim.save_counter++;
        sim.current_save_buildup = sim.save_counter % sim.save_buildup;
        sim.t_elapsed = (double)(clock() - sim.start_t)/CLOCKS_PER_SEC;
        sim.t_estimate = (sim.t_elapsed/((double)sim.loops_performed))*( (double)sim.save_number - (double)sim.loops_performed );

        //Save energy and total angular momentum if option on
        if(opt.save_energy_diag) {
            temp_tensor = Total_angular_momentum(q_vec,p_vec,m_vector,type_vector);
            temp_vec.clear();
            temp_vec.push_back((double)sim.loops_performed);
            temp_vec.push_back(T_current);
            temp_vec.push_back(temp_tensor(0));
            temp_vec.push_back(temp_tensor(1));
            temp_vec.push_back(temp_tensor(2));
            temp_vec.push_back(abs_v(temp_tensor));
            temp_vec.push_back(Total_energy(q_vec,p_vec,m_vector,type_vector));
            save_mat_energy_diag.push_back(temp_vec);
            //Save to file when RAM full or on last itteration
            if((sim.save_buildup % sim.save_counter == 0) || ( i == (sim.save_number + sim.last_align - 1) ) ) {
                save_mat(matrix_created_files[1][0], &save_mat_energy_diag, opt.precision, input_delim); 
                save_mat_energy_diag.clear();
            }
        }
        
        switch(opt.coord) {
            case 1:
                map2centric(q_vec,p_vec,m_vector,sim.central_body_index);
                map2ecliptic(q_vec,p_vec,3);
                break;
        }
        merge_data( &save_mat_state,
                    sim.current_save_buildup,
                    &q_vec,
                    &p_vec,
                    &m_vector,
                    &type_vector,
                    &body_id_vector,
                    T_current);


        //Save to file when RAM full or on last itteration
        if((sim.current_save_buildup == 0) || ( i == (sim.save_number + sim.last_align - 1) ) ) {
            save_tensor(matrix_created_files[2][5], &save_mat_state, opt.precision, input_delim); 
        }

        //TIME DIAG: 3 STATE SAVES
        EXEC_TIMES.push_back((double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC);
        //Save execution time if option on
        if(opt.save_time_diag) {
            temp_vec.clear();
            temp_vec.push_back((double)sim.loops_performed);
            temp_vec.push_back((double)sim.steps_performed);
            temp_vec.insert(temp_vec.end(),EXEC_TIMES.begin(),EXEC_TIMES.end());
            temp_vec.push_back(sim.t_elapsed);//TIME DIAG: 4 TIME ELAPSED
            save_mat_time_diag.push_back(temp_vec);
            //Save to file when RAM full or on last itteration
            if((sim.save_buildup % sim.save_counter == 0) || ( i == (sim.save_number + sim.last_align - 1) ) ) {
                save_mat(matrix_created_files[1][1], &save_mat_time_diag, opt.precision, input_delim);
                save_mat_time_diag.clear();
            }
            
        }

	}

    out << "Cleanning up simulation temporary files";
	//Simulation done, remove level 0
	clean_up(matrix_created_files,0);

	sim.end_t = clock();
	sim.runtime = (double)(sim.end_t - sim.start_t)/CLOCKS_PER_SEC;

    out << "################ - DONE - ###################";
    out << "Simulation complete, summary data below";

    out << "Number of steps performed         : " + Num2Str(sim.steps_performed);
    out << "Number of proposed steps          : " + Num2Str(opt.steps);
    out << "Step adaptation                   : " + Num2Str((int)sim.steps_performed-(int)opt.steps);
    out << "Number of successful steps        : " + Num2Str(sim.nok);
    out << "Number of failed steps            : " + Num2Str(sim.nbad);
    out << "Number of loops                   : " + Num2Str(sim.loops_performed);
    out << "Avrage steps per loop             : " + Num2Str((double)sim.steps_performed/(double)sim.loops_performed);
    out << "Number of saves                   : " + Num2Str(sim.save_counter);
    out << "Number of snapshots               : " + Num2Str(sim.snapshot_counter);
    out << "Simulation runtime                : " + Num2Str(sim.runtime) + " s";
    out << "Loading time                      : " + Num2Str(LOAD_OVERHEAD_TIME) + " s";
    out << "Avrage step execution             : " + Num2Str(sim.runtime/(double)sim.steps_performed) + " s";
    out << "Simulation end time               : " + Num2Str(T_current*s2year) + " y, OR";
    out << "                                  : " + Num2Str(T_current*s2day) + " d";
    out << "Simulation proposed end time      : " + Num2Str(opt.end_time*s2year) + " y, OR";
    out << "                                  : " + Num2Str(opt.end_time*s2day) + " d";


	return 0;
}



void clean_up(std::vector<std::vector<std::string> > matrix_created_files, int type) {
    unsigned int i,j;
    bool selector;
    if(type == -1) {
        selector = TRUE;
    }
    else {
        selector = FALSE;
    }
    for(j=0; j < matrix_created_files.size(); j++ ) {
        if((type == j) || selector) {
            for (i = 0; i < matrix_created_files[j].size(); ++i) {
                if(file_exists(matrix_created_files[j][i])) {
                    remove(matrix_created_files[j][i].c_str());
                    std::cout << "Removed file: " << matrix_created_files[j][i] << std::endl;
                }
            }
        }
    }
}

template <typename T> 
bool EmptyOption(T x) {
    return x == (T)std::numeric_limits<unsigned int>::max();
}
