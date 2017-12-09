#include <iostream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <stdexcept>

#include "FVSearcher.hpp"

using namespace std;

// Info of help
void show_helpinfo() {
    cout <<"A Parallel Free Volume Searcher"<<endl<<endl
    	 <<"Usage: ./fvsearch sstructure_file region_of_interest [OPTIONS]"<<endl
    	 <<"\t-h,-H,--help\tShow help informations"<<endl
    	 <<"\t-o <string> \tOutput file name"<<endl
    	 <<"\t-n <int>    \tSet number of threads to use"<<endl
    	 <<"\t-s <double> \tSpacing for output mesh-grid"<<endl
    	 <<"\t-l <double> \tSpacing for grouping nearby atoms"<<endl
    	 <<"\t-m <int>    \tSize of parallel blocks"<<endl
    	 //<<"\t-vs <int>   \tDANGER"<<endl
    	 <<"\t-test       \tTest mode"<<endl;
}

// Info of invalid arguments
void error_invalid_arguments() {
	cerr <<"[Runtime Error] Not enough or invalid arguments."
	     <<"Please refer to help by running with option -h"<<endl;
}

// Execute the free volume searcher
void execute(FVSearcher fvsearcher, bool test_flag) {

	double t_start, t_end;
	string info;

	// -----------------------------------------------
	// STEP 0 Reading input files
	// -----------------------------------------------
    {
    	// Time record start
	    t_start = omp_get_wtime();
	    cout <<"INFO: Reading input files."<<endl; 

	    // Read input
	    try{
	    	fvsearcher.read_inputs();
	    }catch (runtime_error &e) {
	    	cerr <<"[Runtime Error] "<<e.what()
	    		 <<"Please refer to help by running with option -h"<<endl;
	    	return ;
	    }

		// Show Structure and ROI information
		ROI roi = fvsearcher.get_roi();
		int num_atoms = fvsearcher.get_num_atoms();
	    cout <<"INFO:  - ROI:"<<"x:("<<roi.xb<<","<<roi.xe
	    	 <<"), y:("<<roi.yb<<","<<roi.ye
	    	 <<"), z:("<<roi.zb<<","<<roi.ze<<")"<<endl;
	    cout <<"INFO:  - Number of atoms: "<<num_atoms<<endl;

	    // Time record end
	    info = "Read_inputs";
	    t_end = omp_get_wtime();
	    cout <<"TIME:  - "<<info<<" "<<fixed << std::setprecision(3)<<t_end - t_start<<" seconds"<<endl;
	}

	// -----------------------------------------------
	// STEP 1 Generate array of occupancies
	// -----------------------------------------------
	{
		cout <<"INFO: Step 1 Generate array of occupancies"<<endl;

		int full_sparse[3] = {0,0,1};
		int point_atom[3] = {0,1,0};
		string infos[3] = {"Full_pointwise", "Full_atomwise", "Sparse"};
		string full_infos[3] = {"Building full array of occupancies by pointwise iteration", 
								"Building full array of occupancies by atomwise iteration",
								"Building sparse array of occupancies"};

		if (test_flag) {
			// Three implementations
			for (int i = 0; i < 3; i++) {
				// Time record start
			    t_start = omp_get_wtime();
			    cout <<"INFO: Method "<<(i+1)<<" - "<<full_infos[i]<<endl;

			    // Calculate
			   	fvsearcher.build_array_occupancy(full_sparse[i], point_atom[i]);

			    // Time record end
			    info = infos[i];
			    t_end = omp_get_wtime();
			    //  Show info
			    if (point_atom[i] == 0){
				    cout <<"INFO:  - Maximal van der Waals radius: "<<fvsearcher.get_r_max()<<" Angstrom"<<endl;
			        cout <<"INFO:  - Table of nearby atoms has "<<fvsearcher.get_T_size()<<" blocks"<<endl;
			    }
			    int count_atoms_in_ROI = fvsearcher.get_num_atoms_in_ROI();
			    if (count_atoms_in_ROI <= 1) cout <<"INFO:  - "<<count_atoms_in_ROI<<" atom in ROI"<<endl;
			    else cout <<"INFO:  - "<<count_atoms_in_ROI<<" atoms in ROI"<<endl;
			    cout <<"INFO:  - "<<info<<" validation number "
			    	 <<fvsearcher.validate_array_occupancy(full_sparse[i])<<endl;
			    cout <<"TIME:  - "<<info<<" "<<fixed << std::setprecision(3)<<t_end - t_start<<" seconds"<<endl;
			}
			// Compare memory usage
			cout <<"INFO:  - Full array of occupancy "<<fixed << std::setprecision(3)
				 <<1.0 * sizeof(bool) * fvsearcher.get_M_full_size() / 1024 / 1024<<" MB"<<endl;
			cout <<"INFO:  - Sparse array of occupancy "<<fixed << std::setprecision(3)
			     <<1.0 * sizeof(int) * fvsearcher.validate_array_occupancy(1) / 1024 / 1024<<" MB"<<endl;   
		} else {
			// Default using full_pointwise
			// Time record start
			int i = 0;
		    t_start = omp_get_wtime();
		    cout <<"INFO: Default method - "<<full_infos[i]<<endl;
		    
		    // Calculate
		   	fvsearcher.build_array_occupancy(full_sparse[i], point_atom[i]);

		    // Time record end
		    info = infos[i];
		    t_end = omp_get_wtime();
		    int count_atoms_in_ROI = fvsearcher.get_num_atoms_in_ROI();
		    if (count_atoms_in_ROI <= 1) cout <<"INFO:  - "<<count_atoms_in_ROI<<" atom in ROI"<<endl;
		    else cout <<"INFO:  - "<<count_atoms_in_ROI<<" atoms in ROI"<<endl;
		    cout <<"INFO:  - "<<info<<" validation number "
		    	 <<fvsearcher.validate_array_occupancy(full_sparse[i])<<endl;
		    cout <<"TIME:  - "<<info<<" "<<fixed << std::setprecision(3)<<t_end - t_start<<" seconds"<<endl;
		}
	    

	}

	// -----------------------------------------------
	// STEP 2 Parallel connected-components labeling
	// -----------------------------------------------
	{
		// 2.1 Serial two-pass CCL for each block
		// Time record start
	    t_start = omp_get_wtime();
	    cout <<"INFO: STEP 2 Parallel connected-components labeling"<<endl;
	    cout <<"INFO:  - 2.1 Serial two-pass CCL for each block"<<endl;

	    // Calculate
		fvsearcher.parallel_CCL_serial_block();

		// Time record end
		info = "Serial_CCL";
	 	t_end = omp_get_wtime();
		cout <<"INFO:  - "<<fvsearcher.get_num_labels_provisional()<<" provisional labels assigned"<<endl;
		cout <<"TIME:  - "<<info<<" "<<fixed << std::setprecision(3)<<t_end - t_start<<" seconds"<<endl;
	}

	{
		// 2.2 Retrieve global equivalencies
		// Time record start
	    t_start = omp_get_wtime();
	    cout <<"INFO:  - 2.2 Retrieve global equivalencies"<<endl;

	    // Calculate
		fvsearcher.parallel_CCL_retrieve_equivalencies();

		// Time record end
		info = "Retri_equiv";
	 	t_end = omp_get_wtime();
		cout <<"INFO:  - Totally "<<fvsearcher.get_num_global_equivalencies()<<" global equivalencies found"<<endl;
		cout <<"INFO:  - After applying equivalencies, finally there are "
			 <<fvsearcher.get_num_labels_global()<<" global labels"<<endl;
		cout <<"TIME:  - "<<info<<" "<<fixed << std::setprecision(3)<<t_end - t_start<<" seconds"<<endl;
	}

	{
		// 2.3 Apply global equivalencies
		// Time record start
	    t_start = omp_get_wtime();
	    cout <<"INFO:  - 2.3 Apply global equivalencies"<<endl;

	    // Calculate
		fvsearcher.parallel_CCL_apply_equivalencies();

		// Time record end
		info = "Apply_equiv";
	 	t_end = omp_get_wtime();
	 	cout <<"TIME:  - "<<info<<" "<<fixed << std::setprecision(3)<<t_end - t_start<<" seconds"<<endl;
	}
	// -----------------------------------------------
	// STEP 4 Output
	// -----------------------------------------------
	if (fvsearcher.arg.output_filename != "__NO_OUTPUT__"){
		// Time record start
	    t_start = omp_get_wtime();
	    cout <<"INFO:  STEP 4 Output"<<endl;

	    // Calculate
		fvsearcher.output();

		// Time record end
		info = "Output";
	 	t_end = omp_get_wtime();
	 	cout <<"INFO:  - Output to file "<<fvsearcher.get_output_filename()<<endl;
	 	cout <<"TIME:  - "<<info<<" "<<fixed << std::setprecision(3)<<t_end - t_start<<" seconds"<<endl;
	}
}

int main(int argc, char ** argv) {
	bool test_flag = false;
    // Read arguments from input
    if (argc < 3) {
    	if (argc == 2 and  // Show help information
    		(string(argv[1]) == "-h" or \
    		 string(argv[1]) == "-H" or \
    		 string(argv[1]) == "--help")) {
			show_helpinfo();
			return 0;
		}else{ // Invalid arguments
			error_invalid_arguments();
			return -1;
		}
    }else{
    	// Free volume searcher
    	FVSearcher fvsearcher;
    	FVSArguments arg;

	    arg.structure_filename = argv[1]; // Structural file input name
	    arg.roi_filename = argv[2]; // ROI input file name
	    for (int i = 3; i < argc; i++) {
	    	bool error_flag = false;
	    	if (string(argv[i]) == "-n") { // Number of threads in parallel
	    		if (i + 1 < argc) arg.num_threads = atoi(argv[i+1]);
	    		else error_flag = true;
	    		i++;
	    	}else if (string(argv[i]) == "-s") { // Spacing for mesh grid
	    		if (i + 1 < argc) arg.spacing_grid = atof(argv[i+1]);
	    		else error_flag = true;
	    		i++;
	    	}else if (string(argv[i]) == "-l") { // Spacing for grouping nearby atoms
	    		if (i + 1 < argc) arg.spacing_group = atof(argv[i+1]);
	    		else error_flag = true;
	    		i++;
	    	}else if (string(argv[i]) == "-m") { // Size of parallel blocks
	    		if (i + 1 < argc) arg.size_pblock = atoi(argv[i+1]);
	    		else error_flag = true;
	    		i++;
	    	}else if (string(argv[i]) == "-o") { // Size of parallel blocks
	    		if (i + 1 < argc) arg.output_filename = argv[i+1];
	    		else error_flag = true;
	    		i++;
	    	}else if (string(argv[i]) == "-vs") { // Size of parallel blocks
	    		if (i + 1 < argc) arg.SCALING_FACTOR = atof(argv[i+1]);
	    		else error_flag = true;
	    		i++;
	    	}else if (string(argv[i]) == "-test") { // Size of parallel blocks
	    		test_flag = true;
	    	}
	    	if (error_flag) {
	    		error_invalid_arguments();
				return -1;
	    	}
	    }

	    // Check arguments
	    if (arg.num_threads < 1) { // Number of threads in parallel cannot be zero or negative
	    	cerr <<"[Runtime Error] Invalid number of threads " << arg.num_threads
	     		 <<" Please refer to help by running with option -h"<<endl;
	     	return -1;
	    }
	    if (arg.spacing_grid <= 0) {
	    	cerr <<"[Runtime Error] Invalid spacing for mesh grid " << arg.spacing_grid
	     		 <<" Please refer to help by running with option -h"<<endl;
	     	return -1;
	    }
	    if (arg.spacing_group <= 0) {
	    	cerr <<"[Runtime Error] Invalid spacing for grouping nearby atoms " << arg.spacing_group
	     		 <<" Please refer to help by running with option -h"<<endl;
	     	return -1;
	    }
	    if (arg.size_pblock < 1) {
	    	cerr <<"[Runtime Error] Invalid size of parallel blocks " << arg.size_pblock
	     		 <<" Please refer to help by running with option -h"<<endl;
	     	return -1;
	    }

	    // Show arguments readed
	    cout <<"INFO: Free Volume Searcher begins. "<<endl
	         <<"INFO:  - Structure filename: "<<arg.structure_filename<<endl
	         <<"INFO:  - ROI filename: "<<arg.roi_filename<<endl
	         <<"INFO:  - Output filename: "<<arg.output_filename<<endl
	         <<"INFO:  - Number of threads in parallel: "<<arg.num_threads<<endl
	         <<"INFO:  - Spacing for mesh grid: "<<arg.spacing_grid<<" Angstrom"<<endl
	         <<"INFO:  - Spacing for grouping nearby atoms: "<<arg.spacing_group<<" Angstrom"<<endl
	         <<"INFO:  - Size of parallel blocks: "<<arg.size_pblock<<endl;

	    // Pass arguments to fvsearcher
	    fvsearcher.set_arguments(arg);

	    // Execute
    	double t_start, t_end;
		t_start = omp_get_wtime();
    	execute(fvsearcher, test_flag);
	    t_end = omp_get_wtime();
    	cout <<"INFO: Normal exit. ("<<fixed << std::setprecision(3)<<t_end - t_start<<" seconds total)"<<endl;

	}
}