#ifndef FVSEARCHER_H_
#define FVSEARCHER_H_

#include <omp.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>

#include <vector>
#include <map>
#include <set>
#include <utility> // std::pair

#include <math.h>

#include "common.hpp"

using namespace std;

// Structure of region of interest
struct ROI {double xb, xe;double yb, ye;double zb, ze;};

// Structure of an atom
struct Atom {double x, y, z, r;};

// Arguments for FVSearcher
class FVSArguments{
	public:
	    char * structure_filename; // Structural file input name
	    char * roi_filename; // ROI input file name
	    string output_filename;
	    int num_threads; // Number of threads in parallel
	    double spacing_grid; // Spacing for mesh grid
	    double spacing_group; // Spacing for grouping nearby atoms
	    double size_pblock; // Size of parallel blocks
	    double SCALING_FACTOR;
	    FVSArguments() { // Default argumetns
	    	output_filename = "__NO_OUTPUT__";
			num_threads = 4; 
			spacing_grid = 0.05; 
			spacing_group = 1.0; 
			size_pblock = 100;
			SCALING_FACTOR = 1.0;
	    }
};

// Table of equivalency for serial CCL
class EQ_serial{
	private:
		vector <int> P; // Equivalency tree
		int P_count; // Count of independent labels
	public:
		EQ_serial() {P.push_back(0); }// 1-base

		// Assign new provisional label
		int assign_new_label() {
			int count = P.size();
			P.push_back(count);
			return count;
		}

		// Find the root of a given node
		int find_root(int i) {
			while (P[i] < i) i = P[i];
			return i;
		}

		// Make all nodes in the path of node i to root
		void set_root(int i, int root) {
			while (P[i] < i) {
				int j = P[i];
				P[i] = root;
				i = j;
			}
			P[i] = root;
		}

		// Find the root of a given node i and compress the path in progress
		int find (int i) {
			int root = find_root(i);
			set_root(i, root);
			return root;
		}

		// Union two trees containing nodes i,j
		int merge (int i, int j) {
			int root = find_root(i);
			if (i != j) {
				int rootj = find_root(j);
				if (root > rootj)
					root = rootj;
				set_root(j, root);
			}
			set_root(i, root);
			return root;
		}

		// Flatten to produce final label
		void flatten() {
			int k = 1;
			P[0] = 0;
			for (int i = 1; i < P.size(); i++) {
				if (P[i] != i) 
					P[i] = P[P[i]];
				else
					P[i] = k++;
			}
			P_count = k - 1;
		}

		// Get an element
		int get(int i) {return P[i];}

		// Return size
		int count() {return P_count;}
};

// Table of equivalency for parallel CCL
class EQ_parallel{
	public:
		set < pair <int, int> > edges; // Record edges of global equivalencies
		int * P; // Record calculated global equivalencies

		int num_blocks;
		int * num_labels_provisional_local; // Record the number of local provisional labels in each parallel block
		int num_labels_provisional_global; // Record the number of global provisional labels
		int num_labels_global; // Record the number of final global labels
		map < pair <int, int>, int > local2global; // Record the connection between global and local provisional labels
		map < int, pair <int, int> > global2local; // Record the connection between local and global
	public:
		// Initialization
		EQ_parallel () {}
		void initialize (int num_blocks_) {
			num_blocks = num_blocks_;
			num_labels_provisional_local = new int [num_blocks_];
			num_labels_provisional_global = 0;
			local2global.clear();
			global2local.clear();
			edges.clear();
			P = 0;
		}

		// Get attributes
		int get_num_labels_provisional_global() {return num_labels_provisional_global;}
		int get_num_global_equivalencies() {return edges.size();}
		int get_num_labels_global() {return num_labels_global;}

		// Set num_labels_provisional_local
		void set_num_labels_provisional_local(int flat, int count) {
			num_labels_provisional_local[flat] = count;
		}

		// Buildup global2local and local2global
		void process_local_global_labels() {
	 		global2local.clear();
		    local2global.clear();
		    int global_label = 1;
		    for(int flat = 0; flat < num_blocks; flat++) {
		    	for (int local_label = 1; local_label <= num_labels_provisional_local[flat]; local_label++) {
		    		pair <int, int> local_id (flat, local_label);
		    		global2local[global_label] = local_id;
		    		local2global[local_id] = global_label;
		    		global_label++;
		    	}
		    }
		    num_labels_provisional_global = global_label - 1;
		}

		// Add edge
		void add_edge(int flat0, int local0, int flat1, int local1) {
			int global0 = local2global[pair <int, int> (flat0, local0)];
			int global1 = local2global[pair <int, int> (flat1, local1)];
			edges.insert(pair <int, int> (global0, global1));
		}

		// Process edges to achieve equivalencies
		void process(){
			// First pass
		    vector < pair <int, int> > U;
		    vector < pair <int, int> > E;

		    // Move from set to vector for parallelization
		    for (set < pair <int, int> > :: iterator edge = edges.begin(); 
			    	edge != edges.end(); edge++)
		    	E.push_back(*edge);

		    // Initialize P to be itself
		    P = new int [num_labels_provisional_global + 1]; // 1-base
		    #pragma omp for 
		    for (int i = 1; i <= num_labels_provisional_global; i++){
		    	P[i] = i;
		    }

		    // Process
		    do{
		    	// First pass, merge nodes in each edge
			    #pragma omp for
			    for (int i = 0; i < E.size(); i++){
			    	int x = E[i].first, y = E[i].second;
			    	while (P[x] != P[y]) {
			    		if (P[x] > P[y]){
			    			if (x == P[x]){
			    				#pragma omp atomic write
			    				P[x] = P[y];
			    				#pragma omp critical
			    				U.push_back(pair <int, int> (x, y));
			    				break;
			    			}
			    			int z = P[x];
			    			#pragma omp atomic write
			    			P[x] = P[y];
			    			x = z;
			    		}else{
			    			if (y == P[y]){
			    				#pragma omp atomic write
			    				P[y] = P[x];
			    				#pragma omp critical
			    				U.push_back(pair <int, int> (x, y));
			    				break;
			    			}
			    			int z = P[y];
			    			#pragma omp atomic write
			    			P[y] = P[x];
			    			y = z;
			    		}
			    	}
			    }

			    // Second pass, check for connectivites
			    #pragma omp barrier
			    while(U.size() != 0){
			    	int x,y; pair <int, int> edge;
			    	#pragma omp critical
			    	{
				    	edge = U.back();
				    	U.pop_back();
			    	}
			    	x = edge.first; y = edge.second;
			    	while (P[x] != P[y]){
			    		if (P[x] > P[y]){
			    			if (P[x] == x){
			    				#pragma omp critical
			    				U.push_back(pair <int, int> (x, y));
			    			}
			    			x = P[x];
			    		}else{
			    			if (P[y] == y){
			    				#pragma omp critical
			    				U.push_back(pair <int, int> (x, y));
			    			}
			    			y = P[y];
			    		}
			    	}

			    }
			}while(U.size() != 0); // Check whether all nodes were properly processed

			// Third pass, flatten
			int k = 1;
			P[0] = 0;
			for (int i = 1; i <= num_labels_provisional_global; i++) {
				if (P[i] != i) 
					P[i] = P[P[i]];
				else
					P[i] = k++;
			}
			num_labels_global = k - 1;
		}

		// Return root
		int root(int flat, int local_id){
			return P[local2global[pair <int, int> (flat, local_id)]];
		}
};

// Main class
class FVSearcher{
private:

// Variables

	vector <Atom> A; // Vector storaging all atoms
	ROI roi; // Region of interest
	map <string, double> van_der_Waals_table; // Van der Waals radius of common atoms

	double r_max;// Maximal van der Waals radius
	vector <Atom> * T; // Table of nearby atoms

	bool * M_full; // Array of occupancy (full array)
	set <int> M_sparse; // Array of occupancy (sparse array) 

	int * L_full; // Labeled array
	EQ_parallel EQ; // Process equivalencies
	

// Methods

	// Initialize table of van der Waals radius for common atoms
	void initialize_van_der_Waals() {
		van_der_Waals_table["H"] = 1.2;
		van_der_Waals_table["C"] = 1.7;
		van_der_Waals_table["N"] = 1.55;
		van_der_Waals_table["O"] = 1.52;
		van_der_Waals_table["F"] = 1.47;
		van_der_Waals_table["P"] = 1.8;
		van_der_Waals_table["S"] = 1.8;
		van_der_Waals_table["Cl"] = 1.75;
	}

	// Check the van der Waals radius of a given atom
	double van_der_Waals(string atomtype) {
		if (atomtype.size()  == 0) {
			throw runtime_error("Unrecognized atom type.");
		}else if (atomtype == "CLA") { // Chloride
			return arg.SCALING_FACTOR * van_der_Waals_table.find("Cl")->second;
		}else{ // Other atoms
			std::map<string, double>::iterator item =  \
				van_der_Waals_table.find(atomtype.substr(0,1));
			if (item == van_der_Waals_table.end())
				throw runtime_error("Unrecognized atom type.");
			else
				return arg.SCALING_FACTOR * item->second;
		}
	}

	// Check if the atom occupies the position
	inline bool is_occupied(Atom atom, double x, double y, double z) {
		return (pow(x - atom.x, 2) + 
			pow(y - atom.y, 2) + 
			pow(z - atom.z, 2)) <
			pow(atom.r, 2);
	}

	// Method for building full array of occupancy by interating through points
	void full_pointwise(int flat) {
		// Border of the table
	    int it = int(ceil((roi.xe-roi.xb + 2 * r_max) / arg.spacing_group));
	    int jt = int(ceil((roi.ye-roi.yb + 2 * r_max) / arg.spacing_group));
	    int kt = int(ceil((roi.ze-roi.zb + 2 * r_max) / arg.spacing_group));

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));

		// Flat id to xyz 
		int i,j,k;
		flat2ijk(it,jt,kt, i,j,k, flat);

		// Translate the current block in the table to the subarray in the array of occupancy
		int imbegin = fmax(int(ceil((i * arg.spacing_group - r_max) / arg.spacing_grid)), 0);
		int imend = fmin(int(ceil(((i + 1) * arg.spacing_group - r_max) / arg.spacing_grid)), im);
		int jmbegin = fmax(int(ceil((j * arg.spacing_group - r_max) / arg.spacing_grid)), 0);
		int jmend = fmin(int(ceil(((j + 1) * arg.spacing_group - r_max) / arg.spacing_grid)), jm);
		int kmbegin = fmax(int(ceil((k * arg.spacing_group - r_max) / arg.spacing_grid)), 0);
		int kmend = fmin(int(ceil(((k + 1) * arg.spacing_group - r_max) / arg.spacing_grid)), km);

		// Inspect whether the subarray is within the ROI or not. 
		if (imbegin < imend and jmbegin < jmend and kmbegin < kmend) {
			// Local vector storaging nearby atoms
			vector <Atom> nearby_atoms;

			// Distance of block within which atoms should be considered
			int padding = int(ceil(r_max / arg.spacing_group));
			int itbegin = fmax(i - padding, 0);
			int itend = fmin(i + padding + 1, it);
			int jtbegin = fmax(j - padding, 0);
			int jtend = fmin(j + padding + 1, jt);
			int ktbegin = fmax(k - padding, 0);
			int ktend = fmin(k + padding + 1, kt);

			// Iterate over the nearby blocks to collect atoms
			for (int ii = itbegin; ii < itend; ii++)
				for (int jj = jtbegin; jj < jtend; jj++)
					for (int kk = ktbegin; kk < ktend; kk++) {
						if ((pow(fmax(abs(ii - i) - 1, 0), 2) + 
							 pow(fmax(abs(jj - j) - 1, 0), 2) + 
							 pow(fmax(abs(kk - k) - 1, 0), 2)) * 
							pow(arg.spacing_group, 2)<=  \
							pow(r_max, 2)) {
							int flat_cur;
							ijk2flat(it,jt,kt, ii,jj,kk, flat_cur);
							nearby_atoms.insert(nearby_atoms.end(), 
								T[flat_cur].begin(), T[flat_cur].end());
						}
					}
			// Pointwise update array of occupancy
			for (int ii = imbegin; ii < imend; ii++)
				for (int jj = jmbegin; jj < jmend; jj++)
					for (int kk = kmbegin; kk < kmend; kk++) {
						// Real x,y,z location
						double x = ii * arg.spacing_grid + roi.xb;
						double y = jj * arg.spacing_grid + roi.yb;
						double z = kk * arg.spacing_grid + roi.zb;
						bool not_occupied = true;
						// For each nearby atoms, check if it overlaps the current point
						for (vector <Atom> :: iterator atom = nearby_atoms.begin();
							atom != nearby_atoms.end(); atom++)
							if (is_occupied(*atom, x, y, z)) {
								not_occupied = false;
								break;
							}
						// Update to array
						int flat_cur;
						ijk2flat(im,jm,km, ii,jj,kk, flat_cur);
						M_full[flat_cur] = not_occupied;
					}
		}
	}

	// Method for buidng full array of occupancy by interating through atoms
	void full_atomwise(int flat) {

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));

		// Atoms within the current table block
		vector <Atom> nearby_atoms = T[flat];

		// Atomwise update array of occupancy
		for (vector <Atom> :: iterator atom = nearby_atoms.begin();
			atom != nearby_atoms.end(); atom++) {
			int imbegin = fmax(int(floor((atom->x - atom->r - roi.xb) / arg.spacing_grid)), 0);
			int imend = fmin(int(ceil((atom->x + atom->r - roi.xb) / arg.spacing_grid)), im);
			int jmbegin = fmax(int(floor((atom->y - atom->r - roi.yb) / arg.spacing_grid)), 0);
			int jmend = fmin(int(ceil((atom->y + atom->r - roi.yb) / arg.spacing_grid)), jm);
			int kmbegin = fmax(int(floor((atom->z - atom->r - roi.zb) / arg.spacing_grid)), 0);
			int kmend = fmin(int(ceil((atom->z + atom->r - roi.zb) / arg.spacing_grid)), km);
			for (int ii = imbegin; ii < imend; ii++)
				for (int jj = jmbegin; jj < jmend; jj++)
					for (int kk = kmbegin; kk < kmend; kk++) {
						int flat_cur;
						ijk2flat(im,jm,km, ii,jj,kk, flat_cur);
						if (M_full[flat_cur]) {
							double x = ii * arg.spacing_grid + roi.xb;
							double y = jj * arg.spacing_grid + roi.yb;
							double z = kk * arg.spacing_grid + roi.zb;
							if (is_occupied(*atom, x, y, z)) {
								#pragma omp atomic write
								M_full[flat_cur] = false;
							}
						}
					}
		}
	}

	// Method for building full array of occupancy by interating through atoms
	void full_atomwise_old(int n) {

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));

		// Atom to proceed
		Atom atom = A[n];
					
		// Border of the atom
		int imbegin = fmax(int(floor((atom.x - atom.r - roi.xb) / arg.spacing_grid)), 0);
		int imend = fmin(int(ceil((atom.x + atom.r - roi.xb) / arg.spacing_grid)), im);
		int jmbegin = fmax(int(floor((atom.y - atom.r - roi.yb) / arg.spacing_grid)), 0);
		int jmend = fmin(int(ceil((atom.y + atom.r - roi.yb) / arg.spacing_grid)), jm);
		int kmbegin = fmax(int(floor((atom.z - atom.r - roi.zb) / arg.spacing_grid)), 0);
		int kmend = fmin(int(ceil((atom.z + atom.r - roi.zb) / arg.spacing_grid)), km);

		// Proceed
		for (int ii = imbegin; ii < imend; ii++)
			for (int jj = jmbegin; jj < jmend; jj++)
				for (int kk = kmbegin; kk < kmend; kk++) {
					int flat_cur;
					ijk2flat(im,jm,km, ii,jj,kk, flat_cur);
					if (M_full[flat_cur]) {
						double x = ii * arg.spacing_grid + roi.xb;
						double y = jj * arg.spacing_grid + roi.yb;
						double z = kk * arg.spacing_grid + roi.zb;
						if (is_occupied(atom, x, y, z)) {
							#pragma omp atomic write
							M_full[flat_cur] = false;
						}
					}
				}
	}

	// Method for building sparse array of occupancy by interating through points
	void sparse(int flat) {
		// Border of the table
	    int it = int(ceil((roi.xe-roi.xb + 2 * r_max) / arg.spacing_group));
	    int jt = int(ceil((roi.ye-roi.yb + 2 * r_max) / arg.spacing_group));
	    int kt = int(ceil((roi.ze-roi.zb + 2 * r_max) / arg.spacing_group));

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));

		// Flat id to xyz 
		int i,j,k;
		flat2ijk(it,jt,kt, i,j,k, flat);

		// Translate the current block in the table to the subarray in the array of occupancy
		int imbegin = fmax(int(ceil((i * arg.spacing_group - r_max) / arg.spacing_grid)), 0);
		int imend = fmin(int(ceil(((i + 1) * arg.spacing_group - r_max) / arg.spacing_grid)), im);
		int jmbegin = fmax(int(ceil((j * arg.spacing_group - r_max) / arg.spacing_grid)), 0);
		int jmend = fmin(int(ceil(((j + 1) * arg.spacing_group - r_max) / arg.spacing_grid)), jm);
		int kmbegin = fmax(int(ceil((k * arg.spacing_group - r_max) / arg.spacing_grid)), 0);
		int kmend = fmin(int(ceil(((k + 1) * arg.spacing_group - r_max) / arg.spacing_grid)), km);

		// Inspect whether the subarray is within the ROI or not. 
		if (imbegin < imend and jmbegin < jmend and kmbegin < kmend) {
			// Local vector storaging nearby atoms
			vector <Atom> nearby_atoms;

			// Distance of block within which atoms should be considered
			int padding = int(ceil(r_max / arg.spacing_group));
			int itbegin = fmax(i - padding, 0);
			int itend = fmin(i + padding + 1, it);
			int jtbegin = fmax(j - padding, 0);
			int jtend = fmin(j + padding + 1, jt);
			int ktbegin = fmax(k - padding, 0);
			int ktend = fmin(k + padding + 1, kt);

			// Iterate over the nearby blocks to collect atoms
			for (int ii = itbegin; ii < itend; ii++)
				for (int jj = jtbegin; jj < jtend; jj++)
					for (int kk = ktbegin; kk < ktend; kk++) {
						if ((pow(fmax(abs(ii - i) - 1, 0), 2) + 
							 pow(fmax(abs(jj - j) - 1, 0), 2) + 
							 pow(fmax(abs(kk - k) - 1, 0), 2)) * 
							pow(arg.spacing_group, 2)<=  \
							pow(r_max, 2)) {
							int flat_cur;
							ijk2flat(it,jt,kt, ii,jj,kk, flat_cur);
							nearby_atoms.insert(nearby_atoms.end(), 
								T[flat_cur].begin(), T[flat_cur].end());
						}
					}
			// Pointwise update array of occupancy
			for (int ii = imbegin; ii < imend; ii++)
				for (int jj = jmbegin; jj < jmend; jj++)
					for (int kk = kmbegin; kk < kmend; kk++) {
						// Real x,y,z location
						double x = ii * arg.spacing_grid + roi.xb;
						double y = jj * arg.spacing_grid + roi.yb;
						double z = kk * arg.spacing_grid + roi.zb;
						bool not_occupied = true;
						// For each nearby atoms, check if it overlaps the current point
						for (vector <Atom> :: iterator atom = nearby_atoms.begin();
							atom != nearby_atoms.end(); atom++)
							if (is_occupied(*atom, x, y, z)) {
								not_occupied = false;
								break;
							}
						// Update to array
						if (not_occupied){
							int flat_cur;
							ijk2flat(im,jm,km, ii,jj,kk, flat_cur);
							#pragma omp critical
							M_sparse.insert(flat_cur);
						}
					}
		}
	}

	void CCL_serial_two_pass(int flat) {

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));

	    // Border of blocks
	    int il = int(ceil(1.0 * im / arg.size_pblock));
	    int jl = int(ceil(1.0 * jm / arg.size_pblock));
	    int kl = int(ceil(1.0 * km / arg.size_pblock));

	    // Current block
	    int i,j,k;
	    flat2ijk(il,jl,kl, i,j,k, flat);

	    int ibegin = i * arg.size_pblock;
	    int iend = fmin((i + 1) * arg.size_pblock, im);
	    int jbegin = j * arg.size_pblock;
	    int jend = fmin((j + 1) * arg.size_pblock, jm);
	    int kbegin = k * arg.size_pblock;
	    int kend = fmin((k + 1) * arg.size_pblock, km);

	    // Equivalence
	    EQ_serial eq;

	    // First pass, scanning and giving provisional labels
	    for (int ii = ibegin; ii < iend; ii++)
	    	for (int jj = jbegin; jj < jend; jj++)
	    		for (int kk = kbegin; kk < kend; kk++) {
	    			// Retrieve flat id
	    			int flat_cur;
	    			ijk2flat(im,jm,km, ii,jj,kk, flat_cur);
	    			// If not occupied
	    			if (M_full[flat_cur]) {

	    				// Record nearby labels
		    			vector <int> nl;

		    			// Retrieve nearby labels
		    			if (ii > ibegin and M_full[flat_cur - jm * km]) nl.push_back(L_full[flat_cur - jm * km]);
	    				if (jj > jbegin and M_full[flat_cur - km]) nl.push_back(L_full[flat_cur - km]);
						if (kk > kbegin and M_full[flat_cur - 1]) nl.push_back(L_full[flat_cur - 1]);

						// Apply labels
						if (nl.size() == 0) { 
							// If no neighbor found, apply a new label
							L_full[flat_cur] = eq.assign_new_label();

						} else {

							// Find the minimal label and apply
							int min_label = arg.size_pblock * arg.size_pblock * arg.size_pblock;
							for (int i = 0; i < nl.size(); i++) {
								int root = eq.find_root(nl[i]);
								if (root < min_label) 
									min_label = root;
							}

							L_full[flat_cur] = min_label;

							// Record the equivalence of labels
							for (int i = 0; i < nl.size(); i++)
								eq.merge(nl[i], min_label);
						}
		    		}
	    		}

	    // Finalize equivalency
	    eq.flatten();
	    EQ.set_num_labels_provisional_local(flat, eq.count());

	    // Second pass, applying equivalencies
	    for (int ii = ibegin; ii < iend; ii++)
	    	for (int jj = jbegin; jj < jend; jj++)
	    		for (int kk = kbegin; kk < kend; kk++) {
	    			int flat_cur;
	    			ijk2flat(im,jm,km, ii,jj,kk, flat_cur);
	    			L_full[flat_cur] = eq.get(L_full[flat_cur]);
	    		}
	}

	void CCL_retrieve_equivalencies(int flat) {

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
	    int imjmkm = im * jm * km;

	    // Border of blocks
	    int il = int(ceil(1.0 * im / arg.size_pblock));
	    int jl = int(ceil(1.0 * jm / arg.size_pblock));
	    int kl = int(ceil(1.0 * km / arg.size_pblock));

	    // Location of current block
	    int i,j,k;
	    flat2ijk(il,jl,kl, i,j,k, flat);

	    // Check for equivalencies in between neighbor blocks
    	// x surface
    	if (i > 0) {
    		// Neighbor block
    		int flat_neighbor;
    		ijk2flat(il,jl,kl, i-1,j,k, flat_neighbor);
			
    		// Surface
    		int ib = i * arg.size_pblock - 1;
    		int ie = i * arg.size_pblock;
		    int jbegin = j * arg.size_pblock;
		    int jend = fmin((j + 1) * arg.size_pblock, jm);
		    int kbegin = k * arg.size_pblock;
		    int kend = fmin((k + 1) * arg.size_pblock, km);

		    // Iterate through surface
		    for (int jj = jbegin; jj < jend; jj++)
		    	for (int kk = kbegin; kk < kend; kk++) {
		    		// Retrieve two points
		    		int flatb, flate;
		    		ijk2flat(im,jm,km, ib,jj,kk, flatb);
		    		ijk2flat(im,jm,km, ie,jj,kk, flate);
		    		// If both unoccupied
		    		if (M_full[flatb] and M_full[flate]) {
		    			// Retrieve their labels
			    		int labelb = L_full[flatb];
		    			int labele = L_full[flate];
		    			// Add equivalency
		    			#pragma omp critical
		    			EQ.add_edge(flat_neighbor, labelb, flat, labele);
		    		}
		    	}
    	} 
    	// y surface
    	if (j > 0) {
    		// Neighbor block
    		int flat_neighbor;
    		ijk2flat(il,jl,kl, i,j-1,k, flat_neighbor);
			
    		// Surface
		    int ibegin = i * arg.size_pblock;
		    int iend = fmin((i + 1) * arg.size_pblock, im);
    		int jb = j * arg.size_pblock - 1;
    		int je = j * arg.size_pblock;
		    int kbegin = k * arg.size_pblock;
		    int kend = fmin((k + 1) * arg.size_pblock, km);

		    // Iterate through surface
		    for (int ii = ibegin; ii < iend; ii++)
		    	for (int kk = kbegin; kk < kend; kk++) {
		    		// Retrieve two points
		    		int flatb, flate;
		    		ijk2flat(im,jm,km, ii,jb,kk, flatb);
		    		ijk2flat(im,jm,km, ii,je,kk, flate);
		    		// If both unoccupied
		    		if (M_full[flatb] and M_full[flate]) {
		    			// Retrieve their labels
			    		int labelb = L_full[flatb];
		    			int labele = L_full[flate];
		    			// Add equivalency
		    			#pragma omp critical
		    			EQ.add_edge(flat_neighbor, labelb, flat, labele);
		    		}
		    	}
    	} 
    	if (k > 0) {
    		// Neighbor block
    		int flat_neighbor;
    		ijk2flat(il,jl,kl, i,j,k-1, flat_neighbor);
			
    		// Surface
		    int ibegin = i * arg.size_pblock;
		    int iend = fmin((i + 1) * arg.size_pblock, im);
		    int jbegin = j * arg.size_pblock;
		    int jend = fmin((j + 1) * arg.size_pblock, jm);
    		int kb = k * arg.size_pblock - 1;
    		int ke = k * arg.size_pblock;

		    // Iterate through surface
		    for (int ii = ibegin; ii < iend; ii++)
			    for (int jj = jbegin; jj < jend; jj++) {
		    		// Retrieve two points
		    		int flatb, flate;
		    		ijk2flat(im,jm,km, ii,jj,kb, flatb);
		    		ijk2flat(im,jm,km, ii,jj,ke, flate);
		    		// If both unoccupied
		    		if (M_full[flatb] and M_full[flate]) {
		    			// Retrieve their labels
			    		int labelb = L_full[flatb];
		    			int labele = L_full[flate];
		    			// Add equivalency
		    			#pragma omp critical
		    			EQ.add_edge(flat_neighbor, labelb, flat, labele);
		    		}
		    	}
    	} 
	}

public:
	FVSArguments arg; // Input arguments
	// Init
	FVSearcher () {initialize_van_der_Waals();}
	FVSearcher (FVSArguments arg_) {arg = arg_;	initialize_van_der_Waals();}

	// Retrieve elements in class
	ROI get_roi() {return roi;};
	int get_num_atoms() {return A.size();}
	double get_r_max() {return r_max;}
	int get_T_size() { return \
		int(ceil((roi.xe-roi.xb + 2 * r_max) / arg.spacing_group)) * 
	    int(ceil((roi.ye-roi.yb + 2 * r_max) / arg.spacing_group)) * 
	    int(ceil((roi.ze-roi.zb + 2 * r_max) / arg.spacing_group));
	}
	int get_M_full_size() { return \
	    int(ceil((roi.xe-roi.xb) / arg.spacing_grid)) *
	    int(ceil((roi.ye-roi.yb) / arg.spacing_grid)) *
	    int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
	}
	int get_num_atoms_in_ROI() { // Count how many atoms are in the ROI
	    // Count
	    int count = 0;
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared) reduction(+:count)
	    for(int n = 0; n < A.size(); n++){
	    	Atom atom = A[n];
	    	if (atom.x + atom.r > roi.xb and 
				atom.x - atom.r < roi.xe and 
				atom.y + atom.r > roi.yb and 
				atom.y - atom.r < roi.ye and 
				atom.z + atom.r > roi.zb and 
				atom.z - atom.r < roi.ze)
	    		count++;
	    }
	    return count;
	}
	int get_num_labels_provisional() {return EQ.get_num_labels_provisional_global();}
	int get_num_global_equivalencies() {return EQ.get_num_global_equivalencies();}
	int get_num_labels_global() {return EQ.get_num_labels_global();}
	string get_output_filename() {return arg.output_filename;}

	// Set arguments
	void set_arguments(FVSArguments arg_) {arg = arg_;}

	// Read input files
	void read_inputs() {

    	ifstream rfile;
		string line;

		// Read ROI file
		rfile.open(arg.roi_filename);
		if (rfile.is_open()) {
			for(int i = 0; i < 3; i++) {
				string xyz;
				double begin,end;
				rfile >>xyz>>begin>>end;
				if (xyz == "X") {roi.xb = begin; roi.xe = end;}
				else if (xyz == "Y") {roi.yb = begin; roi.ye = end;}
				else if (xyz == "Z") {roi.zb = begin; roi.ze = end;}
				else{
					rfile.close();
					stringstream ss;
					ss <<"ROI file '"<< arg.roi_filename << "'' was invalid!"<<endl
					   <<"An example of an acceptable ROI file is here"<<endl
					   <<"(BEGIN of the File)"<<endl
					   <<"X -30.0 30.0"<<endl
					   <<"Y -30.0 30.0"<<endl
					   <<"Z -20.0 20.0"<<endl
					   <<"(END of the file)"<<endl;
					throw runtime_error(ss.str());
				}
			}
		}else{
			stringstream ss;
			ss <<"ROI file '"<< arg.roi_filename << "'' does not exist!";
			rfile.close();
			throw runtime_error(ss.str());
		}
		rfile.close();

		// Read input structure

		// Format of protein data bank
		// The line starting with "ATOM" indicating that this line included an atom
		// The atom type was in line[12:16]
		// x,y,z were in line[31,38], line[39,46],[47,54]
		if (ends_with(string(arg.structure_filename).substr(), ".pdb") or \
			ends_with(string(arg.structure_filename).substr(), ".PDB")) {
			rfile.open(arg.structure_filename);
			if (rfile.is_open()) {
				initialize_van_der_Waals();
				while( getline (rfile, line)) {
					if (line.compare(0,4,"ATOM") == 0) {

						string atomtype = trim_string(line.substr(12,4).c_str());
						double x = atof(line.substr(31,7).c_str());
						double y = atof(line.substr(39,7).c_str());
						double z = atof(line.substr(47,7).c_str());

						Atom cur;
						cur.x = x;
						cur.y = y;
						cur.z = z;
						cur.r = van_der_Waals(atomtype);
						A.push_back(cur);
					}
				}
			}else{
				rfile.close();
				stringstream ss;
				ss <<"Structure file '"<< arg.structure_filename << "'' does not exist!";
				throw runtime_error(ss.str());
			}
			rfile.close();

		// Format of a simple colume table of 
		// X Y Z R
		// Such as 
		// 1. 2. 3. 1.2
		// For an atom at x,y,z = (1,2,3), radius = 1.2
		}else if (ends_with(string(arg.structure_filename).substr(), ".xyz") or \
			ends_with(string(arg.structure_filename).substr(), ".XYZ")) {
			rfile.open(arg.structure_filename);
			if (rfile.is_open()) {
				initialize_van_der_Waals();
				while( getline (rfile, line)) {
					double x, y, z, r;
					stringstream ss;
					r = -1;
					ss <<line;
					ss >>x>>y>>z>>r;
					if (r >= 0) {
						Atom cur;
						cur.x = x;
						cur.y = y;
						cur.z = z;
						cur.r = r;
						A.push_back(cur);
					}
				}
			}else{
				rfile.close();
				stringstream ss;
				ss <<"Structure file '"<< arg.structure_filename << "'' does not exist!";
				throw runtime_error(ss.str());
			}
			rfile.close();
		}

	    // Check ROI
	    if (roi.xe <= roi.xb or 
	    	roi.ye <= roi.yb or 
	    	roi.ze <= roi.zb) {
			stringstream ss;
			ss <<"Invalid ROI."<<endl;
			throw runtime_error(ss.str());
	    }

	    // Check atoms
	    if (A.size() == 0) {
	    	stringstream ss;
			ss <<"No atoms in the structural file."<<endl;
			throw runtime_error(ss.str());
	    }
	}

	// Calculate the maximal radius in given atoms
	void calc_maximal_radius() {
		r_max = 0;
		int num_atoms = A.size();
		#pragma omp parallel for num_threads(arg.num_threads) default(shared)
		for (int n = 0; n < num_atoms; n++) {
			// Using a pre-check to avoid frequent entrances to critical region
			if (A[n].r > r_max) {
				#pragma omp critical
				if (A[n].r > r_max)
					r_max = A[n].r;
			}
		}
	}

	// Build the nearby-atom table
	void build_nearby_atom_table() {

		// Border of the ROI
	    int it = int(ceil((roi.xe-roi.xb + 2 * r_max) / arg.spacing_group));
	    int jt = int(ceil((roi.ye-roi.yb + 2 * r_max) / arg.spacing_group));
	    int kt = int(ceil((roi.ze-roi.zb + 2 * r_max) / arg.spacing_group));
	    int itjtkt = it * jt * kt;
	    
	    //  Initialize table of nearby atom table
	    T = new vector <Atom> [it * jt * kt];

	    //  Set synchronous write_locks
	    omp_lock_t * write_lock = new omp_lock_t [it * jt * kt];

		// Initialize all locks
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared)
	    for (int flat = 0; flat < itjtkt; flat++)
			omp_init_lock(write_lock + flat);

	    // Build up the table
	    int num_atoms = A.size();
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared)
		for (int n = 0; n < num_atoms; n++) {

			// Location of current atom in the table
			int i = int(floor((A[n].x - roi.xb + r_max) / arg.spacing_group));
			int j = int(floor((A[n].y - roi.yb + r_max) / arg.spacing_group));
			int k = int(floor((A[n].z - roi.zb + r_max) / arg.spacing_group));

			// Check if the atom is within ROI
			if (i >= 0 and i < it and 
				j >= 0 and j < jt and
				k >= 0 and k < kt) {

				// xyz to flat id
				int flat;
				ijk2flat(it,jt,kt, i,j,k, flat);

				// Lock to avoid race
				omp_set_lock(write_lock + flat);

				// Append the atom to the table
				T[flat].push_back(A[n]);

				// Unlock
				omp_unset_lock(write_lock + flat);
			}
		}
	}

	// Build array of occupancies - full, pointwise
	void build_array_occupancy_full_pointwise(){

		// Calculate maximal van der Waals radius
		calc_maximal_radius();

		// Build nearby-atom table
		build_nearby_atom_table();

		// Border of the table
	    int it = int(ceil((roi.xe-roi.xb + 2 * r_max) / arg.spacing_group));
	    int jt = int(ceil((roi.ye-roi.yb + 2 * r_max) / arg.spacing_group));
	    int kt = int(ceil((roi.ze-roi.zb + 2 * r_max) / arg.spacing_group));
	    int itjtkt = it * jt * kt;

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
	    int imjmkm = im * jm * km;

	    // New arary of occupancy
		M_full = new bool[im * jm * km];

	    // Build array of occupancies
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared) 
	    for (int flat = 0; flat < itjtkt; flat++)
	    	full_pointwise(flat);
	}

	// Build array of occupancies - full, atomwise
	void build_array_occupancy_full_atomwise(){

		// Calculate maximal van der Waals radius
		calc_maximal_radius();

		// Build nearby-atom table
		build_nearby_atom_table();

		// Border of the table
	    int it = int(ceil((roi.xe-roi.xb + 2 * r_max) / arg.spacing_group));
	    int jt = int(ceil((roi.ye-roi.yb + 2 * r_max) / arg.spacing_group));
	    int kt = int(ceil((roi.ze-roi.zb + 2 * r_max) / arg.spacing_group));
	    int itjtkt = it * jt * kt;

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
	    int imjmkm = im * jm * km;

	    // New arary of occupancy
		M_full = new bool[im * jm * km];
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared)
	    for (int flat = 0 ; flat < imjmkm; flat ++)
	    	M_full[flat] = true;

	    // Build array of occupancies
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared) 
	    for (int flat = 0; flat < itjtkt; flat++)
	    	full_atomwise(flat);
	}

	// Build array of occupancies - full, atomwise
	void build_array_occupancy_full_atomwise_old(){

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
	    int imjmkm = im * jm * km;

		// Initialize 
		M_full = new bool[im * jm * km];
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared)
	    for (int flat = 0 ; flat < imjmkm; flat ++)
	    	M_full[flat] = true;

	    // Calculate
		int num_atoms = A.size();
		#pragma omp parallel for num_threads(arg.num_threads) default(shared)
		for (int n = 0; n < num_atoms; n++){
			Atom atom = A[n];
			if (atom.x + atom.r > roi.xb and 
				atom.x - atom.r < roi.xe and 
				atom.y + atom.r > roi.yb and 
				atom.y - atom.r < roi.ye and 
				atom.z + atom.r > roi.zb and 
				atom.z - atom.r < roi.ze)
				full_atomwise(n);
		}
	}

	// Build array of occupancies - sparse, pointwise
	void build_array_occupancy_sparse(){

		// Calculate maximal van der Waals radius
		calc_maximal_radius();

		// Build nearby-atom table
		build_nearby_atom_table();

		// Border of the table
	    int it = int(ceil((roi.xe-roi.xb + 2 * r_max) / arg.spacing_group));
	    int jt = int(ceil((roi.ye-roi.yb + 2 * r_max) / arg.spacing_group));
	    int kt = int(ceil((roi.ze-roi.zb + 2 * r_max) / arg.spacing_group));
	    int itjtkt = it * jt * kt;

	    // New arary of occupancy
		M_sparse.clear();

	    // Build array of occupancies
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared) 
	    for (int flat = 0; flat < itjtkt; flat++)
	    	sparse(flat);
	}

	// Three methods for calculating array of occupancy 
	void build_array_occupancy(int full_sparse, int point_atom) {
		// full_sparse: 0 - full, 1- sparse
		// point_atom: 0 - pointwise, 1 - atomwise (irrelevant if sparse was selected)

	    if (full_sparse == 0) { // 1 and 2, Full array
			if (point_atom == 0) {// 1. Full array, poinwise
				build_array_occupancy_full_pointwise();
			}else if (point_atom == 1) { // 2. Full array, atomwise
				build_array_occupancy_full_atomwise();
			}else {
		    	throw runtime_error("build_array_occupancy: Method not recognized. ");
		    }
	    }else if (full_sparse == 1) { // 3. Sparse array
			build_array_occupancy_sparse();
	    }else {
	    	throw runtime_error("build_array_occupancy: Method not recognized. ");
	    }
	}

	// Validate the array of occupancy
	int validate_array_occupancy(int full_sparse) {
		if (full_sparse == 0) {
			// Border of the array
		    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
		    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
		    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
		    int imjmkm = im * jm * km;

		    // Count the number of atoms
		    int count = 0;
		    #pragma omp parallel for num_threads(arg.num_threads) default(shared) reduction(+:count)
		    for(int flat = 0; flat < imjmkm; flat++)
		    	if (M_full[flat]) count ++;
		    return count;
		}else if(full_sparse == 1) {
			return M_sparse.size();
		}else {
	    	throw runtime_error("build_array_occupancy: Method not recognized. ");
	    }
	}

	// Parallel connected-components labeling
	// Step 1 serial two-pass CCL for each block
	void parallel_CCL_serial_block() {

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
	    int imjmkm = im * jm * km;

	    // Border of blocks
	    int il = int(ceil(1.0 * im / arg.size_pblock));
	    int jl = int(ceil(1.0 * jm / arg.size_pblock));
	    int kl = int(ceil(1.0 * km / arg.size_pblock));
	    int iljlkl = il * jl * kl;

    	// New array of label, initializing all to 0
    	L_full = new int [im * jm * km];
    	#pragma omp parallel for num_threads(arg.num_threads) default(shared)
	    for (int flat = 0; flat < imjmkm; flat ++)
	    	L_full[flat] = 0;

		// Record equivalencies among blocks
	    EQ.initialize(iljlkl);

	    // Do CCL on each sub blocks
		#pragma omp parallel for num_threads(arg.num_threads) default(shared)
	    for (int flat = 0; flat < iljlkl; flat++) {
	    	CCL_serial_two_pass(flat);
	    }

	    // Build up map between local and global provisional labels
	    EQ.process_local_global_labels();
	}

	// Step 2 Retrieve equivalencies
	void parallel_CCL_retrieve_equivalencies() {

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
	    int imjmkm = im * jm * km;

	    // Border of blocks
	    int il = int(ceil(1.0 * im / arg.size_pblock));
	    int jl = int(ceil(1.0 * jm / arg.size_pblock));
	    int kl = int(ceil(1.0 * km / arg.size_pblock));
	    int iljlkl = il * jl * kl;

	    // Extract equivalencies
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared)
	    for (int flat = 0; flat < iljlkl; flat++) {
	    	CCL_retrieve_equivalencies(flat);
	    }

	    // Sort out equivalencies
	    #pragma omp paralle num_threads(arg.num_threads) default(shared)
	    EQ.process();
	}

	// Step 3 apply equivalencies
	void parallel_CCL_apply_equivalencies() {

		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
	    int imjmkm = im * jm * km;

	    // Border of blocks
	    int il = int(ceil(1.0 * im / arg.size_pblock));
	    int jl = int(ceil(1.0 * jm / arg.size_pblock));
	    int kl = int(ceil(1.0 * km / arg.size_pblock));

	    // Extract equivalencies
	    #pragma omp parallel for num_threads(arg.num_threads) default(shared)
	    for (int flat = 0; flat < imjmkm; flat++) {
	    	// Element id
	    	int i,j,k;
	    	flat2ijk(im,jm,km, i,j,k, flat);

	    	// Block id
	    	int ii = i / arg.size_pblock;
	    	int jj = j / arg.size_pblock;
	    	int kk = k / arg.size_pblock;
	    	int flat_block_id;
			ijk2flat(il,jl,kl, ii,jj,kk, flat_block_id);

	    	L_full[flat] = EQ.root(flat_block_id, L_full[flat]);
	    }
	}

	// Step 4 output
	void output(){
		// Border of the array
	    int im = int(ceil((roi.xe-roi.xb) / arg.spacing_grid));
	    int jm = int(ceil((roi.ye-roi.yb) / arg.spacing_grid));
	    int km = int(ceil((roi.ze-roi.zb) / arg.spacing_grid));
	    int imjmkm = im * jm * km;

	    // Output
		ofstream wfile;

		// ROI info
		wfile.open((arg.output_filename+".roi").c_str());
		// First line record ROI
		wfile <<roi.xb<<","<<roi.xe<<","<<roi.yb<<","<<roi.ye<<","<<roi.zb<<","<<roi.ze<<","<<endl;
		// Second line record spacing for mesh grid, x, y, z dimension of the matrix
		wfile <<arg.spacing_grid<<","<<im<<","<<jm<<","<<km<<","<<endl;
		wfile.close();

		// Array of labels
		wfile.open((arg.output_filename+".ccl").c_str(), ios::out | ios::binary );
		wfile.write((char*) L_full, sizeof(int) * imjmkm);
		wfile <<endl;
	}

	void exec(){
	    read_inputs();
		calc_maximal_radius();
	    build_nearby_atom_table();
		build_array_occupancy(0, 1);
		parallel_CCL_serial_block();
		parallel_CCL_retrieve_equivalencies();
		parallel_CCL_apply_equivalencies();
		output();
	}
};

#endif