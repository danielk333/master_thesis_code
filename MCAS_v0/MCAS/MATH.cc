/*
 * math_functions.cc
 *
 *  Created on: Sep 15, 2014
 *      Author: dankas
 */

#include <vector>
#include <math.h>
#include <algorithm>
#include <sstream>

 
#include "define.hh"
#include "functions.hh"
#include "class.hh"

double mean_mat(std::vector<std::vector<double> > M) {
	double ret;
	unsigned int i,j;
	unsigned long long int counter;
	ret = 0; counter = 0;
	for(i=0; i < M.size(); i++) {
		for(j=0; j < M[i].size(); j++) {
			counter++;
			ret = ret + M[i][j];
		}
	}
	ret = ret/counter;

	return ret;
}

double stdv_mat(std::vector<std::vector<double> > M, double M_mean) {
	double ret;
	unsigned int i,j;
	unsigned long long int counter;
	ret = 0; counter = 0;
	for(i=0; i < M.size(); i++) {
		for(j=0; j < M[i].size(); j++) {
			counter++;
			ret = ret + M[i][j]*M[i][j];
		}
	}
	ret = ret/counter - M_mean;

	return ret;
}

std::vector<std::vector<double> > fill_mat(unsigned int a, unsigned int b, double X) {
	std::vector<std::vector<double> > ret;
	std::vector<double> temp;
	unsigned int i,j;

	for(i=0; i < a; i++) {
		temp.clear();
		for(j=0; j < b; j++) {
			temp.push_back(X);
		}
		ret.push_back(temp);
	}

	return ret;
}

double abs_double(double a) {
	return abs(a);
}

double JDtoJ(double JD) {
	return 2000.0 + (JD - 2451545.0)/365.25;
}
double DATEtoJD(double year, double month, double day, double hour, double minute, double second) {
	double a = floor((14.0 - month)/12.0);
	double y = year + 4800.0 - a;
	double m = month + 12.0*a - 3.0;
	double JDN = day + floor((153.0*m + 2.0)/5.0) + 365.0*y + floor(y/4.0) - floor(y/100.0) + floor(y/400.0) - 32045.0;

	double JD = JDN + (hour - 12.0)/24.0 + minute/1440.0 + second/86400.0;
	return JD;
}

std::vector<double> rot_to_plane(std::vector<double> u ,std::vector<double> n) {
	std::vector<double> s;

	double r = abs_v(n);
	double phi = acos(n[2]/r);
	double theta = atan2(n[1],n[0]);

	s = rot_y(rot_z(u,-theta),-phi);

	return s;
}

std::vector<double> rot_z(std::vector<double> u ,double theta) {
	std::vector<double> s;
	
	s.push_back(u[0]*cos(theta) - u[1]*sin(theta));
	s.push_back(u[0]*sin(theta) + u[1]*cos(theta));
	s.push_back(u[2]);

	return s;
}

std::vector<double> rot_y(std::vector<double> u ,double theta) {
	std::vector<double> s;
	
	s.push_back(u[0]*cos(theta) + u[2]*sin(theta));
	s.push_back(u[1]);
	s.push_back(-u[0]*sin(theta) + u[2]*cos(theta));

	return s;
}

std::vector<double> cross_product(std::vector<double> u,std::vector<double> v) {
	std::vector<double> s;
	s.push_back(u[1]*v[2]-u[2]*v[1]);
	s.push_back(u[2]*v[0]-u[0]*v[2]);
	s.push_back(u[0]*v[1]-u[1]*v[0]);

	return s;
}

double dot_product(std::vector<double> u,std::vector<double> v) {
	return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}

double abs_v(std::vector<double> u) {
	return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}
/*
std::vector<double> multiply_v(double a, std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(u[i]*a);
	}
	return ret_vec;
}

std::vector<double> add_v(double a, std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(u[i]+a);
	}
	return ret_vec;
}

std::vector<std::vector<double> > multiply_mat(double a, std::vector<std::vector<double> > M) {
	std::vector<double> temp_vec;
	std::vector<std::vector<double> > return_vec;
	unsigned int i,j;
	for(i=0; i < M.size(); i++) {
		temp_vec.clear();
		for(j=0; j < M[i].size(); j++) {
			temp_vec.push_back(a*((M[i])[j]));
		}
		return_vec.push_back(temp_vec);
	}

	return return_vec;
}

std::vector<std::vector<double> > add_v_mat(std::vector<double> a, std::vector<std::vector<double> > M) {
	std::vector<double> temp_vec;
	std::vector<std::vector<double> > return_vec;
	unsigned int i,j;
	for(i=0; i < M.size(); i++) {
		temp_vec.clear();
		for(j=0; j < M[i].size(); j++) {
			temp_vec.push_back(a[j] + ((M[i])[j]));
		}
		return_vec.push_back(temp_vec);
	}

	return return_vec;
}

std::vector<double> add_v_v(std::vector<double> a, std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(u[i]+a[i]);
	}
	return ret_vec;
}

std::vector<double> sub_v_v(std::vector<double> a, std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(u[i]-a[i]);
	}
	return ret_vec;
}
*/
std::vector<double> abs_element_v(std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(abs(u[i]));
	}
	return ret_vec;
}

std::vector<double> extract_col(std::vector<std::vector<double> > M, unsigned int id) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<M.size();i++){
		ret_vec.push_back(M[i][id]);
	}
	return ret_vec;
}


std::string IntToString(unsigned int a) {
    std::ostringstream temp;
    temp<<a;
    return temp.str();
}


std::vector<unsigned int> histogram(const std::vector<double> *data, unsigned int n) {
	if(n == 0) {
		n = round(sqrt((*data).size()));
	}
	unsigned int i,j;
	std::vector<double> temp_data_container = *data;
	std::vector<double>::iterator MAX_ID, MIN_ID;
	MIN_ID = std::min_element(temp_data_container.begin(),temp_data_container.end());
	MAX_ID = std::max_element(temp_data_container.begin(),temp_data_container.end());
	double Delta = (*MAX_ID - *MIN_ID)/n;

	std::vector<unsigned int> index_vec;
	std::vector<unsigned int> ret;
	for(i=0; i < temp_data_container.size(); i++) {
		index_vec.push_back(i);
	}

	for(i=0; i < n; i++) {
		ret.push_back(0);
		for(j=index_vec.size(); j-- > 0; ) {
			if(temp_data_container[index_vec[j]] >= (*MIN_ID + ((double)i)*Delta) && temp_data_container[index_vec[j]] <= (*MIN_ID + ((double)(i+1))*Delta)) {
				ret[i]++;
				index_vec.erase(index_vec.begin()+j);
			}
		}
	}
	std::cout << "Histogram created with: " << std::endl;
	std::cout << "MIN = " << *MIN_ID << ", MAX = " << *MAX_ID  << ", Delta = " << Delta  << ", n = " << n << ", N = " << sum_v(ret) << std::endl;
	return ret;
}


std::vector<double> histogram_bins(const std::vector<double> *data, unsigned int n) {
	if(n == 0) {
		n = round(sqrt((*data).size()));
	}
	double min = *std::min_element((*data).begin(),(*data).end());
	double max = *std::max_element((*data).begin(),(*data).end());
	double Delta = (max - min)/n;

	std::vector<double> ret;

	unsigned int i;
	for(i=0; i < n; i++) {
		ret.push_back(min + Delta*(0.5 + (double)i));
	}
	return ret;
}

unsigned int draw_from_dist(std::vector<double> dist) {
	unsigned int i;
	double rng = (double)rand() / RAND_MAX;
	double cum = 0;
	std::vector<double> temp = dist;
	double dist_SUM = sum_v(temp);

	if(dist_SUM > (1.0+1e-6) || dist_SUM < (1.0-1e-6)) {
		for(i=0; i < temp.size(); i++) {
			temp[i] = temp[i]/dist_SUM;
		}
	}

	for(i=0; i < temp.size(); i++) {
		if((rng < (cum+temp[i])) && (rng > cum)) {
			return i;
		}
		else {
			cum += temp[i];
		}
	}

	return 0;
}

double draw_from_smoth_dist(std::vector<double> dist, std::vector<double> bins) {
	unsigned int i;
	double rng = (double)rand() / RAND_MAX;
	double rng2 = (double)rand() / RAND_MAX;
	double cum = 0;
	std::vector<double> temp = dist;
	double dist_SUM = sum_v(temp);

	if(dist_SUM > (1.0+1e-6) || dist_SUM < (1.0-1e-6)) {
		for(i=0; i < temp.size(); i++) {
			temp[i] = temp[i]/dist_SUM;
		}
	}

	for(i=0; i < temp.size(); i++) {
		if((rng < (cum+temp[i])) && (rng > cum)) {
			return (bins[i-1] + bins[i])*0.5 + 0.5*rng2*(bins[i+1] - bins[i-1]);
		}
		else {
			cum += temp[i];
		}
	}

	return 0;
}

std::vector<double> generate_N_hist_list(std::vector<std::vector<double> > sample, std::vector<std::vector<double> > bins, std::vector<double> dx) {
	unsigned int DIM,N,side_N,bin_N;

	unsigned int i,j,k;

	std::vector<double> HIST;

	DIM = sample[0].size();

	bin_N = bins.size();
	side_N = (unsigned int)round(pow((double)bin_N,1.0/(double)DIM));

	std::vector<std::vector<double> > NEW_SAMP;
	std::vector<unsigned int> erase_ind;
	bool in_bin;

	NEW_SAMP = sample;

	for(i=0; i < bin_N; i++) {
		N = NEW_SAMP.size();
		HIST.push_back(0.0);
		erase_ind.clear();
		for(j=0; j < N; j++) {
			in_bin = true;
			for(k=0; k < DIM; k++) {
				in_bin = in_bin && ( (sample[j][k] >= (bins[i][k] - dx[k]*0.5) ) && (sample[j][k] <= (bins[i][k] + dx[k]*0.5) ) );
			}
			if(in_bin) {
				HIST[i] = HIST[i] + 1.0;
				erase_ind.push_back(j);
			}
		}
		if(erase_ind.size() > 0) {
			for(j=0; j < erase_ind.size(); j++) {
				NEW_SAMP.erase(NEW_SAMP.begin() + erase_ind[erase_ind.size()-1-j]);
			}
		}
		
	}

	return HIST;
}

std::vector<std::vector<double> > generate_N_bin_list(std::vector<std::vector<double> > sample, unsigned int N_usr) {
	unsigned int DIM,N,bin_N;
	double side_N;
	unsigned int i,j;

	std::vector<std::vector<double> > BINS;
	std::vector<double> temp_cord;

	// each row is data point, all data points have same dimension
	DIM = sample[0].size();
	N = sample.size();

	if(N_usr != 0) {
		bin_N = N_usr;
	}
	else if(DIM <= 3) {
		bin_N = (unsigned int)round(sqrt(sqrt(N)));
		//bin_N = (unsigned int)round(pow(N,1.0/pow(2,DIM)));
	}
	else if(DIM <= 5) {
		bin_N = (unsigned int)round(sqrt(N));
		//bin_N = (unsigned int)round(pow(N,1.0/DIM));
	}
	else {
		bin_N = (unsigned int)round(sqrt(N));
		//bin_N = (unsigned int)round(pow(N,1.0/sqrt(DIM)));
	}

	bin_N = find_closest_XpowN((double)bin_N,DIM);
	side_N = pow((double)bin_N,1.0/(double)DIM);
	
	//construct hyper cube
	//find min max
	std::vector<double> low_edges;
	std::vector<double> high_edges;
	std::vector<double> dx;
	std::vector<double> temp_data;
	std::vector<unsigned int> mods;
	unsigned int ind,temp_ind;
	double temp_val;

	for(i=0; i < DIM; i++) {
		temp_data = extract_col(sample, i);
		temp_val = *std::min_element(temp_data.begin(),temp_data.end());
		low_edges.push_back(temp_val);
		
		temp_val = *std::max_element(temp_data.begin(),temp_data.end());
		high_edges.push_back(temp_val);
	}
	temp_data.clear();

	for(j=0; j < DIM; j++) {
		dx.push_back((high_edges[j] - low_edges[j])/side_N);
		mods.push_back((unsigned int)pow(side_N,(double)j));
	}

	for(i=0; i < bin_N; i++) {
		temp_cord.clear();
		for(j=0; j < DIM; j++) {
			ind = (unsigned int)floor((double)i/(double)mods[j]);
			temp_ind = ind % (unsigned int)side_N;
			temp_val = low_edges[j] + dx[j]*(0.5 + (double)( temp_ind ));
			temp_cord.push_back(temp_val);
		}
		BINS.push_back(temp_cord);
	}

	return BINS;
}

unsigned int find_closest_pow2(unsigned int n) {
	bool val_found = false;
	unsigned int ret = 0;
	double N = (double)n;
	double Q = 1;
	double next_Q = 2;
	double diff_down, diff_up;
	while(!val_found) {
		diff_down = abs(N-Q);
		diff_up = abs(N-next_Q);
		if(diff_up < diff_down) {
			Q = next_Q;
			next_Q = next_Q*2.0;
		}
		else {
			val_found = true;
			ret = (unsigned int)Q;
		}
	}

	return ret;
}

unsigned int find_closest_XpowN(double X, unsigned int n) {
	bool val_found = false;
	unsigned int ret = 0;
	double N = (double)n;
	double Q = 0;
	double next_Q = 1;
	double diff_down, diff_up;
	double P = 1.0;
	while(!val_found) {
		diff_down = abs(X-Q);
		diff_up = abs(X-next_Q);
		if(diff_up < diff_down) {
			P = P + 1.0;
			Q = next_Q;
			next_Q = pow(P,N);
		}
		else {
			val_found = true;
			ret = (unsigned int)Q;
		}
	}

	return ret;
}

std::vector<double> draw_from_N_dist(std::vector<double> dist, std::vector<std::vector<double> > bins, std::vector<double> dx) {

	unsigned int i,j;
	double rng = (double)rand() / RAND_MAX;
	double rng2;
	double cum = 0;
	std::vector<double> temp = dist;
	double dist_SUM = sum_v(temp);
	unsigned int D = dx.size();
	std::vector<double> COORD;

	if(dist_SUM > (1.0+1e-6) || dist_SUM < (1.0-1e-6)) {
		for(i=0; i < temp.size(); i++) {
			temp[i] = temp[i]/dist_SUM;
		}
	}

	for(i=0; i < temp.size(); i++) {
		if((rng < (cum+temp[i])) && (rng > cum)) {
			for(j=0; j < D; j++) {
				rng2 = (double)rand() / RAND_MAX;
				COORD.push_back( bins[i][j] - dx[j]*0.5 + rng2*dx[j] );
			}
			return COORD;
		}
		else {
			cum += temp[i];
		}
	}

	return COORD;

}

double draw_from_uniform(double a, double b) {
	return ((double)rand()/RAND_MAX)*(b - a) + a;
}

std::vector<std::vector<double> > matrix_transpose(const std::vector<std::vector<double> > M) {
	std::vector<std::vector<double> > Mt;
	std::vector<double> temp;
	std::vector<unsigned int> size;
	unsigned int i,j,max_size;
	for(i = 0; i < M.size(); i++) {
		size.push_back(M[i].size());
	}

	std::vector<unsigned int>::iterator TEMP_ID;
	TEMP_ID = std::max_element(size.begin(),size.end());
	max_size = *TEMP_ID;
	for(i = 0; i < max_size; i++) {
		temp.clear();
		for (j = 0; j < M.size(); j++) {
			if(size[j] > i) {
				temp.push_back(M[j][i]);
			}
		}
		Mt.push_back(temp);
	}

	return Mt;
}

std::vector<std::vector<double> > swap_cols(const std::vector<std::vector<double> > M, unsigned int a, unsigned int b) {
	std::vector<std::vector<double> > Mt;
	unsigned int i;
	Mt = M;

	for(i = 0; i < M.size(); i++) {
		Mt[i][b] = M[i][a];
		Mt[i][a] = M[i][b];
	}

	return Mt;
}

bool is_finite(double x) {
    return (x <= DBL_MAX && x >= -DBL_MAX); 
} 


int factorial(int n) {
	int temp = 1,i;
	for(i=n;i > 0; i--) {
		temp = temp*i;
	}
	return temp;
}

int stumpff(int k, int n, double x) {
	double temp = 0;
	int i;

	for(i=0; i < n; i++) {
		temp = temp + pow(-x,i)/factorial(k+2*i);
	}

	return temp;
}


std::vector<double> normal_distribution_density_vector(double mu, double sigma, double sigma_range, unsigned int n) {
	std::vector<double> ret_vec;
	unsigned int i;
	double min_x = mu - sigma_range*sigma;
	double max_x = mu + sigma_range*sigma;
	double step = (max_x - min_x)/n;

	for(i=0;i<n;i++){
		ret_vec.push_back(normal_distribution_density_function(mu,sigma,min_x + i*step));
	}

	return ret_vec;
}

std::vector<double> normal_distribution_bin_mid_vector(double mu, double sigma, double sigma_range, unsigned int n) {
	std::vector<double> ret_vec;
	unsigned int i;
	double min_x = mu - sigma_range*sigma;
	double max_x = mu + sigma_range*sigma;
	double step = (max_x - min_x)/n;

	for(i=0;i<n;i++){
		ret_vec.push_back(min_x + i*step);
	}

	return ret_vec;
}

std::vector<double> normal_distribution_density_vector_above_zero(double mu, double sigma, double sigma_range, unsigned int n) {
	std::vector<double> ret_vec;
	unsigned int i;
	double min_x = mu - sigma_range*sigma;
	if(min_x < 0) {
		min_x = 0;
	}

	double max_x = mu + sigma_range*sigma;
	double step = (max_x - min_x)/n;

	for(i=0;i<n;i++){
		ret_vec.push_back(normal_distribution_density_function(mu,sigma,min_x + i*step));
	}

	return ret_vec;
}

std::vector<double> normal_distribution_bin_mid_vector_above_zero(double mu, double sigma, double sigma_range, unsigned int n) {
	std::vector<double> ret_vec;
	unsigned int i;
	double min_x = mu - sigma_range*sigma;
	if(min_x < 0) {
		min_x = 0;
	}
	double max_x = mu + sigma_range*sigma;
	double step = (max_x - min_x)/n;

	for(i=0;i<n;i++){
		ret_vec.push_back(min_x + i*step);
	}

	return ret_vec;
}

double normal_distribution_density_function(double mu, double sigma, double x) {
	return (1/(sigma*sqrt(2*PI))*exp(-pow(x-mu,2)/(2*pow(sigma,2))));
}
