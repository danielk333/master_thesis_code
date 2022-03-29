#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "resources.hh"
#include "functions.hh"

void save_mat(std::string out_file_name, std::vector<std::vector<double> > *M, unsigned int precision, std::string delim) {

        unsigned int i;
        unsigned int j;

        std::ofstream out_file(out_file_name.c_str(), std::ios::app);
        out_file << std::setprecision(precision) << std::scientific;
        if(out_file.fail()) {
            throw("Faild to open file to write matrix");
        }

        for(j = 0; j < (*M).size(); j++) {
        	for(i = 0; i < (*M)[j].size(); i++) {
        		out_file << (((*M)[j])[i]);
                if(i < ((*M)[j].size() - 1)) {
                    out_file << delim;
                }
                else {
                    out_file << "\r\n";
                }
        	}
        }

        out_file.close();
}


//If the dimensionality is > 2 (cube or higher)
//The function will write the data as slices of the first 2 dimensions matricies, e.g.
// a 4d-tensor X of ones that is 2x2x2x2 large will be written as
// [2x2]x2x2 = [2x2]x4 = 2x8 (hypercube)
//    1 1
//    1 1 = X(:,:,1,1)
//    X(:,:,1,2)
//    X(:,:,2,1)
//    X(:,:,2,2)
void save_tensor(std::string out_file_name, tensor<double> *M, unsigned int precision, std::string delim) {

        std::ofstream out_file(out_file_name.c_str(), std::ios::app);
        out_file << std::setprecision(precision) << std::scientific;
        if(out_file.fail()) {
            throw("Faild to open file to write tensor");
        }
        int n0 = (*M).size(0);
        long int Z = (*M).size();

        double * T_ITER;
        T_ITER = (*M).begin();

        for(int i = 0; i < Z; i++) {
            out_file << (*T_ITER);
            T_ITER++;
            if((i+1) % n0 != 0) {
                out_file << delim;
            }
            else {
                out_file << "\r\n";
            }
        }


        out_file.close();
}

void merge_data(tensor<double> *collect_tens,
                          int save_cnt,
                          tensor<double> *q_vec, 
                          tensor<double> *p_vec, 
                          tensor<double> *m_vector, 
                          tensor<int> *type_vector, 
                          tensor<int> *body_id_vector,
                          double T) {

    //Loop objects over column index (or row number) size(1)
    for(int i = 0; i < (*collect_tens).size(1); i++) {
    	(*collect_tens)(0,i,save_cnt) = T;
    	(*collect_tens)(1,i,save_cnt) = (double)(*body_id_vector)(i);
    	(*collect_tens)(2,i,save_cnt) = (*q_vec)(i,0);
    	(*collect_tens)(3,i,save_cnt) = (*q_vec)(i,1);
    	(*collect_tens)(4,i,save_cnt) = (*q_vec)(i,2);
    	(*collect_tens)(5,i,save_cnt) = (*p_vec)(i,0);
    	(*collect_tens)(6,i,save_cnt) = (*p_vec)(i,1);
    	(*collect_tens)(7,i,save_cnt) = (*p_vec)(i,2);
    	(*collect_tens)(8,i,save_cnt) = (*m_vector)(i);
    	(*collect_tens)(9,i,save_cnt) = (double)(*type_vector)(i);
    }
}
