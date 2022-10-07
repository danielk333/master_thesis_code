/*
 * ASSOCIATION.cc
 *
 *  Created on: Aug 20, 2015
 *      Author: dankas
 */

 #include <math.h>

#include "define.hh"
#include "functions.hh"


//build function to associate sample trough pointers to data structs


//size determines method of calc, if keep all in mem or save in between

//cross validation


void clustering_profile_sweep_file(std::vector<std::vector<double> > *association_profile, std::string file_name, output_stream *out, std::string clustering_method) {

	if(clustering_method.compare("SL") == 0) { //single linkage

	}
	else if(clustering_method.compare("CL") == 0) { //centroid linkage

	}

}


std::vector<std::vector<unsigned int> > clustering_matrix(std::vector<std::vector<double> > *d_matrix, double d_crit, std::string clustering_method) {
	unsigned int i,j;
   unsigned int objs = ((*d_matrix).size()) + 1;

   std::vector<std::vector<unsigned int> > clusters;
   std::vector<unsigned int> empty_vec;

   std::vector<node> nodes;
   unsigned int cluster_counter;
   node temp_node;
   for(i=0; i < objs; i++) {
      temp_node.set_node(d_matrix,i);
      nodes.push_back(temp_node);
   }

   cluster_counter = 0;
   for(i=0; i < objs; i++) {
      node_associate(&nodes, &cluster_counter, d_crit, i, -1);
   }
   for(i=0; i < cluster_counter; i++) {
      clusters.push_back(empty_vec);
   }
   for(i=0; i < objs; i++) {
      clusters[nodes[i].cluster_id-1].push_back(i);
   }

   return clusters;

	if(clustering_method.compare("SL") == 0) { //single linkage

	}
	else if(clustering_method.compare("CL") == 0) { //centroid linkage

	}
}

std::vector<std::vector<unsigned int> > single_linkage_clustering_matrix(std::vector<std::vector<double> > *d_matrix, double d_crit) {
   unsigned int i,j;
   unsigned int objs;
   std::vector<std::vector<unsigned int> > clusters;
   std::vector<unsigned int> empty_vec;

   std::vector<node> nodes;
   unsigned int cluster_counter;
   node temp_node;

   if(((*d_matrix).size()) == 1 && (*d_matrix)[0][0] == 0) {
    objs = 1;
    clusters.clear();
    clusters.push_back(empty_vec);
    clusters[0].push_back(0);
   }
   else {
    objs = ((*d_matrix).size()) + 1;

    for(i=0; i < objs; i++) {
      temp_node.set_node(d_matrix,i);
      nodes.push_back(temp_node);
    }

    cluster_counter = 0;
    for(i=0; i < objs; i++) {
        node_associate(&nodes, &cluster_counter, d_crit, i, -1);
    }
    for(i=0; i < cluster_counter; i++) {
        clusters.push_back(empty_vec);
    }
    for(i=0; i < objs; i++) {
        clusters[nodes[i].cluster_id-1].push_back(i);
    }
   }

   return clusters;
}

std::vector<std::vector<unsigned int> > single_linkage_clustering_files(std::string matrix_folder, std::vector<std::string> file_list, std::vector<std::vector<unsigned int> > index_matrix, double d_crit) {
   unsigned int i,j;
   std::vector<std::vector<unsigned int> > clusters;
   unsigned int objs = index_matrix.back().back();

   std::vector<unsigned int> empty_vec;

   std::vector<node> nodes;
   unsigned int cluster_counter;
   node temp_node;
   for(i=0; i < objs; i++) {
      temp_node.id = i;
      temp_node.member_of_cluster = FALSE;
      nodes.push_back(temp_node);
   }

   cluster_counter = 0;
   for(i=0; i < objs; i++) {
      node_associate_file(&nodes, &cluster_counter, d_crit, i, -1, matrix_folder, file_list, index_matrix);
   }
   for(i=0; i < cluster_counter; i++) {
      clusters.push_back(empty_vec);
   }
   for(i=0; i < objs; i++) {
      clusters[nodes[i].cluster_id-1].push_back(i);
   }

   return clusters;
}

int node_associate_file(std::vector<node> *nodes, unsigned int *cluster_counter, double d_crit, unsigned int node_id, int mother_cluster, std::string matrix_folder, std::vector<std::string> file_list, const std::vector<std::vector<unsigned int> > &index_matrix) {

   if((*nodes)[node_id].member_of_cluster) {
      return (*nodes)[node_id].cluster_id;
   }
   else if(mother_cluster == -1) {
      (*cluster_counter)++;
      (*nodes)[node_id].member_of_cluster = TRUE;
      (*nodes)[node_id].cluster_id = (*cluster_counter);
   }
   else if(mother_cluster >= 0) {
      (*nodes)[node_id].member_of_cluster = TRUE;
      (*nodes)[node_id].cluster_id = (unsigned int)mother_cluster;
   }
   
   unsigned int i,j,ind_d;
   int ret;
   std::vector<std::vector<double> > temp_matrix;
   std::string temp_file_path;
   double temp_d;
   ind_d = node_id - 1;
   i = 0;
   bool loop_go = TRUE;

   while(loop_go) {
      temp_matrix.clear();
      temp_file_path = matrix_folder + "/" + file_list[i];
      load_file(&temp_matrix,temp_file_path,DATA_MATRIX_SILENT);

      for(j=index_matrix[i][0]; j < index_matrix[i][1]-1; j++) {
         //std::cout << "Associating index j " << j << " of node_id " << node_id << " with temp_d " << temp_d << " index (";
         if(j < node_id) {
            temp_d = temp_matrix[j-index_matrix[i][0]][ind_d];
            //std::cout << j-index_matrix[i][0] << "," << ind_d << ")" << std::endl;
            ind_d--;
         }
         else if(j >= node_id) {
            temp_d = temp_matrix[node_id-index_matrix[i][0]][j-node_id];
            //std::cout << node_id-index_matrix[i][0] << "," << j-node_id << ")" << std::endl;
         }

         if(temp_d < d_crit) {
            //std::cout << "Neighbor found" << std::endl;
            ret = node_associate_file(nodes, cluster_counter, d_crit, j, (*nodes)[node_id].cluster_id, matrix_folder, file_list, index_matrix);
            if(ret > (*cluster_counter)) {
               std::cout << "Something went wrong in cluster association" << std::endl;
               return -1;
            }
         }
      }
      
      i++;
      if(i < index_matrix.size()) {
         if(node_id >= index_matrix[i][1]) {
            loop_go = FALSE;
         }
      }
      else {
         loop_go = FALSE;
      }
   }

   return (*nodes)[node_id].cluster_id;
}

int node_associate(std::vector<node> *nodes, unsigned int *cluster_counter, double d_crit, unsigned int node_id, int mother_cluster) {
   if((*nodes)[node_id].member_of_cluster) {
      return (*nodes)[node_id].cluster_id;
   }
   else if(mother_cluster == -1) {
      (*cluster_counter)++;
      (*nodes)[node_id].member_of_cluster = TRUE;
      (*nodes)[node_id].cluster_id = (*cluster_counter);
   }
   else if(mother_cluster >= 0) {
      (*nodes)[node_id].member_of_cluster = TRUE;
      (*nodes)[node_id].cluster_id = (unsigned int)mother_cluster;
   }

   unsigned int i;
   int ret;

   for(i=0; i < (*nodes)[node_id].d_vec.size(); i++) {
      if((*nodes)[node_id].d_vec[i] < d_crit) {
         ret = node_associate(nodes, cluster_counter, d_crit, i, (*nodes)[node_id].cluster_id);
         if(ret > (*cluster_counter)) {
            std::cout << "Something went wrong in cluster association" << std::endl;
            return -1;
         }
      }
   }

   return (*nodes)[node_id].cluster_id;
}

/*

            out << "# -- Creating database association profile " << MURMHED_profile_files[j]; out.endl();
            list_of_batches_d_matrix.clear();
            list_of_batches_d_matrix = list_files(MURMHED_D_mat_files[j]);

            index_matrix.clear();
            for(k=0; k < list_of_batches_d_matrix.size(); k++) {
                list_of_batches_d_matrix[k].erase(0,MURMHED_D_mat_files[j].size()+1);
                out << "Index range recorded for file " << list_of_batches_d_matrix[k]; out.endl();

                char_index_1 = list_of_batches_d_matrix[k].find("_", 0);
                char_index_2 = list_of_batches_d_matrix[k].find("_", char_index_1+1);

                index_matrix.push_back(EMPTY_UINT_VEC);
                buffer_temp = list_of_batches_d_matrix[k].substr(0,char_index_1);
                index_matrix[k].push_back((unsigned int)strtod(buffer_temp.c_str(), NULL));

                buffer_temp = list_of_batches_d_matrix[k].substr(char_index_1+1,char_index_2);
                index_matrix[k].push_back((unsigned int)strtod(buffer_temp.c_str(), NULL));
            }
            out << "Completed index matrix:"; out.endl();
            print_matrix(index_matrix);
            number_of_objects_db_load = index_matrix.back().back()-1;

            out << "Running cluster analysis"; out.endl();


                min_d_crit_step = 0.00001;
    max_d_crit_step = 1;
    association_adaptive_step_threshold = 0.2;//0.01;
    step_multiplyer_inc = 1.7;
    d_crit = 0;
    d_crit_step = 0.001;
    fraction_associated = 0;
    fraction_associated_next = 0;
    
    association_vec.clear();
    out << "# Starting association loop."; out.endl();
    association_loop_counter = 0;
    fraction_associated_old = 0;
    d_crit_vec.clear();
    while(fraction_associated < 1) {
        association_loop_counter++;
        
        cluster_matrix.clear();
        fractions.clear();

        cluster_matrix = cluster_analysis_batches(MURMHED_D_mat_files[j], list_of_batches_d_matrix, index_matrix, d_crit);
        for(i=0; i < cluster_matrix.size(); i++) {
            if((double)cluster_matrix[i].size() > 1) {
                fractions.push_back((double)cluster_matrix[i].size());
            }
        }
        if(fractions.size() > 0) {
            fractions = fractions/((double)number_of_objects_db_load);
            fraction_associated = sum_v(fractions);
        }
        else {
            fraction_associated = 0;
        }
        
        if(fraction_associated - fraction_associated_old != 0) {
            d_crit_vec.push_back(d_crit);
            association_vec.push_back(fraction_associated);
        }

        fraction_associated_old = fraction_associated;
        
        d_crit_old = d_crit;
        d_crit += d_crit_step;

        visited_zero_progress_counter = 0;
        visited_max_progress_counter = 0;
        do {
            cluster_matrix_next.clear();
            fractions_next.clear();

            cluster_matrix_next = cluster_analysis_batches(MURMHED_D_mat_files[j], list_of_batches_d_matrix, index_matrix, d_crit);
            for(i=0; i < cluster_matrix_next.size(); i++) {
                if((double)cluster_matrix_next[i].size() > 1) {
                    fractions_next.push_back((double)cluster_matrix_next[i].size());
                }
            }
            fractions_next = fractions_next/((double)number_of_objects_db_load);

            fraction_associated_next = sum_v(fractions_next);
            out << "# Checking next step: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
            if((fraction_associated_next - fraction_associated) > association_adaptive_step_threshold && opt.adaptive_step) {
                visited_max_progress_counter++;
                next_step_ok = FALSE;
                d_crit = d_crit_old;
                if((d_crit_step*0.5) > min_d_crit_step) {
                    d_crit_step = d_crit_step*0.5;
                }
                else {
                    out << "# Minimum step reached, continuing: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                    next_step_ok = TRUE;
                }
                d_crit += d_crit_step;
            }
            else if((fraction_associated_next - fraction_associated) == 0 && opt.adaptive_step) {
                visited_zero_progress_counter++;
                next_step_ok = FALSE;
                d_crit = d_crit_old;
                if((d_crit_step*step_multiplyer_inc) < max_d_crit_step) {
                    d_crit_step = d_crit_step*step_multiplyer_inc;
                }
                else {
                    out << "# Maximum step reached, continuing: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                    next_step_ok = TRUE;
                }
                d_crit += d_crit_step;
            }
            else {
                next_step_ok = TRUE;
            }

            if(visited_zero_progress_counter > visited_max_progress_counter && visited_max_progress_counter > 10) {
                out << "# Cluster too small, continuing: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                next_step_ok = TRUE;
            }
        } while(!next_step_ok && fraction_associated != 1);

        out << "Sample associated " << fraction_associated*100 << " percent with " << d_crit_old << " threshold."; out.endl();
    }
    
    out << "# Association loop exited."; out.endl();
    
    out << "# Calculating new means."; out.endl();
    
    //j*2       d_crit vectors
    //j*2+1     frac vectors

    association_profile_db[j*2].insert( association_profile_db[j*2].end(), d_crit_vec.begin(), d_crit_vec.end() );
    association_profile_db[j*2+1].insert( association_profile_db[j*2+1].end(), association_vec.begin(), association_vec.end() );

    order_D_c_vec.clear();
    order_D_c_vec.resize(association_profile_db[j*2].size());

    D_c_n = 0;
    for (double_vector_iterator it = association_profile_db[j*2].begin(); it != association_profile_db[j*2].end(); ++it, ++D_c_n) {
        order_D_c_vec[D_c_n] = std::make_pair(D_c_n, it);
    }

    std::sort(order_D_c_vec.begin(), order_D_c_vec.end(), double_vec_ordering_struct());

    association_profile_db[j*2] = sort_from_ref((association_profile_db[j*2]),order_D_c_vec);
    association_profile_db[j*2+1] = sort_from_ref((association_profile_db[j*2+1]),order_D_c_vec);
    bin_member_counter = 1;
    for(i=association_profile_db[j*2].size()-1; i > 0; i--) {
        if(association_profile_db[j*2][i] == association_profile_db[j*2][i-1]) {
            bin_member_counter++;
            association_profile_db[j*2+1][i-1] = association_profile_db[j*2+1][i-1] + association_profile_db[j*2+1][i];

            association_profile_db[j*2].erase(association_profile_db[j*2].begin()+i);
            association_profile_db[j*2+1].erase(association_profile_db[j*2+1].begin()+i);
        }
        else {
            std::cout << "Merging duplicate data points to mean: bin total=" << association_profile_db[j*2+1][i-1] << ", n=" << bin_member_counter <<  ", result=" << association_profile_db[j*2+1][i-1]/((double)bin_member_counter) <<  ", d_c=" << association_profile_db[j*2][i-1] << std::endl;
            if(bin_member_counter != 1) {
                association_profile_db[j*2+1][i-1] = association_profile_db[j*2+1][i-1]/((double)bin_member_counter);
                bin_member_counter = 1;
            }
        }
    }

    out << "########################################################################"; out.endl();

       }     
    }

    out << "# Saving data."; out.endl();
for(i=0; i < FUNCTIONS.size(); i++) {
    if(!(file_exists(MURMHED_profile_files[i]))) {
    association_profile_db_save.clear();
    association_profile_db_save.push_back(association_profile_db[i*2]);
    association_profile_db_save.push_back(association_profile_db[i*2+1]);

    association_profile_db_save_means.clear();
    association_profile_db_save_means.push_back(EMPTY_DOUBLE_VEC);
    association_profile_db_save_means.push_back(EMPTY_DOUBLE_VEC);

    if(opt.d_crit_bins == 0) {
        d_crit_bin_n = round(sqrt((double)(association_profile_db_save[0].size())));
    }
    else {
        d_crit_bin_n = opt.d_crit_bins;
    }
    d_crit_bin_width = (association_profile_db_save[0].back() - association_profile_db_save[0][0])/((double)d_crit_bin_n);


    for(j=0; j < d_crit_bin_n; j++) {
        association_profile_db_save_means[0].push_back(association_profile_db_save[0][0]+(j+0.5)*d_crit_bin_width);
        association_profile_db_save_means[1].push_back(-1);
    }

    //std::cout << "PRINTING MATRIX DEBUG MODE: " << std::endl;
    //print_matrix(association_profile_db_save_means);
    //std::cout << std::endl << std::endl;

    out << "## Starting bin avrage calculation"; out.endl();
    out << "Bin number: " << d_crit_bin_n << ", max val: " << association_profile_db_save[0].back(); out.endl();
    k = 0;
    for(j=0; j < association_profile_db_save[1].size(); j++) {
        if(association_profile_db_save[0][j] > (association_profile_db_save[0][0] + (k+1)*d_crit_bin_width)) {
            out << "-- Ending bin at " << association_profile_db_save_means[1][k] << ", to " << association_profile_db_save_means[1][k]/((double)bin_member_counter) << ", using n=" << bin_member_counter; out.endl();
            association_profile_db_save_means[1][k] = association_profile_db_save_means[1][k]/((double)bin_member_counter);
            k++;

            while(association_profile_db_save[0][j] > (association_profile_db_save[0][0] + (k+1)*d_crit_bin_width)) {
                out << "-- Skipping bin at " << (association_profile_db_save[0][0] + (k+1)*d_crit_bin_width); out.endl();
                association_profile_db_save_means[1][k] = 0;
                k++;
            }
        }

        //std::cout << "current bin val: " << association_profile_db_save_means[1][k] << ", j: " << j <<  ", k: " << k << std::endl;
        if(association_profile_db_save_means[1][k] == -1) {
            bin_member_counter = 1;
            out << "-- Starting new bin at: (" << association_profile_db_save[0][j] << "," << association_profile_db_save[1][j] << "), bin center and width: (" << association_profile_db_save_means[0][k] << "," << d_crit_bin_width << ")"; out.endl();
            association_profile_db_save_means[1][k] = association_profile_db_save[1][j];
        }
        else {
            out << "-- Adding value to bin: (" << association_profile_db_save[0][j] << "," << association_profile_db_save[1][j] << "), bin value: " << association_profile_db_save_means[1][k]; out.endl();
            bin_member_counter++;
            association_profile_db_save_means[1][k] = association_profile_db_save_means[1][k] + association_profile_db_save[1][j];
        }

    }
    out << "-- Ending bin at " << association_profile_db_save_means[1][k] << ", to " << association_profile_db_save_means[1][k]/((double)bin_member_counter) << ", using n=" << bin_member_counter; out.endl();
    association_profile_db_save_means[1][k] = association_profile_db_save_means[1][k]/((double)bin_member_counter);

    for(j=0; j < (d_crit_bin_n-1); j++) {
        if(association_profile_db_save_means[1][j] == 0 && j > 0) {
            empty_bin_counter=j;
            while(association_profile_db_save_means[1][empty_bin_counter] == 0) {
                empty_bin_counter++;
            }
            out << "-- Discovered empty bins, extrapolating from " << association_profile_db_save_means[1][j-1] << " to " << association_profile_db_save_means[1][empty_bin_counter] << ", missing bins=" << empty_bin_counter-j;out.endl();
            for(k=j; k < empty_bin_counter; k++) {
                association_profile_db_save_means[1][k] = association_profile_db_save_means[1][j-1] + ( ((double)(k-j+1))/((double)(empty_bin_counter-j+1)) )*(association_profile_db_save_means[1][empty_bin_counter] - association_profile_db_save_means[1][j-1]);
                out << "-- Assigning (" << association_profile_db_save_means[0][k] << "," << association_profile_db_save_means[1][k] << "), ";
            }
            out.endl();
        }
    }

    out << "- Transposing data matrix."; out.endl();
    association_profile_db_save_means = matrix_transpose(association_profile_db_save_means);

    save_mat(MURMHED_profile_files[i], &association_profile_db_save_means);
    out << "- Saved file: " << MURMHED_profile_files[i]; out.endl();
    }
}
}


*/