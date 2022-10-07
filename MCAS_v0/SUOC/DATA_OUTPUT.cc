#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "define.hh"
#include "functions.hh"


int save_mat(std::string out_file_name, std::vector<std::vector<double> > *M) {

        unsigned int i;
        unsigned int j;

        std::ofstream out_pos(out_file_name.c_str(), std::ios::app);
        out_pos << std::setprecision(10) << std::scientific;
        if(out_pos.fail()) {
            return DATA_SAVE_FAILED;
        }

        for(j = 0; j < (*M).size(); j++) {
            for(i = 0; i < (*M)[j].size(); i++) {
                out_pos << (((*M)[j])[i]);
                if(i < ((*M)[j].size() - 1)) {
                    out_pos << " ";
                }
                else if((*M).size() > 1) {
                    out_pos << "\r\n";
                }
                else if(j < ((*M).size() - 1)) {
                    out_pos << " ";
                }
            }
        }

        out_pos.close();

        return CALC_OK;
}
