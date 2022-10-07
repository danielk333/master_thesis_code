/*
 * classes.hh
 *
 *  Created on: Jul 1, 2015
 *      Author: dankas
 */

#include <vector>

#ifndef CLASS_HH_
#define CLASS_HH_

class output_stream {
    public:
    std::ostream& out1_;
    std::ostream& out2_;
    unsigned int type;

    output_stream(std::ostream& out1) : out1_(out1),out2_(out1) {}
    output_stream(std::ostream& out1,std::ostream& out2) : out1_(out1), out2_(out2) {}

    void endl();
};


#endif /* CLASS_HH_ */
