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

class OBSERVATION_record {
    public:
    double a,e,i,omega,Omega;
    double ra,dec,v_g,lambda;

    OBSERVATION_record() {};

    OBSERVATION_record operator=(const OBSERVATION_record& obs);

    void print(void);
};

class node {
    public:
    std::vector<double> d_vec;
    bool member_of_cluster;
    unsigned int cluster_id;
    unsigned int id;

    node() {};
    node(std::vector<std::vector<double> > *d_mat,unsigned int ID_TAG);

    void set_node(std::vector<std::vector<double> > *d_mat,unsigned int ID_TAG);
};

class statistical_set {
    public:
    std::vector<double> data;
    std::vector<unsigned int> hist;
    std::vector<double> hist_bins;
    std::vector<double> PDF;

    statistical_set() {};
    statistical_set(std::vector<double> DATA);
    statistical_set(const statistical_set& dist);

    double mean(void);
    void calculate_hist(unsigned int n);
    void calculate_hist_bins(unsigned int n);
    void calculate_PDF(void);
    void clear_data(void);

    void convert_to_PDF_and_clear(unsigned int n);

    std::vector<std::vector<double> > save_histogram_to_file();
};

/*
class position {
    public:
    double x_;
    double y_;
    double z_;

    position() {};
    position(double x,double y,double z);
    position(const position& pos);

    void update(double x,double y,double z);
    void update(const position& pos);
    void add(double dx,double dy,double dz);
    std::vector<double> get();
    double r(void);
    double theta(void);
    double phi(void);

    position operator=(const position& pos);
    position operator+(const position& pos) const;
    position operator-(const position& pos) const;
    double operator[](const int index);
};

class velocity {
    public:
    double vx_;
    double vy_;
    double vz_;

    velocity() {};
    velocity(double vx, double vy, double vz);
    velocity(const velocity& input);
    std::vector<double> get();

    double v();
    void update(double vx, double vy, double vz);
    void update(const velocity& input);
    void add(double dvx, double dvy, double dvz);

    velocity operator=(const velocity& v);
    velocity operator+(const velocity& v) const;
    velocity operator-(const velocity& v) const;
    double operator[](const int index);

};*/

#endif /* CLASS_HH_ */
