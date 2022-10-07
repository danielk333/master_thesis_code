
#include "DATA_TYPE.hh"

#include <cstddef>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <math.h>

#include <initializer_list>

int main (int argc, char *argv[]) {

    int D[] = {9,3};
    tensor<double> q(D,2,2.0);
    std::vector<std::vector<double> > Q;


    unsigned int i,j,N=1e6;

    Q.resize(9);
    for(i = 0; i < 9; i++) {
        Q[i].resize(3,2.0);
    }

    double time_ex;
    //Look at the time!
    clock_t EXEC_TIMES_tick;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        q(rand() % 9,2) = 3.0;
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    std::cout << "Get data tensor " << N << " X q(rand() % 9,2): " << time_ex/(double)N << " sec" << std::endl << std::endl;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        Q[rand() % 9][2] = 3.0;
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    std::cout << "Get data vector " << N << " X Q[rand() % 9][2]: " << time_ex/(double)N << " sec" << std::endl << std::endl;



    return 0;
}

/*
void construct_tens_1() {
    tensor<double> a;
}

void construct_vec_1() {
    std::vector<double> a;
}

void construct_tens_2() {
    int n[] = {100,100};
    tensor<double> a(n,2);
}

void construct_vec_2() {
    std::vector<std::vector<double> > a;
    a.resize(100);
    for(int i=0; i < 100; i++) {
        a[i].resize(100);
    }
}

void construct_tens_3() {
	int n[] = {100,100};
    tensor<double> a(n,2,0.0);
}

void construct_vec_3() {
    std::vector<std::vector<double> > a;
    a.resize(100);
    for(int i=0; i < 100; i++) {
        a[i].resize(100,0.0);
    }
}

void construct_tens_4(tensor<double> &b) {
    tensor<double> a(b);
}

void construct_vec_4(std::vector<std::vector<std::vector<double> > > &b) {
    std::vector<std::vector<std::vector<double> > > a(b);
}

int main (int argc, char *argv[]) {

    std::vector<double> results;
    std::vector<std::string> str;
	//Dedicated data structure itterators
	unsigned int i,j,N=3e2;
    double time_ex,tmp_time;
	//Look at the time!
    clock_t EXEC_TIMES_tick;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        construct_tens_1();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Constructor a()");
    std::cout << "Constructor tensor " << N << " X a(): " << time_ex/(double)N << " sec" << std::endl << std::endl;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        construct_vec_1();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    std::cout << "Constructor vector " << N << " X a(): " << time_ex/(double)N << " sec" << std::endl << std::endl;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        construct_tens_2();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Constructor a({100,100},2)");
    std::cout << "Constructor tensor " << N << " X a({100,100},2): " << time_ex/(double)N << " sec" << std::endl << std::endl;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        construct_vec_2();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    std::cout << "Constructor vector " << N << " X a(100,100): " << time_ex/(double)N << " sec" << std::endl << std::endl;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        construct_tens_3();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Constructor a({100,100},2,0.0)");
    std::cout << "Constructor tensor " << N << " X a({100,100},2,0.0): " << time_ex/(double)N << " sec" << std::endl << std::endl;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        construct_vec_3();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    std::cout << "Constructor vector " << N << " X a(100,100) = 0.0: " << time_ex/(double)N << " sec" << std::endl << std::endl;

    int n[] = {100,100,100};
    tensor<double> bt(n,3,0.0);
    tensor<double> bt2(n,3);
    std::vector<std::vector<std::vector<double> > > bv;
    std::vector<std::vector<std::vector<double> > > bv2;
    bv.resize(100);
    for(int i=0; i < 100; i++) {
        bv[i].resize(100);
        for(int j=0; j < 100; j++) {
            bv[i][j].resize(100,0.0);
        }
    }

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        construct_tens_4(bt);
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Constructor a(b), b({100,100,100},3,0.0)");
    std::cout << "Constructor tensor " << N << " X a(b), b({100,100,100},3,0.0): " << time_ex/(double)N << " sec" << std::endl << std::endl;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        construct_vec_4(bv);
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    std::cout << "Constructor vector " << N << " X a(b), b(100,100,100) = 0.0: " << time_ex/(double)N << " sec" << std::endl << std::endl;
    


    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        bt2 = bt;
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Assign a=b, b({100,100,100},3,0.0)");
    std::cout << "Assign tensor " << N << " X a = b, b({100,100,100},3,0.0): " << time_ex/(double)N << " sec" << std::endl;
    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        bt2.copy(bt);
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    tmp_time = time_ex;
    str.push_back("Assign (same size) a.copy(b), b({100,100,100},3,0.0)");
    std::cout << "Copy tensor " << N << " X a.copy(b), b({100,100,100},3,0.0): " << time_ex/(double)N << " sec" << std::endl << std::endl;


    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        bv2 = bv;
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    results.push_back(tmp_time);
    results.back() /= time_ex;
    std::cout << "Assign vector " << N << " X a = b, b(100,100,100) = 0.0: " << time_ex/(double)N << " sec" << std::endl << std::endl;
    
    N=N*1e4;

    double x;
    int getn[] = {6,4,4};
    long int getndir = 6 + 4*100 + 4*100*100;
    long int getndir2 = bt.vec2ind(getn);
    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        x = bt[getn];
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Get element: double x = a[{6,4,4}]");
    std::cout << "Get element from tensor " << N << " X double x = a[{6,4,4}]: " << time_ex/(double)N << " sec" << std::endl;
    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        x = bt[getndir];
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    tmp_time = time_ex;
    str.push_back("Get element: double x = a[6 + 4*100 + 4*100^2]");
    std::cout << "Get element from tensor " << N << " X double x = a[6 + 4*100 + 4*100^2]: " << time_ex/(double)N << " sec" << std::endl << std::endl;


    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        x = bv[6][4][4];
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    results.push_back(tmp_time);
    results.back() /= time_ex;
    std::cout << "Get element from vector " << N << " X double x = a[6][4][4]: " << time_ex/(double)N << " sec" << std::endl << std::endl;
    N=N/1e4;

    N=N*1e4;

    
    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        bt[getndir2] = x;
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Set element: double a[6 + 4*100 + 4*100^2] = x");
    std::cout << "Set element from tensor " << N << " X double a[6 + 4*100 + 4*100^2] = x: " << time_ex/(double)N << " sec" << std::endl;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        bv[6][4][4] = x;
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    std::cout << "Set element from vector " << N << " X double a[6][4][4] = x: " << time_ex/(double)N << " sec" << std::endl << std::endl;
    
    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        getn[0] = rand() % 100;
        getn[1] = rand() % 100;
        getn[2] = rand() % 100;
        bt[getn] = rand();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Itterate random element a[{}]");
    std::cout << "Itterate random element and set from tensor " << N << " X double a[{rand_n1,rand_n2,rand_n3}] = rand: " << time_ex/(double)N << " sec" << std::endl;
    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        getndir2 = rand() % 100 + (rand() % 100)*100 + (rand() % 100)*10000;
        bt[getndir2] = rand();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    tmp_time = time_ex;
    str.push_back("Itterate random element a[]");
    std::cout << "Itterate random element and set from tensor " << N << " X double a[rand_n1 + rand_n2*100 + rand_n3*100^2] = rand: " << time_ex/(double)N << " sec" << std::endl;

    EXEC_TIMES_tick = clock();
    for(i=0; i < N; i++) {
        bv[rand() % 100][rand() % 100][rand() % 100] = rand();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    results.push_back(tmp_time);
    results.back() /= time_ex;
    std::cout << "Itterate random element and set from vector " << N << " X double a[rand_n1][rand_n2][rand_n3] = rand: " << time_ex/(double)N << " sec" << std::endl << std::endl;
    
    EXEC_TIMES_tick = clock();
    for(long int K=0; K < bt.size(); K++) {
        bt[K] = rand();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Itterate all elements (one loop) a[:]");
    std::cout << "Itterate all elements and set from tensor " << N << " X double a[:] (one loop) = rand: " << time_ex << " sec" << std::endl;
    EXEC_TIMES_tick = clock();
    for(int k0 = 0; k0 < bt.size(0); k0++) {
        getn[0] = k0;
        for(int k1 = 0; k1 < bt.size(1); k1++) {
            getn[1] = k1;
            for(int k2 = 0; k2 < bt.size(2); k2++) {
                getn[2] = k2;
                bt[getn] = rand();
            }
        }
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    tmp_time = time_ex;
    str.push_back("Itterate all elements (3 loops) a[{:}]");
    std::cout << "Itterate all elements and set from tensor " << N << " X double a[{:}] (3 loops) = rand: " << time_ex << " sec" << std::endl;

    EXEC_TIMES_tick = clock();
    for(int k0 = 0; k0 < bv.size(); k0++) {
        for(int k1 = 0; k1 < bv[k0].size(); k1++) {
            for(int k2 = 0; k2 < bv[k0][k1].size(); k2++) {
                bv[k0][k1][k2] = rand();
            }
        }
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    results.push_back(tmp_time);
    results.back() /= time_ex;
    std::cout << "Itterate random element and set from vector " << N << " X double a[:][:][:] = rand: " << time_ex << " sec" << std::endl << std::endl;
    
    double * DITER;

    EXEC_TIMES_tick = clock();
    for(DITER = bt.begin(); DITER <= bt.end(); DITER++) {
        *DITER = (double)rand();
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Iterate all elements (1 loop iterator) a[:]");
    std::cout << "Iterate all elements and set from tensor (1 loop iterator) a[:] = rand: " << time_ex << " sec" << std::endl;

    EXEC_TIMES_tick = clock();
    for(int k0 = 0; k0 < bv.size(); k0++) {
        for(int k1 = 0; k1 < bv[k0].size(); k1++) {
            for(int k2 = 0; k2 < bv[k0][k1].size(); k2++) {
                bv[k0][k1][k2] = rand();
            }
        }
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    std::cout << "Iterate all elements and set from vector a[:][:][:] = rand: " << time_ex << " sec" << std::endl << std::endl;
    
    
    EXEC_TIMES_tick = clock();
    for(int k = 0; k < N; k++) {
        DITER = bt.begin() + bt.dim_C(0)*5 + bt.dim_C(1)*6;
        for(int k2 = 0; k2 < bt.size(2); k2++) {
            *DITER = (double)rand();
            DITER+=bt.dim_C(2);
        }
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Iterate one depth elements (1 loop iterator) a[{5,6,:}]");
    std::cout << "Iterate one depth elements and set from tensor (1 loop iterator) a[:] = rand: " << time_ex/(double)N << " sec" << std::endl;

    EXEC_TIMES_tick = clock();
    for(int k = 0; k < N; k++) {
        for(int k2 = 0; k2 < bv[5][6].size(); k2++) {
            bv[5][6][k2] = rand();
        }
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    std::cout << "Iterate one depth elements and set from vector a[5][6][:] = rand: " << time_ex/(double)N << " sec" << std::endl << std::endl;
    
    EXEC_TIMES_tick = clock();
    for(int k = 0; k < N; k++) {
        for(int k2 = 0; k2 < bt.size(2); k2++) {
            bt(5,6,k2) = (double)rand();
        }
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.push_back(time_ex);
    str.push_back("Iterate one depth elements (1 loop iterator) a(5,6,k2)");
    std::cout << "Iterate one depth elements and set from tensor (1 loop iterator) a(5,6,:) = rand: " << time_ex/(double)N << " sec" << std::endl;

    EXEC_TIMES_tick = clock();
    for(int k = 0; k < N; k++) {
        for(int k2 = 0; k2 < bv[5][6].size(); k2++) {
            bv[5][6][k2] = rand();
        }
    }
    time_ex = (double)(clock() - EXEC_TIMES_tick)/CLOCKS_PER_SEC;
    results.back() /= time_ex;
    std::cout << "Iterate one depth elements and set from vector a[5][6][:] = rand: " << time_ex/(double)N << " sec" << std::endl << std::endl;
    

    //example:
    //ITERATE TROUGH DEPTH 12 at ROW 5 COLUMN 6
    
    //DITER = bt.begin() + bt.dim_C(0)*5 + bt.dim_C(1)*6 + bt.dim_C(2)*0;
    //OR
    //DITER = bt.begin() + bt.vec2ind({5,6,0});
    //for(int k = 0; k < bv.size(2); k++) {
    //    *DITER = (double)rand();
    //    DITER+=bt.dim_C(2);
    //}

    

    N=N/1e4;

    std::cout << "RESULTS: " << std::endl;
    for(int k = 0; k < results.size(); k++) {
        std::cout << str[k] << ", Tensor " << results[k]*100 << " % of vector" << std::endl;
    }

	return 0;
}
*/
