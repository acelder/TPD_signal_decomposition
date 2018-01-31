

#include "unit_tests.h"
#include <iostream>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

int main( int argc , char **argv ) {

    clock_t t;
    t = clock();

    run_unit_tests();






    t = clock() - t;
    cout << "process finished in " << (float)t/CLOCKS_PER_SEC << "seconds." << endl;



    return 0;
}