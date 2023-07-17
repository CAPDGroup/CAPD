# CAPD2
Experimental repository for CAPD library

## Original source code repositories:

https://svn.capdnet.ii.uj.edu.pl/capd/  
https://svn.capdnet.ii.uj.edu.pl/capdDynSys4

## Building the library

Clone the repository:

    git clone https://github.com/AleksanderPasiut/capd2
    
Enter the repository, create the build folder, configure the repository and then build:

    cd capd2
    mkdir build
    cd build
    cmake ..
    make

Above commands will build the library only without tests or examples. Krak component also will not be built. In order to build krak, it is necessary to append `-DBUILD_KRAK=true` option when calling `cmake` command. In order to build examples, it is necessary to append `-DBUILD_EXAMPLE_EXECUTABLES=true` option when calling `cmake` command. In order to build and launch tests, look into section "Building and executing tests".

## Installing the library

In order to install the library simply call

    make install

when in `capd2/build` folder. You can specify the directory where the library should be installed with option `CMAKE_INSTALL_PREFIX` that must be set when calling `cmake` command. For Linux users: it might be necessary to use `sudo make install` instead of `make install` if you specify the installation directory that requires super user access rights.

## Building and executing tests

In order to build and execute tests, it necessary to run `cmake` command with option `-DBUILD_TEST_EXECUTABLES=true`, then build the tests and execute them with commands:

    make
    make tests
   
when in `capd2/build` directory.
