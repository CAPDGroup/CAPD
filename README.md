# CAPD
Official website for CAPD::DynSys library

http://capd.ii.uj.edu.pl

See also CAPD documentation 

http://capd.ii.uj.edu.pl/html/

## Older stable releases (5.3.0 an earlier)

https://sourceforge.net/projects/capd/files/

## Quick guide on how to build the library

Check system requirements

http://capd.ii.uj.edu.pl/html/capd_requirements.html

For detailed description see

http://capd.ii.uj.edu.pl/html/capd_compilation.html

Clone the repository:

    git clone https://github.com/CAPDGroup/CAPD
    
Enter the repository, create the build folder, configure the repository and then build:

    cd CAPD
    mkdir build
    cd build
    cmake ..
    make

Above commands will build the library only (without tests or examples). 

Options:

* `-DCAPD_ENABLE_MULTIPRECISION=false` - disable multiprecision support,
* `-DCAPD_BUILD_ALL=true` - include tests and example programs into build, 
* `-DCAPD_BUILD_EXAMPLES=true` - include example programs into build. In order to build and launch tests, look into section "Building and executing tests".

## Installing the library

In order to install the library simply call

    make install

when in `CAPD/build` folder. You can specify the directory where the library should be installed with option `CMAKE_INSTALL_PREFIX` that must be set when calling `cmake` command. For Linux users: it might be necessary to use `sudo make install` instead of `make install` if you specify the installation directory that requires super user access rights.

## Building and executing tests

In order to build and execute tests, it necessary to run `cmake` command with option `-DCAPD_BUILD_TESTS=true`, then build the tests and execute them with commands:

```bash
cmake .. -DCAPD_BUILD_TESTS=true
make
make test
```
   
when in `CAPD/build` directory.
