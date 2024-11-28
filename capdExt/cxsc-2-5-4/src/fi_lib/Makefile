#========== you shouldn't modify anything below ===========================

# /*********************************************************************/
# /*       fi_lib  --- A fast interval library (Makefile)              */
# /*********************************************************************/

export FILIB = .
include fi_lib.obj

.SUFFIXES:
.SUFFIXES: .cpp .hpp .o

# Default and clean up rules
#---------------------------
default:
	@echo
	@echo "C-XSC - C++ library for eXtended Scientific Computation"
	@echo "Please use the Makefile in the upper directory!"
	@echo
	
# Delete all *.o *.a and executable files
#----------------------------------------
clean:	
	@rm -f *.a *.o


lib: $(OBJECTS)

%.o: %.cpp fi_lib.hpp
	$(CXX) -c $(CXXFLAGS) -I../rts -I../asm -I../.. $(FPIC) $(CXXOPTIMIZE) $<


all: lib 

.PHONY: default lib all

