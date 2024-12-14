#========== you shouldn't modify anything below ===========================
export RTS = .
include rts.obj

.SUFFIXES:
.SUFFIXES: .c .h .hpp .o

default:
	@echo
	@echo "C-XSC - C++ library for eXtended Scientific Computation"
	@echo "Please use the Makefile in the upper directory!"
	@echo

all: .depend lib

lib: $(OBJ_RTS) 

FORCE: # I need the ar done !

%.o: %.c
	$(CC) $(CCFLAGS) $(FPIC) $(CCOPTIMIZE) -c $<

clean:
	-$(RM) $(OBJ_RTS) rts.a

depend .depend:
	$(CC) -MM $(CCFLAGS) $(CCOPTIMIZE) *.c > .depend
	 

.PHONY: default clean all lib


