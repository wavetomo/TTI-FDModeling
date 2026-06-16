CC = icc

CUR_DIR = $(shell pwd)
MAKE_INSTALL_PREFIX = $(CUR_DIR)/..
INSTALL_BIN_DIR = $(MAKE_INSTALL_PREFIX)/bin
INSTALL_INCLUDE_DIR = $(MAKE_INSTALL_PREFIX)/include
INSTALL_LIB_DIR = $(MAKE_INSTALL_PREFIX)/lib

INC = -I $(CWPROOT)/include -I $(INSTALL_INCLUDE_DIR)
LIB = -L $(CWPROOT)/lib -L $(INSTALL_LIB_DIR)

LIBFLAG = -Wall -pedantic -std=c11 -qopenmp -O3 -g -fPIC
EXEFLAG = -Wall -pedantic -std=c11 -qopenmp -O3 -g -fPIE

# LIBFLAG = -Wall -std=c11 -O0 -g
# EXEFLAG = -Wall -std=c11 -O0 -g

export CC CXX MPICC MPICXX NVCC \
INSTALL_BIN_DIR INSTALL_INCLUDE_DIR INSTALL_LIB_DIR \
INC LIB \
LIBFLAG EXEFLAG

all:
	mkdir -p \
	$(INSTALL_BIN_DIR) \
	$(INSTALL_INCLUDE_DIR) \
	$(INSTALL_LIB_DIR)
	cd ./include           ; $(MAKE)
	cd ./lib               ; $(MAKE)
	cd ./main              ; $(MAKE)

.PHONY: clean
clean:
	cd ./lib               ; $(MAKE) clean
	cd ./main              ; $(MAKE) clean
	rm -rf \
	$(INSTALL_BIN_DIR) \
	$(INSTALL_INCLUDE_DIR) \
	$(INSTALL_LIB_DIR)

.PHONY: remake
remake:
	$(MAKE) clean	
	$(MAKE) -j10