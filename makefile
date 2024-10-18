CXX:=g++-13
CC:=gcc-13
CFLAGS:= -march=native -O3 -fopenmp -DGSL -DHAVE_INLINE

IBASIC:=-isystem/usr/local/include/
LINK:=-L/usr/local/lib/
LIBS:=-lm -lgsl -lcblas

SRC:= src/
SRC_DIR:= $(addsuffix objs/, $(SRC))
TEST_SRC:=tests/

IDIR:=$(IBASIC) -I$(SRC)

APV:= apv_ode

OBJS:= $(addprefix $(SRC_DIR), $(addsuffix .o, $(APV)))

# rules
$(SRC_DIR)%.o: $(SRC)%.cc | $(SRC_DIR)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) -c $< -o $@

apv_test: $(TEST_SRC)apv_test.cc $(OBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

$(SRC_DIR):
	mkdir -p $@

clean:
	rm -f *_test

clean_objs:
	rm -f $(SRC_DIR)*.o

clean_all: clean clean_objs

clean_dat:
	rm -f ./dat/*.apvdat
