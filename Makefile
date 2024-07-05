
CC=gcc
CFLAGS=-std=c99 -O3 -fno-math-errno -flto
WARNINGS=-Wall -Wextra -pedantic -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wconversion -Wredundant-decls -Wmissing-declarations

BIN=bin
SRC=.

ifneq ($(shell uname -s),Darwin)
  CFLAGS += -s -O0
  WARNINGS += -Wunsuffixed-float-constants
endif

%.o: %.c
	$(CC) $(CFLAGS) -MT $@ -MMD -MP -c -o $@ $< $(WARNINGS)

# all: nm-abs nm-ackley nm-hartmann3 nm-hartmann6 nm-rosenbrock nm-himmelblau nm-beale nm-sphere sh-chmod nm-magnetic objclean 
all: nm-magnetic sh-chmod objclean 

nm-%: $(SRC)/%.o $(SRC)/nelder_mead.o $(SRC)/main.o
	@mkdir -p $(BIN)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ -lm

test: all
	# @./run.sh -b nm-abs -p -7.60
	# @./run.sh -b nm-ackley -p -2.10,-3.04,4.50 
	# @./run.sh -b nm-sphere -p -2.10,-3.04,4.50
	# @./run.sh -b nm-hartmann3 -p 0.5,0.5,0.5
	# @#./run.sh -b nm-hartmann6 -p 0.5,0.5,0.5,0.5,0.5,0.5
	# @./run.sh -b nm-rosenbrock -p 1.0,0.0
	# @./run.sh -b nm-rosenbrock -p -1.0,0.0
	# @./run.sh -b nm-rosenbrock -p -1.0,1.0
	# @./run.sh -b nm-himmelblau -p 3.0,3.0
	# @./run.sh -b nm-himmelblau -p 3.0,-3.0  
	# @./run.sh -b nm-himmelblau -p -3.0,3.0
	# @./run.sh -b nm-himmelblau -p -3.0,-3.0
	# @./run.sh -b nm-beale -p 0.0,0.0
	# @./run.sh -b nm-magnetic -p 0.16,0.15,0.2 -x 1e-4 -f 1e-4
	@./run.sh -b nm-magnetic -p 0.16,0.15,0.2


sh-chmod:
	@chmod +x *.sh

clean: objclean
	@rm -rf $(BIN)

objclean:
	@rm -f $(SRC)/*.o
	@rm -f $(SRC)/*.d

ctags:
	@/usr/bin/ctags $(SRC)/*.h $(SRC)/*.c

.PHONY: test ctags clean objclean

-include $(SRC)/*.d
