PROJECT := polyfit

CC := gcc
CPP := g++

IDIR   := $(CPPUTEST_HOME)/include
CFLAGS := -I$(IDIR) -Wall

LDIR := -L$(CPPUTEST_HOME)/lib
UNITTESTEXE := runUnitTests

LIBS := -lCppUTest -lCppUTestExt

CPP_SRC = test_runner.cpp \
          test_polyfit.cpp

OBJS := test_runner.o \
        test_polyfit.o

all: $(PROJECT)
	./$(UNITTESTEXE)

$(PROJECT):
	$(CPP) -c $(CPP_SRC) $(CFLAGS)
	$(CPP) -o $(UNITTESTEXE) $(OBJS) $(LIBS) $(LDIR)

.PHONY: clean
clean:
	$(RM) $(OBJS) $(PROJECT) $(UNITTESTEXE)


