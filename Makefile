########################################################################
# Makefile for Quick sort (parallell algorithm)
########################################################################

CC         =  mpicc
CCFLAGS    =  -O3
LIBS       =  -lmpi

main:         	main.c
	$(CC) $(CCFLAGS) -o main main.c $(LIBS)

main1:			main1.c
	$(CC) $(CCFLAGS) -o main1 main1.c $(LIBS)