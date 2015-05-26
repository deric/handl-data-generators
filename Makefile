CC = g++ -Wno-deprecated -lm -Wall -pedantic -O3
LD = g++ -Wno-deprecated
TARGET = mult_generator
CFLAGS = -lm

all: mult_generator elly

mult_generator : mult-generation.o gasdev.o random.o random_data.o
	$(LD) $(CFLAGS) -o mult_generator mult-generation.o gasdev.o random.o random_data.o

mult-generation.o : mult-generation.C
	$(CC) $(CFLAGS) -c mult-generation.C
gasdev.o : gasdev.C
	$(CC) $(CFLAGS) -c gasdev.C
random.o : random.C
	$(CC) $(CFLAGS) -c random.C
random_data.o : random_data.C
	$(CC) $(CFLAGS) -c random_data.C

elly : ellipsoid.o
	$(LD) $(CFLAGS) -o elly ellipsoid.o
elly.o : ellipsoid.C
	$(CC) $(CFLAGS) -c ellipsoid.C

clean:
	@/bin/rm *.o
