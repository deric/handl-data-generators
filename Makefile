CC = g++ -Wno-deprecated 
LD = g++ -Wno-deprecated
TARGET = mult_generator
CFLAGS = 

mult_generator : mult-generation.o gasdev.o random.o random_data.o
	$(LD) $(CFLAGS) -o mult_generator mult-generation.o gasdev.o random.o  random_data.o

mult-generation.o : mult-generation.C
	$(CC) $(CFLAGS) -c mult-generation.C
gasdev.o : gasdev.C
	$(CC) $(CFLAGS) -c gasdev.C
random.o : random.C
	$(CC) $(CFLAGS) -c random.C
random_data.o : random_Data.C
	$(CC) $(CFLAGS) -c random_data.C


clean: 
	@/bin/rm *.o