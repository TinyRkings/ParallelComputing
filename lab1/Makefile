amr_csr_serial: readfile.o
	gcc readfile.o -o amr_csr_serial -lpthread -lrt -O3
readfile.o: readfile.c
	gcc -c readfile.c

clean:
	rm -rf *.o *.out amr_csr_serial

