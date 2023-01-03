.PHONY: clean
eigen: main.o eigen.o hessenberg.o schur.o matrix.o
	gcc -o eigen -Wall main.o eigen.o hessenberg.o schur.o matrix.o -lm
main.o: main.c
	gcc -c -Wall main.c
eigen.o: eigen.c
	gcc -c -Wall eigen.c
hessenberg.o: hessenberg.c
	gcc -c -Wall hessenberg.c
schur.o: schur.c
	gcc -c -Wall schur.c
matrix.o: matrix.c
	gcc -c -Wall matrix.c
clean:
	rm *.o
