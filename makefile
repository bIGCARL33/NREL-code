#Solar Energy

all: #all in one command
	make clean
	make build
	make pipe
	make run

build: #compiles the code
	gcc spa_tester.c spa.c bird.c -lm -o spa_tester

clean: #removes output and executables
	rm -rf spa_tester
	rm -rf *.out

debug:
	gcc -g  -std=c99 spa_tester.c spa.c bird.c -lm
	gdb a.out

pipe: #display in output file
	./spa_tester >output.out

run: #display in terminal
	./spa_tester