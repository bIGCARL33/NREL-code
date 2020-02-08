all:
	make clean
	make build
	make pipe
	make run

build:
	gcc spa_tester.c spa.c sampa.c bird.c -lm -o spa_tester

clean:
	rm -rf spa_tester
	rm -rf *.out

pipe: #display in output file
	./spa_tester >output.out

run: #display in terminal
	./spa_tester