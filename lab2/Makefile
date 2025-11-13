PYTHON_VERSION=3.12

main: main.cpp
	g++ -g3 -std=c++17 "$^" -o "$@"

plots: main genplots.py
	./main
	python3 genplots.py

clean:
	rm -f ./main