#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

static double D[M]; // initialized delay line to zero
static long ptr = 0; // read-write offset

static double delayLine(double x) {

	double y = D[ptr]; //read operation
	D[ptr++] = x; // write operation
	if (ptr >= M) { ptr -= M;  }
 // wrap ptr around
	ptr %= M; // module operator syntax
	return y;

}