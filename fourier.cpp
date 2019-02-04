/* ---=== DISCRETE FOURIER TRANSFORM ===---
*	Author:		Colin Day
*	Object:		Given discrete (ie non-continuous, 
*				points of data) information, the DFT
*				or Discrete Fourier Transform gives back
*				the frequencies of which the wave is composed.
*
*				This code accomplishes this with a wave class
*				and various functions inside. 
*
*				main takes three arguments:
*					- inputfile:	The name of the source file for wave file data.
									Should have one floating point entry per line.
					- outputfile:	The name of the file to which the top frequencies
									will be written.
					- request:		Number of frequencies for the program to 
									give. Range between 1 and the number of total 
									data points in the wave (incluseive). Written to
									file in order of likelyhood/confidence.
*/

#include <vector>
#include <functional>
#include <iostream>
#include <complex>
#include <fstream>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <set>
#include <cstdlib>

constexpr std::complex<double> i(0, 1);
const double pi = atan(1)*4;

using std::vector;
using std::function;

class wave {
private:
	/* --== simpoints ==--
	*  Private collection of all the wave point values.
	*  Time not given explicitly but indexes are used in 
	*  place of time.
	*/
	vector<double> simpoints{};

	/* --== GET_INDEX ==--
	*  Given a value returns its index in the wave, ie
	*  finds its time. Does NOT accound for multiple points having
	*  the same value (as std::find does not appear to). See
	*  fourier_freq for usage.
	*
	*  params:
	*    - value: the value in the wave to be searched for
	*/
	double GET_INDEX(double value) {
		return std::find(simpoints.begin(), simpoints.end(), value) - simpoints.begin();
	}

	/* --== FOURIER_ELEMENT ==--
	*  Performs the appropriate summation for some element
	*  using the un-transformed data. 
	*
	*  params:
	*    - int k: the index value for the current element to be
	*             in the new wave set of data. Ranges from 0-size()-1,
	*             see fourier() for usage and ref [1] for details
	*/
	double FOURIER_ELEMENT(int k) {
		const double N = size();
		const std::complex<double> top = -2 * pi * k * i;
		std::complex<double> rhs{};
		for (double m = 0; m < size(); m++) {
			rhs += simpoints.at(m)*exp((top*m) / N);
		}
		return std::imag(rhs);
	}

public:
	/* --== size ==--
	*  Gives the number of elements in the wave
	*
	*  params: none
	*/
	unsigned int size(){ return simpoints.size(); }
	
	/* --== wave constructor: vector ==--
	*  Constructs a wave object out of a vector
	*  of doubles
	*
	*  params: 
	*    - vector<double> invec: the input vector
	*/
	wave(vector<double> invec):simpoints(invec) { }

	/* --== wave constructor: function and size ==--
	*  Constructs a wave object out of a function and
	*  the maximum of that function. Increments by integer.
	*
	*  params:
	*    - function<double(double)> func: the function that dictates 
	*									  wave behaviour
	*
	*    - size:						  the maxmimum + 1 the for the 
	*									  input function
	*/
	wave(function<double(double)> func, double size) {
		for (int i = 0; i < size; i++)
			simpoints.push_back(func(i));
	}

	/* --== pair_print ==--
	*  Prints the contents of the wave along with indecies
	*  to an output stream in a formatted way. See the output
	*  stream operator overload for more functional printing. 
	*
	*  params:
	*    - std::ostream& os: the desired output stream
	*/
	void pair_print(std::ostream& os) {
		for (int i = 0; i < size();i++) {
			os << "(" << i <<"," << simpoints.at(i) << ")" << std::endl;
		}
	}

	/* --== fourier ==--
	*  Returns the discrete fourier transform of the
	*  wave as another wave object.
	* 
	*  params: none
	*/
	wave fourier() {
		double N = size();
		function<double(double)> g = [=](double m) {
			return simpoints.at(m);
		};
		vector<double> out{};
		
		for (float k = 0; k < N; k++) {
			double rhs{};
			rhs = FOURIER_ELEMENT(k);
			out.push_back(rhs);
		}
		return wave(out);
	}

	/* --== freq ==--
	*  Gives a number of top frequencies as found by the 
	*  discrete fourier transform of the wave. Returns them
	*  in vector form. 
	*
	*  params:
	*    - unsigned int requested_num: the number of frequencies to return
	*/
	vector<double> fourier_freq(unsigned int requested_num) {
		if (requested_num == 0) 
			throw std::runtime_error("wave::fourier_freq: Cannot return empty set of frequencies, requested_num too low");
		if (requested_num > size())
			throw std::runtime_error("wave::fourier_freq: Cannot give more frequencies than there are points in the wave, try a smaller request");
		std::set<double> minimums{};
		std::set<double> ordered_points;
		wave f_wave = fourier();
		for (auto i : f_wave.simpoints) ordered_points.insert(i);
		while (minimums.size() < requested_num) {
			minimums.insert(*ordered_points.rbegin());
			ordered_points.erase(*ordered_points.rbegin());
		}
		vector<double> indices{};
		for (auto m : minimums) {
			double ind = f_wave.GET_INDEX(m);
			indices.push_back((1 - ind / size()) * 2 * pi);
		}	
		return indices;
	}
	friend std::ostream& operator<<(std::ostream& os, wave w);
};

std::ostream& operator<<(std::ostream& os, wave w) {
	for (auto el : w.simpoints)
		os << el << std::endl;
	return os;
}

int main(int argc, char** argv) {
//	function<double(double)> f = [=](double x) {return sin(2*x)+sin(4*x)+sin(8*x); };
//	wave s(f, 100);
	if (argc != 4) {
		std::cout << "Entered wrong number of arguements." << std::endl;
		std::cout << "Usage:" << std::endl;
		std::cout << "fourier.exe <inputfile.txt> <outputfile.txt> <number_of_freq>" << std::endl;
		return -1;
	}
	else {
		std::ifstream inputfile(argv[1]);
		std::ofstream outputfile(argv[2]);
		if (!(inputfile&&outputfile)) {
			std::cout << "Could not open file(s)" << std::endl;
			return -1;
		}
		auto request = strtol(argv[3], nullptr, 10);
		if ((int)argv[3] < 1) {
			std::cout << "Number of frequencies too small, need positive and above-zero number" << std::endl;
			return -1;
		}
		vector<double> input_from_text{};
		double datapoint{};
		while (inputfile >> datapoint)
			input_from_text.push_back(datapoint);
		wave given(input_from_text);
		vector<double> frequencies{};
		try {
			frequencies = given.fourier_freq(request);
		}catch(std::runtime_error e){
			std::cout << e.what();
			return -1;
		}
		for (auto freq : frequencies)
			outputfile << freq << std::endl;
		inputfile.close();
		outputfile.close();
		std::cout << "Process completed, output frequencies written in order of likelyhood in " << argv[2];
		return 0;
	}
}