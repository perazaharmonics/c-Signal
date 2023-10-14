#ifndef COMPUTE_AF_H  // Include guard to prevent double inclusion
#define COMPUTE_AF_H

#include <vector>
#include <complex>

// If needed, you can include the definition of M_PI here or ensure it's defined elsewhere in your project
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 

typedef std::vector<double> Signal;
typedef std::vector<std::vector<double>> AFMatrix;  // Ensure it's double, not complex

AFMatrix computeAF(const Signal& x, const std::vector<double>& tau_values, const std::vector<double>& f_values);

// Declaration of the computeAF function
AFMatrix computeAF(const Signal& x, const std::vector<double>& tau_values, const std::vector<double>& f_values);

#endif // COMPUTE_AF_H
