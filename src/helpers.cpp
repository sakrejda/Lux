bool close_enough(double a, double b) {
	if (abs(a-b) < pow(10,-3))
		return true;
	else
		return false;
}

bool has_imaginary(std::complex<double> x) {
	return x.imag() == 0;
}

