#ifndef FLIPPER_H
#define FLIPPER_H

#include <iostream>

template <class T> class FlipFlop {
	T * x;
	T * y;
	bool which; // true -> 'point at x', false -> 'point at y'

public:
	FlipFlop(T * a, T * b);
	FlipFlop();
	void flip();
	T * focus_pointer();
	T * other_pointer();
};




template <class T> FlipFlop<T>::FlipFlop(T * a, T * b) {
	x = a;
	y = b;
	which = true;
}

template <class T> FlipFlop<T>::FlipFlop() {
  T* a = NULL;
  T* b = NULL;
  x = a;
  y = b;
  which = true;
}


template <class T> void FlipFlop<T>::flip() {
	if (which == true) {
		which = false;
	} else if (which == false) {
		which = true;
	}
}

template <class T> T * FlipFlop<T>::focus_pointer() {
	if (which == true) {
		return x;
	} else if (which == false) {
		return y;
	}
	return 0;
}

template <class T> T * FlipFlop<T>::other_pointer() {
	if (which == false) {
		return x;
	} else if (which == true) {
		return y;
	}
	return 0;
}

//
//int main () {
//	int a = 1;
//	int b = 8;
//	int* ap = &a;
//	int* bp = &b;
//
//	FlipFlop<int> ff (ap, bp);
//	std::cout << *(ff.focus_pointer()) << '\n';
//	std::cout << *(ff.other_pointer()) << '\n';
//	ff.flip();
//	std::cout << *(ff.focus_pointer()) << '\n';
//	std::cout << *(ff.other_pointer()) << '\n';
//
//	return 0;
//}
//
//
#endif
