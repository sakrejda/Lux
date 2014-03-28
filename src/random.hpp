#ifndef RANDOM_H
#define RANDOM_H

#include <vector>
#include <map>

#include <stdexcept>
#include <sstream>

class Random {

public:
    virtual void draw() = 0;
    virtual double lpdf(double X) = 0;
    virtual double lpdf() = 0;

private:

};

#endif
