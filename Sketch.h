#ifndef SKETCH_H
#define SKETCH_H

class Sketch {
public:
    virtual void increment_count(std::string element) = 0;
    virtual uint64_t retrieve_count(std::string element) = 0;
};

#endif
