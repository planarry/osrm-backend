#ifndef COMMON_HEADER_H
#define COMMON_HEADER_H

#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>
#include <list>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <map>
#include <utility>
#include <ostream>
#include <memory>

namespace ublas = boost::numeric::ublas;

typedef unsigned PointID;
typedef unsigned CoreID;
typedef unsigned ChainID;
typedef std::pair<double, double> Coordinate;
typedef std::vector<ChainID>::iterator TourBound;
typedef std::shared_ptr<std::pair<TourBound, TourBound>> TourBoundsPair;

const int UNREACHED_WEIGHT = std::numeric_limits<int>::max();
const unsigned SPECIAL_ID = std::numeric_limits<unsigned>::max();
const unsigned NEAREST_RADIUS = 2;
const unsigned MIN_CORE_RADIUS = 3;
const unsigned MOVING_MEAN_RADIUS = 0;
const float MOVING_MEAN_FACTOR = 1.9;
const unsigned MAX_TOUR_LENGTH = 20;

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
    return out;
}

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::list<T>& v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
    return out;
}


class ChainBase;
class Core;

#endif