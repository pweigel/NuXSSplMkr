#include <vector>
#include <stdexcept>

namespace nuxssplmkr {

class BilinearInterpolator {
public:
    BilinearInterpolator() {}
    BilinearInterpolator(const std::vector<double>& data, int numPointsX, int numPointsY, double minX, double maxX, double minY, double maxY)
        : data_(data), numPointsX_(numPointsX), numPointsY_(numPointsY), minX_(minX), maxX_(maxX), minY_(minY), maxY_(maxY),
          xStep_((maxX - minX) / (numPointsX - 1)), yStep_((maxY - minY) / (numPointsY - 1))
    {
        if (data_.size() != static_cast<size_t>(numPointsX * numPointsY)) {
            throw std::invalid_argument("Invalid data size");
        }
    }

    double interpolate(double x, double y) const {
        // if (x < minX_ || x > maxX_ || y < minY_ || y > maxY_) {
        //     throw std::out_of_range("Interpolation point is outside the grid limits");
        // }

        int index = linearIndex(x, y);
        double xFraction = (x - minX_) / xStep_;
        double yFraction = (y - minY_) / yStep_;

        double value11 = getValue(index);
        double value12 = getValue(index + numPointsX_);
        double value21 = getValue(index + 1);
        double value22 = getValue(index + numPointsX_ + 1);

        // Bilinear interpolation formula
        return (1 - xFraction) * ((1 - yFraction) * value11 + yFraction * value12) +
               xFraction * ((1 - yFraction) * value21 + yFraction * value22);
    }

private:
    int linearIndex(double x, double y) const {
        int xIndex = static_cast<int>((x - minX_) / xStep_);
        int yIndex = static_cast<int>((y - minY_) / yStep_);
        return yIndex * numPointsX_ + xIndex;
    }

    double getValue(int index) const {
        return data_[index];
    }

    std::vector<double> data_;
    int numPointsX_;
    int numPointsY_;
    double minX_;
    double maxX_;
    double minY_;
    double maxY_;
    double xStep_;
    double yStep_;
};

}