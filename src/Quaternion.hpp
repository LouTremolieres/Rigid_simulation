#ifndef QUATERNION_H
#define QUATERNION_H

#include "Vector3.hpp"
#include "Matrix3x3.hpp"

class Quaternion {
    
private:
    float x = 1.0f;
    float y{};
    float z{};
    float w{};

public:
    Quaternion(float x, float y, float z, float w) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = w;
    }

    Quaternion() {}

    Quaternion product(Quaternion q) {
        Quaternion res;
        Vec3f w1 = Vec3f(this->y, this->z, this->w);
        Vec3f w2 = Vec3f(q.y , q.z, q.w);

        res.x = this->x * q.x - w1.dotProduct(w2);
        Vec3f newW = this->x * w2 + q.x * w1 + w1.crossProduct(w2);
        res.y = newW[0];
        res.z = newW[1];
        res.w = newW[2];

        return res;
    }

    Quaternion sum(Quaternion q) {
        Quaternion res;
        res.x = this->x + q.x;
        res.y = this->y + q.y;
        res.z = this->z + q.z;
        res.w = this->w + q.w;

        return res;
    }

    Quaternion scalarMultiply(float a) {
        Quaternion res;
        res.x = this->x*a;
        res.y = this->y*a;
        res.z = this->z*a;
        res.w = this->w*a;

        return res;
    }

    Mat3f convertToR() {
        Mat3f R;
        R.v00 = 1 - 2*z*z - 2*w*w;
        R.v01 = 2*y*z - 2*x*w;
        R.v02 = 2*y*w + 2*x*z;
        R.v10 = 2*y*z + 2*x*w;
        R.v11 = 1 - 2*y*y - 2*w*w;
        R.v12 = 2*z*w - 2*x*y;
        R.v20 = 2*y*w - 2*x*z;
        R.v21 = 2*z*w + 2*x*y;
        R.v22 = 1 - 2*y*y - 2*z*z;

        return R;
    }

    Quaternion normalize() {
        float norm = sqrtf(x*x + y*y + z*z + w*w);
        if(norm==0.0f) {
            return *this;
        }
        this->scalarMultiply(1/norm);
        return *this;
    }
};

#endif //QUATERNION_H