#ifndef VEC3_HPP
#define VEC3_HPP

#include <cmath>
#include <cassert>
#include <cstdlib>
// code from computer graphics module
struct Vec3f
{
	float x, y, z;

	constexpr
		float& operator[] (std::size_t aI) noexcept
	{
		assert(aI < 3);
		return aI[&x]; // This is a bit sketchy, but concise and efficient.
	}
	constexpr
		float operator[] (std::size_t aI) const noexcept
	{
		assert(aI < 3);
		return aI[&x]; // This is a bit sketchy.
	}
};



float clip( float n, float lower, float upper )
{
    n = ( n > lower ) * n + !( n > lower ) * lower;
    return ( n < upper ) * n + !( n < upper ) * upper;
}

constexpr
Vec3f operator+(Vec3f aVec) noexcept
{
	return aVec;
}
constexpr
Vec3f operator-(Vec3f aVec) noexcept
{
	return { -aVec.x, -aVec.y, -aVec.z };
}

constexpr
Vec3f operator+(Vec3f aLeft, Vec3f aRight) noexcept
{
	return Vec3f{
		aLeft.x + aRight.x,
		aLeft.y + aRight.y,
		aLeft.z + aRight.z
	};
}


constexpr
bool operator==(Vec3f aLeft, Vec3f aRight) noexcept
{
	if (aLeft.x != aRight.x)
		return false;
	if (aLeft.y != aRight.y)
		return false;
	if (aLeft.z != aRight.z)
		return false;

	return true;

}
constexpr
Vec3f operator-(Vec3f aLeft, Vec3f aRight) noexcept
{
	return Vec3f{
		aLeft.x - aRight.x,
		aLeft.y - aRight.y,
		aLeft.z - aRight.z
	};
}


constexpr
Vec3f operator*(float aScalar, Vec3f aVec) noexcept
{
	return Vec3f{
		aScalar * aVec.x,
		aScalar * aVec.y,
		aScalar * aVec.z
	};
}
constexpr
Vec3f operator*(Vec3f aVec, float aScalar) noexcept
{
	return aScalar * aVec;
}

constexpr
Vec3f operator/(Vec3f aVec, float aScalar) noexcept
{
	return Vec3f{
		aVec.x / aScalar,
		aVec.y / aScalar,
		aVec.z / aScalar
	};
}

constexpr
Vec3f operator/( float aScalar,Vec3f aVec) noexcept
{
	return Vec3f{
		aScalar /aVec.x ,
		aScalar/aVec.y,
		aScalar/aVec.z
	};
}

constexpr
Vec3f& operator+=(Vec3f& aLeft, Vec3f aRight) noexcept
{
	aLeft.x += aRight.x;
	aLeft.y += aRight.y;
	aLeft.z += aRight.z;
	return aLeft;
}
constexpr
Vec3f& operator-=(Vec3f& aLeft, Vec3f aRight) noexcept
{
	aLeft.x -= aRight.x;
	aLeft.y -= aRight.y;
	aLeft.z -= aRight.z;
	return aLeft;
}

constexpr
Vec3f& operator*=(Vec3f& aLeft, float aRight) noexcept
{
	aLeft.x *= aRight;
	aLeft.y *= aRight;
	aLeft.z *= aRight;
	return aLeft;
}
constexpr
Vec3f& operator/=(Vec3f& aLeft, float aRight) noexcept
{
	aLeft.x /= aRight;
	aLeft.y /= aRight;
	aLeft.z /= aRight;
	return aLeft;
}


// A few common functions:

constexpr
float dot(Vec3f aLeft, Vec3f aRight) noexcept
{
	return aLeft.x * aRight.x
		+ aLeft.y * aRight.y
		+ aLeft.z * aRight.z
		;
}

constexpr
Vec3f cross_product(Vec3f aLeft, Vec3f aRight)
{
	return{
		aLeft[1] * aRight[2] - aLeft[2] * aRight[1],
		aLeft[2] * aRight[0] - aLeft[0] * aRight[2],
		aLeft[0] * aRight[1] - aLeft[1] * aRight[0]
	};

}

inline
float length(Vec3f aVec) noexcept
{
	// The standard function std::sqrt() is not marked as constexpr. length()
	// calls std::sqrt() unconditionally, so length() cannot be marked
	// constexpr itself.
	return std::sqrt(dot(aVec, aVec));
}


inline
Vec3f normalise(Vec3f aVec) noexcept
{
	float l = length(aVec);
	return { aVec.x / l, aVec.y / l, aVec.z / l };
}
#endif // VEC3_HPP
