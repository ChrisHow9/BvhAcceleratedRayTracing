#include <cmath>
#include <cassert>
#include <cstdlib>

#include "vec3.hpp"
#include "vec4.hpp"

/** Mat44f: 4x4 matrix with floats
 *
 * See vec2f.hpp for discussion. Similar to the implementation, the Mat44f is
 * intentionally kept simple and somewhat bare bones.
 *
 * The matrix is stored in row-major order (careful when passing it to OpenGL).
 *
 * The overloaded operator () allows access to individual elements. Example:
 *    Mat44f m = ...;
 *    float m12 = m(1,2);
 *    m(0,3) = 3.f;
 *
 * The matrix is arranged as:
 *
 *   ⎛ 0,0  0,1  0,2  0,3 ⎞
 *   ⎜ 1,0  1,1  1,2  1,3 ⎟
 *   ⎜ 2,0  2,1  2,2  2,3 ⎟
 *   ⎝ 3,0  3,1  3,2  3,3 ⎠
 */
struct Mat44f
{
	float v[16];

	constexpr
	float& operator() (std::size_t aI, std::size_t aJ) noexcept
	{
		assert( aI < 4 && aJ < 4 );
		return v[aI*4 + aJ];
	}
	constexpr
	float const& operator() (std::size_t aI, std::size_t aJ) const noexcept
	{
		assert( aI < 4 && aJ < 4 );
		return v[aI*4 + aJ];
	}
};

// Identity matrix
constexpr Mat44f kIdentity44f = { {
	1.f, 0.f, 0.f, 0.f,
	0.f, 1.f, 0.f, 0.f,
	0.f, 0.f, 1.f, 0.f,
	0.f, 0.f, 0.f, 1.f
} };

// Common operators for Mat44f.
// Note that you will need to implement these yourself.

constexpr
Mat44f operator*( Mat44f const& aLeft, Mat44f const& aRight ) noexcept
{
	Mat44f kRet = {
		0.f, 0.f, 0.f, 0.f,
		0.f, 0.f, 0.f, 0.f,
		0.f, 0.f, 0.f, 0.f,
		0.f, 0.f, 0.f, 0.f
	};
	for( int i = 0; i < 4 ; i++)
	{
		for( int j = 0 ; j < 4 ; j++)
		{
			for( int k = 0 ; k  < 4 ; k++)
			{
				kRet(i,j) += aLeft(i, k) * aRight(k, j);
			}	
		}
	}
	return kRet;
}

constexpr
Vec4f operator*( Mat44f const& aLeft, Vec4f const& aRight ) noexcept
{
	Vec4f kRet = {0.f, 0.f, 0.f, 0.f};
	for (int i = 0; i < 4; i++)
	 {
        for (int j = 0; j < 4; j++)
        {
            kRet[i] += aLeft(i, j) * aRight[j];
        }
    }
	return kRet;
}

// Functions:

inline
Mat44f make_rotation_x( float aAngle ) noexcept
{
	Mat44f kRet = {
		1.f, 0.f, 0.f, 0.f,
		0.f, +cosf(aAngle), -sinf(aAngle), 0.f,
		0.f, +sinf(aAngle), +cosf(aAngle), 0.f,
		0.f, 0.f, 0.f, 1.f
	};
	return kRet;
}


inline
Mat44f make_rotation_y( float aAngle ) noexcept
{
	Mat44f kRet = {
		+cosf(aAngle), 0.f, +sinf(aAngle), 0.f,
		0.f, 		   1.f, 0.f,           0.f,
		-sinf(aAngle), 0.f, +cosf(aAngle), 0.f,
		0.f,           0.f, 0.f,           1.f
	};
	return kRet;
}

inline
Mat44f make_rotation_z( float aAngle ) noexcept
{
	Mat44f kRet = {
		+cosf(aAngle), -sinf(aAngle), 0.f, 0.f,
		+sinf(aAngle), +cosf(aAngle), 0.f, 0.f,
		0.f,           0.f,           1.f, 0.f,
		0.f,           0.f,           0.f, 1.f
	};
	return kRet;
}

inline
Mat44f make_translation( Vec3f aTranslation ) noexcept
{
	Mat44f kRet = { 
		1.f, 0.f, 0.f, aTranslation.x,
		0.f, 1.f, 0.f, aTranslation.y,
		0.f, 0.f, 1.f, aTranslation.z,
		0.f, 0.f, 0.f, 1.f
	};
	return kRet;
}

inline
Mat44f make_scaling( float aSX, float aSY, float aSZ ) noexcept
{
	Mat44f kRet = { 
		aSX, 0.f, 0.f, 0.f,
		0.f, aSY, 0.f, 0.f,
		0.f, 0.f, aSZ, 0.f,
		0.f, 0.f, 0.f, 1.f
	};
	return kRet;
}


inline
Mat44f make_perspective_projection( float aFovInRadians, float aAspect, float aNear, float aFar ) noexcept
{
	float s = 1 / tanf(aFovInRadians / 2);
	float sx = s / aAspect;
	float a = -1 * ( (aFar + aNear) / (aFar - aNear) );
	float b = -2 * ( (aFar * aNear) / (aFar - aNear) );

	Mat44f kRet = {
		sx, 0.f, 0.f, 0.f,
		0.f, s, 0.f, 0.f,
		0.f, 0.f, a, b,
		0.f, 0.f, -1.f, 1.f
	};
	return kRet;
}

constexpr
Mat44f rotation_camera(Vec3f X, Vec3f Y, Vec3f Z) noexcept
{
	return {
		X[0], X[1], X[2], 0.f,
		Y[0], Y[1], Y[2], 0.f,
		Z[0], Z[1], Z[2], 0.f,
		0.f,  0.f,  0.f,  1.f
	};
}

inline 
Mat44f lookAt( Vec3f aPosition, Vec3f aTarget, Vec3f aUp )
{
	Vec3f cameraRevDir = normalise( aPosition - aTarget );

	Vec3f cameraXaxis = normalise( cross_product( aUp, cameraRevDir ) );

	Vec3f cameraZaxis = cross_product( cameraRevDir, cameraXaxis );

	return rotation_camera( cameraXaxis, cameraZaxis, cameraRevDir ) * make_translation( -aPosition );
}
