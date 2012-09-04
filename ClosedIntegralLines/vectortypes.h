#include <math.h>

typedef unsigned int uint;
typedef unsigned char uchar;


struct float2
{
  float x, y;
};
struct float3
{
  float x, y, z;
};
struct float4
{
  float x, y, z, w;
};
struct double2
{
  double x, y;
};
struct double3
{
  double x, y, z;
};
struct double4
{
  double x, y, z, w;
};
struct int2
{
  int x, y;
};
struct int3
{
  int x, y, z;
};
struct int4
{
  int x, y, z, w;
};


////////////////////////////////////////////////////////////////////////////////
// constructors
////////////////////////////////////////////////////////////////////////////////
inline float2 make_float2(float r, float s)
{
	float2 o;
	o.x = r; o.y = s;
    return o;
}
inline float3 make_float3(float r, float s, float t)
{
	float3 o;
	o.x = r; o.y = s; o.z = t;
    return o;
}
inline float4 make_float4(float r, float s, float t, float u)
{
	float4 o;
	o.x = r; o.y = s; o.z = t; o.w = u;
    return o;
}

inline double2 make_double2(double r, double s)
{
	double2 o;
	o.x = r; o.y = s;
    return o;
}
inline double3 make_double3(double r, double s, double t)
{
	double3 o;
	o.x = r; o.y = s; o.z = t;
    return o;
}
inline double4 make_double4(double r, double s, double t, double u)
{
	double4 o;
	o.x = r; o.y = s; o.z = t; o.w = u;
    return o;
}

inline int2 make_int2(int r, int s)
{
	int2 o;
	o.x = r; o.y = s;
    return o;
}
inline int3 make_int3(int r, int s, int t)
{
	int3 o;
	o.x = r; o.y = s; o.z = t;
    return o;
}
inline int4 make_int4(int r, int s, int t, int u)
{
	int4 o;
	o.x = r; o.y = s; o.z = t; o.w = u;
    return o;
}





inline float2 make_float2(float s)
{
    return make_float2(s, s);
}
inline float2 make_float2(float3 a)
{
    return make_float2(a.x, a.y);
}
inline float2 make_float2(int2 a)
{
    return make_float2(float(a.x), float(a.y));
}

inline int2 make_int2(int s)
{
    return make_int2(s, s);
}
inline int2 make_int2(int3 a)
{
    return make_int2(a.x, a.y);
}
inline int2 make_int2(float2 a)
{
    return make_int2(int(a.x), int(a.y));
}

inline float3 make_float3(float s)
{
    return make_float3(s, s, s);
}
inline float3 make_float3(float2 a)
{
    return make_float3(a.x, a.y, 0.0f);
}
inline float3 make_float3(float2 a, float s)
{
    return make_float3(a.x, a.y, s);
}
inline float3 make_float3(float4 a)
{
    return make_float3(a.x, a.y, a.z);
}
inline float3 make_float3(int3 a)
{
    return make_float3(float(a.x), float(a.y), float(a.z));
}


inline int3 make_int3(int s)
{
    return make_int3(s, s, s);
}
inline int3 make_int3(int2 a)
{
    return make_int3(a.x, a.y, 0);
}
inline int3 make_int3(int2 a, int s)
{
    return make_int3(a.x, a.y, s);
}
inline int3 make_int3(float3 a)
{
    return make_int3(int(a.x), int(a.y), int(a.z));
}

inline float4 make_float4(float s)
{
    return make_float4(s, s, s, s);
}
inline float4 make_float4(float3 a)
{
    return make_float4(a.x, a.y, a.z, 0.0f);
}
inline float4 make_float4(float3 a, float w)
{
    return make_float4(a.x, a.y, a.z, w);
}
inline float4 make_float4(int4 a)
{
    return make_float4(float(a.x), float(a.y), float(a.z), float(a.w));
}

inline int4 make_int4(int s)
{
    return make_int4(s, s, s, s);
}
inline int4 make_int4(int3 a)
{
    return make_int4(a.x, a.y, a.z, 0);
}
inline int4 make_int4(int3 a, int w)
{
    return make_int4(a.x, a.y, a.z, w);
}
inline int4 make_int4(float4 a)
{
    return make_int4(int(a.x), int(a.y), int(a.z), int(a.w));
}


////////////////////////////////////////////////////////////////////////////////
// negate
////////////////////////////////////////////////////////////////////////////////

inline float2 operator-(float2 &a)
{
    return make_float2(-a.x, -a.y);
}
inline int2 operator-(int2 &a)
{
    return make_int2(-a.x, -a.y);
}
inline float3 operator-(float3 &a)
{
    return make_float3(-a.x, -a.y, -a.z);
}
inline int3 operator-(int3 &a)
{
    return make_int3(-a.x, -a.y, -a.z);
}
inline float4 operator-(float4 &a)
{
    return make_float4(-a.x, -a.y, -a.z, -a.w);
}
inline int4 operator-(int4 &a)
{
    return make_int4(-a.x, -a.y, -a.z, -a.w);
}

////////////////////////////////////////////////////////////////////////////////
// addition
////////////////////////////////////////////////////////////////////////////////

inline float2 operator+(float2 a, float2 b)
{
    return make_float2(a.x + b.x, a.y + b.y);
}
inline void operator+=(float2 &a, float2 b)
{
    a.x += b.x; a.y += b.y;
}
inline float2 operator+(float2 a, float b)
{
    return make_float2(a.x + b, a.y + b);
}
inline float2 operator+(float b, float2 a)
{
    return make_float2(a.x + b, a.y + b);
}
inline void operator+=(float2 &a, float b)
{
    a.x += b; a.y += b;
}

inline int2 operator+(int2 a, int2 b)
{
    return make_int2(a.x + b.x, a.y + b.y);
}
inline void operator+=(int2 &a, int2 b)
{
    a.x += b.x; a.y += b.y;
}
inline int2 operator+(int2 a, int b)
{
    return make_int2(a.x + b, a.y + b);
}
inline int2 operator+(int b, int2 a)
{
    return make_int2(a.x + b, a.y + b);
}
inline void operator+=(int2 &a, int b)
{
    a.x += b; a.y += b;
}


inline float3 operator+(float3 a, float3 b)
{
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline void operator+=(float3 &a, float3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline float3 operator+(float3 a, float b)
{
    return make_float3(a.x + b, a.y + b, a.z + b);
}
inline void operator+=(float3 &a, float b)
{
    a.x += b; a.y += b; a.z += b;
}

inline int3 operator+(int3 a, int3 b)
{
    return make_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline void operator+=(int3 &a, int3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline int3 operator+(int3 a, int b)
{
    return make_int3(a.x + b, a.y + b, a.z + b);
}
inline void operator+=(int3 &a, int b)
{
    a.x += b; a.y += b; a.z += b;
}

inline int3 operator+(int b, int3 a)
{
    return make_int3(a.x + b, a.y + b, a.z + b);
}
inline float3 operator+(float b, float3 a)
{
    return make_float3(a.x + b, a.y + b, a.z + b);
}

inline float4 operator+(float4 a, float4 b)
{
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline void operator+=(float4 &a, float4 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline float4 operator+(float4 a, float b)
{
    return make_float4(a.x + b, a.y + b, a.z + b, a.w + b);
}
inline float4 operator+(float b, float4 a)
{
    return make_float4(a.x + b, a.y + b, a.z + b, a.w + b);
}
inline void operator+=(float4 &a, float b)
{
    a.x += b; a.y += b; a.z += b; a.w += b;
}

inline int4 operator+(int4 a, int4 b)
{
    return make_int4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline void operator+=(int4 &a, int4 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline int4 operator+(int4 a, int b)
{
    return make_int4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline int4 operator+(int b, int4 a)
{
    return make_int4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline void operator+=(int4 &a, int b)
{
    a.x += b; a.y += b; a.z += b; a.w += b;
}


////////////////////////////////////////////////////////////////////////////////
// subtract
////////////////////////////////////////////////////////////////////////////////

inline float2 operator-(float2 a, float2 b)
{
    return make_float2(a.x - b.x, a.y - b.y);
}
inline void operator-=(float2 &a, float2 b)
{
    a.x -= b.x; a.y -= b.y;
}
inline float2 operator-(float2 a, float b)
{
    return make_float2(a.x - b, a.y - b);
}
inline float2 operator-(float b, float2 a)
{
    return make_float2(b - a.x, b - a.y);
}
inline void operator-=(float2 &a, float b)
{
    a.x -= b; a.y -= b;
}

inline int2 operator-(int2 a, int2 b)
{
    return make_int2(a.x - b.x, a.y - b.y);
}
inline void operator-=(int2 &a, int2 b)
{
    a.x -= b.x; a.y -= b.y;
}
inline int2 operator-(int2 a, int b)
{
    return make_int2(a.x - b, a.y - b);
}
inline int2 operator-(int b, int2 a)
{
    return make_int2(b - a.x, b - a.y);
}
inline void operator-=(int2 &a, int b)
{
    a.x -= b; a.y -= b;
}

inline float3 operator-(float3 a, float3 b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline void operator-=(float3 &a, float3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline float3 operator-(float3 a, float b)
{
    return make_float3(a.x - b, a.y - b, a.z - b);
}
inline float3 operator-(float b, float3 a)
{
    return make_float3(b - a.x, b - a.y, b - a.z);
}
inline void operator-=(float3 &a, float b)
{
    a.x -= b; a.y -= b; a.z -= b;
}

inline int3 operator-(int3 a, int3 b)
{
    return make_int3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline void operator-=(int3 &a, int3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline int3 operator-(int3 a, int b)
{
    return make_int3(a.x - b, a.y - b, a.z - b);
}
inline int3 operator-(int b, int3 a)
{
    return make_int3(b - a.x, b - a.y, b - a.z);
}
inline void operator-=(int3 &a, int b)
{
    a.x -= b; a.y -= b; a.z -= b;
}

inline float4 operator-(float4 a, float4 b)
{
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline void operator-=(float4 &a, float4 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline float4 operator-(float4 a, float b)
{
    return make_float4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline void operator-=(float4 &a, float b)
{
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}

inline int4 operator-(int4 a, int4 b)
{
    return make_int4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline void operator-=(int4 &a, int4 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline int4 operator-(int4 a, int b)
{
    return make_int4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline int4 operator-(int b, int4 a)
{
    return make_int4(b - a.x, b - a.y, b - a.z, b - a.w);
}
inline void operator-=(int4 &a, int b)
{
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}


////////////////////////////////////////////////////////////////////////////////
// multiply
////////////////////////////////////////////////////////////////////////////////

inline float2 operator*(float2 a, float2 b)
{
    return make_float2(a.x * b.x, a.y * b.y);
}
inline void operator*=(float2 &a, float2 b)
{
    a.x *= b.x; a.y *= b.y;
}
inline float2 operator*(float2 a, float b)
{
    return make_float2(a.x * b, a.y * b);
}
inline float2 operator*(float b, float2 a)
{
    return make_float2(b * a.x, b * a.y);
}
inline void operator*=(float2 &a, float b)
{
    a.x *= b; a.y *= b;
}

inline int2 operator*(int2 a, int2 b)
{
    return make_int2(a.x * b.x, a.y * b.y);
}
inline void operator*=(int2 &a, int2 b)
{
    a.x *= b.x; a.y *= b.y;
}
inline int2 operator*(int2 a, int b)
{
    return make_int2(a.x * b, a.y * b);
}
inline int2 operator*(int b, int2 a)
{
    return make_int2(b * a.x, b * a.y);
}
inline void operator*=(int2 &a, int b)
{
    a.x *= b; a.y *= b;
}

inline float3 operator*(float3 a, float3 b)
{
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline void operator*=(float3 &a, float3 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline float3 operator*(float3 a, float b)
{
    return make_float3(a.x * b, a.y * b, a.z * b);
}
inline float3 operator*(float b, float3 a)
{
    return make_float3(b * a.x, b * a.y, b * a.z);
}
inline void operator*=(float3 &a, float b)
{
    a.x *= b; a.y *= b; a.z *= b;
}

inline int3 operator*(int3 a, int3 b)
{
    return make_int3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline void operator*=(int3 &a, int3 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline int3 operator*(int3 a, int b)
{
    return make_int3(a.x * b, a.y * b, a.z * b);
}
inline int3 operator*(int b, int3 a)
{
    return make_int3(b * a.x, b * a.y, b * a.z);
}
inline void operator*=(int3 &a, int b)
{
    a.x *= b; a.y *= b; a.z *= b;
}


inline float4 operator*(float4 a, float4 b)
{
    return make_float4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline void operator*=(float4 &a, float4 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline float4 operator*(float4 a, float b)
{
    return make_float4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline float4 operator*(float b, float4 a)
{
    return make_float4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline void operator*=(float4 &a, float b)
{
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

inline int4 operator*(int4 a, int4 b)
{
    return make_int4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline void operator*=(int4 &a, int4 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline int4 operator*(int4 a, int b)
{
    return make_int4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline int4 operator*(int b, int4 a)
{
    return make_int4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline void operator*=(int4 &a, int b)
{
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}


////////////////////////////////////////////////////////////////////////////////
// divide
////////////////////////////////////////////////////////////////////////////////

inline float2 operator/(float2 a, float2 b)
{
    return make_float2(a.x / b.x, a.y / b.y);
}
inline void operator/=(float2 &a, float2 b)
{
    a.x /= b.x; a.y /= b.y;
}
inline float2 operator/(float2 a, float b)
{
    return make_float2(a.x / b, a.y / b);
}
inline void operator/=(float2 &a, float b)
{
    a.x /= b; a.y /= b;
}
inline float2 operator/(float b, float2 a)
{
    return make_float2(b / a.x, b / a.y);
}

inline float3 operator/(float3 a, float3 b)
{
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline void operator/=(float3 &a, float3 b)
{
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
}
inline float3 operator/(float3 a, float b)
{
    return make_float3(a.x / b, a.y / b, a.z / b);
}
inline void operator/=(float3 &a, float b)
{
    a.x /= b; a.y /= b; a.z /= b;
}
inline float3 operator/(float b, float3 a)
{
    return make_float3(b / a.x, b / a.y, b / a.z);
}

inline float4 operator/(float4 a, float4 b)
{
    return make_float4(a.x / b.x, a.y / b.y, a.z / b.z,  a.w / b.w);
}
inline void operator/=(float4 &a, float4 b)
{
    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;
}
inline float4 operator/(float4 a, float b)
{
    return make_float4(a.x / b, a.y / b, a.z / b,  a.w / b);
}
inline void operator/=(float4 &a, float b)
{
    a.x /= b; a.y /= b; a.z /= b; a.w /= b;
}
inline float4 operator/(float b, float4 a){
    return make_float4(b / a.x, b / a.y, b / a.z, b / a.w);
}

////////////////////////////////////////////////////////////////////////////////
// min
////////////////////////////////////////////////////////////////////////////////
/*
inline int2 min(int2 a, int2 b)
{
    return make_int2(min(a.x,b.x), min(a.y,b.y));
}
inline int3 min(int3 a, int3 b)
{
    return make_int3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}
inline int4 min(int4 a, int4 b)
{
    return make_int4(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z), min(a.w,b.w));
}


////////////////////////////////////////////////////////////////////////////////
// max
////////////////////////////////////////////////////////////////////////////////


inline int2 max(int2 a, int2 b)
{
    return make_int2(max(a.x,b.x), max(a.y,b.y));
}
inline int4 max(int4 a, int4 b)
{
    return make_int4(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z), max(a.w,b.w));
}*/

////////////////////////////////////////////////////////////////////////////////
// dot product
////////////////////////////////////////////////////////////////////////////////

inline float dot(float2 a, float2 b)
{ 
    return a.x * b.x + a.y * b.y;
}
inline float dot(float3 a, float3 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline float dot(float4 a, float4 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline int dot(int2 a, int2 b)
{ 
    return a.x * b.x + a.y * b.y;
}
inline int dot(int3 a, int3 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline int dot(int4 a, int4 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}


////////////////////////////////////////////////////////////////////////////////
// length
////////////////////////////////////////////////////////////////////////////////

inline float length(float2 v)
{
    return sqrtf(dot(v, v));
}
inline float length(float3 v)
{
    return sqrtf(dot(v, v));
}
inline float length(float4 v)
{
    return sqrtf(dot(v, v));
}

////////////////////////////////////////////////////////////////////////////////
// normalize
////////////////////////////////////////////////////////////////////////////////

inline float2 normalize(float2 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}
inline float3 normalize(float3 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}
inline float4 normalize(float4 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}
