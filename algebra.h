#pragma once

#include <tuple>

using namespace std;

enum {_XC, _YC, _ZC, _WC};
enum mtrx_type {MTRX_IDTT, MTRX_ZERO};

class vec3_t {
	public:
		float operator[](const int32_t id) const {
			return data[id];
		};

		float &operator[](const int32_t id) {
			return data[id];
		};

		vec3_t(): 
			data{0.0, 0.0, 0.0} {};
		
		vec3_t(const float x, const float y, const float z): 
			data{x, y, z} {};
		
		vec3_t(const vec3_t &v): 
			data{v[_XC], v[_YC], v[_ZC]} {};

		~vec3_t() {};
	
	private:
		float data[3];
};

class mtrx3_t {
	public:
		int mrange = 3;

		float operator[](const int32_t id) const {
			return data[id];
		};

		float &operator[](const int32_t id) {
			return data[id];
		};

		mtrx3_t(): 
			data{1.0f, 0.0f, 0.0f,
				 0.0f, 1.0f, 0.0f,
				 0.0f, 0.0f, 1.0f} {};

		mtrx3_t(float a00, float a01, float a02,
                float a10, float a11, float a12,
                float a20, float a21, float a22):
			data{a00, a01, a02,
                 a10, a11, a12,
                 a20, a21, a22} {};

		mtrx3_t(const mtrx3_t &m):
			data{m[0], m[1], m[2],
				 m[3], m[4], m[5],
				 m[6], m[7], m[8]} {};
		
		~mtrx3_t() {};
	
	private:
		float data[9];
};

class qtnn_t {
	public:
		float operator[](const int32_t id) const {
			return data[id];
		};

		float &operator[](const int32_t id) {
			return data[id];
		};

		qtnn_t(): 
			data{0.0, 0.0, 0.0} {};
		
		qtnn_t(const float x, const float y, const float z, const float w): 
			data{x, y, z, w} {};
		
		qtnn_t(const qtnn_t &q): 
			data{q[_XC], q[_YC], q[_ZC], q[_WC]} {};

		~qtnn_t () {};
	
	private:
		float data[4];
};

class mtrx4_t {
	public:
		float operator[](const int32_t id) const {
			return data[id];
		};

		float &operator[](const int32_t id) {
			return data[id];
		};

		mtrx4_t(): 
			data{1.0f, 0.0f, 0.0f, 0.0f,
				 0.0f, 1.0f, 0.0f, 0.0f,
				 0.0f, 0.0f, 1.0f, 0.0f,
				 0.0f, 0.0f, 0.0f, 1.0f} {};

		mtrx4_t(float a00, float a01, float a02, float a03,
                float a10, float a11, float a12, float a13,
                float a20, float a21, float a22, float a23,
                float a30, float a31, float a32, float a33):
			data {a00, a01, a02, a03,
                  a10, a11, a12, a13,
                  a20, a21, a22, a23,
                  a30, a31, a32, a33} {};

		mtrx4_t(const mtrx4_t &m):
			data{m[0],  m[1],  m[2],  m[3], 
				 m[4],  m[5],  m[6],  m[7], 
				 m[8],  m[9],  m[10], m[11],
				 m[12], m[13], m[14], m[15]} {};
		
		~mtrx4_t() {};
	
	private:
		float data[16];
};

class vec4_t {
	public:
		float operator[](const int32_t id) const {
			return data[id];
		};

		float &operator[](const int32_t id) {
			return data[id];
		};

		vec4_t(): 
			data{0.0, 0.0, 0.0, 0.0} {};
		
		vec4_t(const float x, const float y, const float z, const float w): 
			data{x, y, z, w} {};
		
		vec4_t(const vec4_t &v): 
			data{v[_XC], v[_YC], v[_ZC], v[_WC]} {};

		~vec4_t() {};
	
	private:
		float data[4];
};

const float f_eps = 0.00001f;

/*	multidimensional array mapping, array[i][j]
	row-wise (C, C++):
	(0	1)
	(2	3)

	column-wise (Fortran, Matlab):
	(0	2)
	(1	3)
*/

constexpr int32_t id_rw(int32_t i, int32_t j, int32_t n) {
	return (i*n + j);
};

constexpr int32_t id_cw(int32_t i, int32_t j, int32_t n) {
	return (j*n + i);
};

vec3_t  vec3_copy(const vec3_t &v);
void 	vec3_show(const vec3_t &v);
vec3_t  vec3_zero();
vec3_t 	vec3_set(float x, float y, float z); 
float 	vec3_lenght(const vec3_t &v);
void 	vec3_normalize_self(vec3_t &v);
vec3_t	vec3_get_normalize(const vec3_t &v);
vec3_t	vec3_scale(const vec3_t &v,const float scale);
void	vec3_invert_self(vec3_t &v);
vec3_t	vec3_get_invert(const vec3_t v);
float	vec3_dot(const vec3_t &a, const vec3_t &b);
vec3_t	vec3_sum(const vec3_t &a, const vec3_t &b);
vec3_t	vec3_sub(const vec3_t &a, const vec3_t &b);
vec3_t  vec3_cross(const vec3_t &a, const vec3_t &b);

mtrx3_t	mtrx_copy(const mtrx3_t &m);
mtrx4_t	mtrx_copy(const mtrx4_t &m);
mtrx3_t mtrx_set_euler(const float yaw, const float pitch, const float roll);
mtrx3_t mtrx_set_axisangl(const vec3_t &ax, const float phi);
void	mtrx_show(const mtrx3_t &m);
void	mtrx_show(const mtrx4_t &m);
float	mtrx_det(const mtrx3_t &m);

template <typename mtrxT_t, int mrange>
mtrxT_t mtrx_mult(const mtrxT_t &a, const mtrxT_t &b);

vec3_t	mtrx_mult_vec3(const mtrx3_t &m, const vec3_t &v);

template <typename mtrxT_t, int mrange>
tuple<mtrxT_t, mtrxT_t> mtrx_lu(const mtrxT_t &m);

tuple<mtrx3_t, vec3_t> mtrx_ldlt(const mtrx3_t &m);
tuple<mtrx4_t, vec3_t> mtrx_ldlt(const mtrx4_t &m);
mtrx3_t	mtrx_get_transpose(const mtrx3_t &m);
mtrx4_t	mtrx_get_transpose(const mtrx4_t &m);
void	mtrx_tranpose_self(mtrx3_t &m);
void	mtrx_tranpose_self(mtrx4_t &m);
mtrx3_t	mtrx_get_inv(const mtrx3_t &m);
mtrx4_t mtrx_get_inv(const mtrx4_t &m);
mtrx4_t mtrx_get_inv_gauss(const mtrx4_t &m); //empty

qtnn_t  qtnn_zero();
qtnn_t  qtnn_copy(const qtnn_t &q);    
void 	qtnn_show(const qtnn_t &q); 
float	qtnn_lenght(const qtnn_t &q);
void	qtnn_normalize_self(qtnn_t &q);
qtnn_t	qtnn_get_normalize(const qtnn_t &q);
void	qtnn_invert_self(qtnn_t &q);
qtnn_t	qtnn_get_invert(const qtnn_t &q);
qtnn_t	qtnn_scale(const qtnn_t &q, float scale);
qtnn_t	qtnn_sum(const qtnn_t &a, const qtnn_t &b);
qtnn_t 	qtnn_sub(const qtnn_t &a, const qtnn_t &b);
float   qtnn_dot(const qtnn_t &a, const qtnn_t &b);
qtnn_t  qtnn_mult(const qtnn_t &a, const qtnn_t &b); 
qtnn_t  qtnn_mult_vec3(const qtnn_t a, const qtnn_t &b);
qtnn_t	qtnn_from_vec3(const vec3_t);
qtnn_t  qtnn_from_axisangl(const vec3_t &a, float phi); 
qtnn_t	qtnn_from_euler(float yaw, float pitch, float roll); 
vec3_t  qtnn_to_vec3(const qtnn_t &q);
vec3_t  qtnn_transform_vec3(const qtnn_t &a, const vec3_t &b);
