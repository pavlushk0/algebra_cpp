#pragma once

#include <array>
#include <tuple>

using namespace std;

enum {_XC, _YC, _ZC, _WC};

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

		mtrx4_t(const mtrx4_t &m):
			data{m[0],  m[1],  m[2],  m[3], 
				 m[4],  m[5],  m[6],  m[7], 
				 m[8],  m[9],  m[10], m[11],
				 m[12], m[13], m[14], m[15]} {};
		
		~mtrx4_t() {};
	
	private:
		float data[16];
};

/*
	Redefinition by std::array

typedef array<float, 3> vec3_t;
typedef array<float, 9> mtrx3_t;
typedef array<float, 16> mtrx4_t;
typedef array<float, 4> qtnn_t;
*/

const float f_eps = 0.00001f;

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
	
constexpr int32_t id_rw(int32_t i, int32_t j, int32_t n);
constexpr int32_t id_cw(int32_t i, int32_t j, int32_t n);

mtrx3_t	mtrx3_copy(const mtrx3_t &m);
mtrx3_t mtrx3_zero();
mtrx3_t mtrx3_set(float a00, float a01, float a02,
                  float a10, float a11, float a12,
                  float a20, float a21, float a22);
mtrx3_t mtrx3_set_euler(const float yaw, const float pitch, const float roll);
mtrx3_t mtrx3_set_axisangl(const vec3_t &ax, const float phi);
void	mtrx3_show(const mtrx3_t &m);
mtrx3_t mtrx3_idtt();
float	mtrx3_det(const mtrx3_t &m);
mtrx3_t mtrx3_mult(const mtrx3_t &a, const mtrx3_t &b);
vec3_t	mtrx3_mult_vec3(const mtrx3_t &m, const vec3_t &v);
tuple<mtrx3_t, mtrx3_t> mtrx3_lu(const mtrx3_t &m);
tuple<mtrx3_t, vec3_t> mtrx3_ldlt(const mtrx3_t &m);
mtrx3_t	mtrx3_get_transpose(const mtrx3_t &m);
void	mtrx3_tranpose_self(mtrx3_t &m);
mtrx3_t	mtrx3_get_inv(const mtrx3_t &m);

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
