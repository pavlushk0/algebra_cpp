
#ifndef __algebra_h_
#define __algebra_h_

#include <array>
#include <tuple>

using namespace std;

typedef array<float, 3> vec3_t;
typedef array<float, 9> mtrx3_t;
typedef array<float, 16> mtrx4_t;
typedef array<float, 4> qtnn_t;

const int8_t _XC  = 0;
const int8_t _YC  = 1;
const int8_t _ZC  = 2;
const int8_t _WC  = 3;

void 	vec3_show(const vec3_t v); 
vec3_t 	vec3_set(float x, float y, float z);
float 	vec3_lenght(const vec3_t v);
void 	vec3_normalize_self(vec3_t v);
vec3_t	vec3_get_normalize(const vec3_t v);
vec3_t	vec3_scale(const vec3_t v,const float scale);
void	vec3_invert_self(vec3_t v);
vec3_t	vec3_get_invert(const vec3_t v);
float	vec3_dot(const vec3_t a, const vec3_t b);
vec3_t	vec3_sum(const vec3_t a, const vec3_t b);
vec3_t	vec3_sub(const vec3_t a, const vec3_t b);
vec3_t  vec3_cross(const vec3_t a, const vec3_t b);
	
int32_t id_rw(int32_t i, int32_t j, int32_t n);
int32_t id_cw(int32_t i, int32_t j, int32_t n);

mtrx3_t mtrx3_set(mtrx3_t m);
mtrx3_t mtrx3_set_euler(const float yaw, const float pitch, const float roll);
mtrx3_t mtrx3_set_axisangl(const vec3_t ax, const float phi);
void	mtrx3_show(mtrx3_t m);
mtrx3_t mtrx3_get_idtt();
float	mtrx3_det(mtrx3_t m);
mtrx3_t mtrx3_mult(const mtrx3_t a, const mtrx3_t b);
vec3_t	mtrx3_mult_vec3(const mtrx3_t m, const vec3_t v);
//std::tuple<mtrx3_t, mtrx3_t> mtrx3_lu(const mtrx3_t m);
tuple<mtrx3_t, mtrx3_t> mtrx3_lu(const mtrx3_t m);
tuple<mtrx3_t, vec3_t> mtrx3_ldlt(const mtrx3_t m);
mtrx3_t	mtrx3_get_transpose(const mtrx3_t m);
void	mtrx3_tranpose_self(mtrx3_t m);
	

#endif
