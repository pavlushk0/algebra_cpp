#include <iostream>
#include <tuple>
#include <cmath>
#include "algebra.h"

using namespace std;

void vec_show(const vec3_t &v) {
	cout << " " << v[_XC] << " " << v[_YC] << " " << v[_ZC];
}

float vec_lenght(vec3_t &v) {
	return sqrt(v[_XC]*v[_XC] +
				v[_YC]*v[_YC] +
				v[_ZC]*v[_ZC]);

}

void vec_normalize_self(vec3_t &v) {
	float len;

	len = vec_lenght(v);

	if (len != 0.0) {
		v[_ZC] = v[_ZC] / len;
		v[_XC] = v[_XC] / len;
		v[_YC] = v[_YC] / len;
	}
}

vec3_t vec_get_normalize(vec3_t &v) {
	vec3_t rt = v;

	vec_normalize_self(rt);

	return rt;
}

vec3_t vec_scale(vec3_t &v, float scale) {
	vec3_t rt;

	v[0] *= scale;
	v[1] *= scale;
	v[2] *= scale;

	return rt;
}

void vec_invert_self(vec3_t &v) {
	v[_XC] = -v[_XC];
	v[_YC] = -v[_YC];
	v[_ZC] = -v[_ZC];
}

vec3_t vec_get_invert(const vec3_t &v) {
	vec3_t rt = vec3_t(v);

	vec_invert_self(rt);

	return rt;
}

float vec_dot(vec3_t& a, vec3_t &b) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

vec3_t  vec_sum(vec3_t &a, vec3_t &b) {
	vec3_t rt;

	rt[0] = a[0] + b[0];
	rt[1] = a[1] + b[1];
	rt[2] = a[2] + b[2];

	return rt;
}

vec3_t vec_sub(vec3_t &a, vec3_t &b) {
	vec3_t rt;

	rt[0] = a[0] - b[0];
	rt[1] = a[1] - b[1];
	rt[2] = a[2] - b[2];

	return rt;
}

vec3_t vec_cross(vec3_t &a, vec3_t &b) {
	vec3_t rt;
	
	rt[0] = a[_YC]*b[_ZC] - a[_ZC]*b[_YC];
	rt[1] = a[_ZC]*b[_XC] - a[_XC]*b[_ZC];
	rt[2] = a[_XC]*b[_YC] - a[_YC]*b[_XC];

	return rt;
}
