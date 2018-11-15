
#include <cmath>
#include <algorithm>
#include "algebra_OO.h"

float lerp(float a, float b, float t) {
	return a + t*(b - a);
}

void swap_float(float &a, float &b) {
	float c;
	
	c = b;
	b = a;
	a = c;
}

inline int id_rw(int i, int j, int n) {
	return (i*n + j); 
}

inline int id_cw(int i, int j, int n) {
	return (j*n + i); 
}

/*------------------------------------------------------------------------*/
/*																		  */
/* Three dimensional vector definition									  */
/*																		  */
/*------------------------------------------------------------------------*/

vec3_c::~vec3_c() {

}

void vec3_c::set(const int index, const float x) {
	if ((index < 0) && (index > 3)) return;
	else data[index] = x;
}

void vec3_c::set(const vec3_c & vec) {
	this->data[_XC] = vec.data[_XC];
	this->data[_YC] = vec.data[_YC];
	this->data[_ZC] = vec.data[_ZC];
}

float vec3_c::get(const int index) const {
	if ((index < 0) && (index > 3)) return 0.0;
	else return data[index];
}

vec3_c& vec3_c::operator=(const vec3_c &v) {
	if (this == &v) {
		return (*this);
	}

    for (int i = 0; i < 3; i++) this->data[i] = v.data[i];

	return (*this);
}

float vec3_c::operator[](const int index) const {
	return get(index);
}

float& vec3_c::operator[](const int index) {
	return this->data[index];
}

vec3_c vec3_c::operator*(const float scalar) {
	return scale(scalar);
}

float vec3_c::operator*(const vec3_c &v) {
	return dot(v);
}

vec3_c vec3_c::operator^(const vec3_c &v) {
	return cross(v);
}

vec3_c vec3_c::operator+(const vec3_c &v) {
	return add(v);
}

vec3_c	vec3_c::operator-(const vec3_c &v) {
	return sub(v);
}

float vec3_c::dot(const vec3_c &b) const {
	return data[_XC]*b.data[_XC] + data[_YC]*b.data[_YC] + data[_ZC]*b.data[_ZC];
}

vec3_c vec3_c::cross(const vec3_c &b) const {
	vec3_c result;

	result.data[_XC] = data[_YC]*b.data[_ZC] - data[_ZC]*b.data[_YC];
	result.data[_YC] = data[_ZC]*b.data[_XC] - data[_XC]*b.data[_ZC];
	result.data[_ZC] = data[_XC]*b.data[_YC] - data[_YC]*b.data[_XC];

	return result;
}

float vec3_c::lenght() const {
	return sqrt(data[_XC]*data[_XC] + data[_YC]*data[_YC] + data[_ZC]*data[_ZC]);
}

vec3_c vec3_c::normalize() const {
	vec3_c ret(*this);

	float l = this->lenght();
	
	if (l != 0.0) {
		ret.data[_XC] = ret.data[_XC] / l;
		ret.data[_YC] = ret.data[_YC] / l;
		ret.data[_ZC] = ret.data[_ZC] / l;
	} else return *this;

	return ret;
}

void vec3_c::normalize_self() {
	float l = this->lenght();
	
	if (l != 0.0) {
		data[_XC] = data[_XC] / l;
		data[_YC] = data[_YC] / l;
		data[_ZC] = data[_ZC] / l;
	} else return;
}

vec3_c vec3_c::invert() const {
	vec3_c ret;

	ret.data[_XC] = -this->data[_XC];
	ret.data[_YC] = -this->data[_YC];
	ret.data[_ZC] = -this->data[_ZC];

	return ret;
}

void vec3_c::invert_self() {
	this->data[_XC] = -this->data[_XC];
	this->data[_YC] = -this->data[_YC];
	this->data[_ZC] = -this->data[_ZC];
}

vec3_c vec3_c::scale(const float &t) {
	return vec3_c(data[_XC] * t, data[_YC] * t, data[_ZC] * t);
}

vec3_c vec3_c::add(const vec3_c &b) {
	return vec3_c(data[_XC] + b.data[_XC], data[_YC] + b.data[_YC], data[_ZC] + b.data[_ZC]);
}

vec3_c vec3_c::sub(const vec3_c &b) {
	return vec3_c(data[_XC] - b.data[_XC], data[_YC] - b.data[_YC], data[_ZC] - b.data[_ZC]);
};	

/*------------------------------------------------------------------------*/
/*																		  */
/* 3x3 matrix definition		     									  */
/*																		  */
/* 0   1   2															  */
/* 3   4   5															  */
/* 6   7   8															  */
/*------------------------------------------------------------------------*/

mtrx3_c::~mtrx3_c() {
	
}

void mtrx3_c::from_euler(const float yaw, const float pitch, const float roll) {
	float cosy = cos(yaw);
	float siny = sin(yaw);
	float cosp = cos(pitch);
	float sinp = sin(pitch);
	float cosr = cos(roll);
	float sinr = sin(roll);

	data[0] =  cosy*cosr - siny*cosp*sinr;
	data[1] = -cosy*sinr - siny*cosp*cosr;
	data[2] =  siny*sinp;

	data[3] =  siny*cosr + cosy*cosp*sinr;
	data[4] = -siny*sinr + cosy*cosp*cosr;
	data[5] = -cosy*sinp;

	data[6] = sinp*sinr;
	data[7] = sinp*cosr;
	data[8] = cosp;
}

void mtrx3_c::from_axis_angl(const vec3_c &ax, const float phi) {
	float cosphi = cos(phi);
	float sinphi = sin(phi);
	float vxvy = ax[_XC]*ax[_YC];
	float vxvz = ax[_XC]*ax[_ZC];
	float vyvz = ax[_YC]*ax[_ZC];
	float vx = ax[_XC];
	float vy = ax[_YC];
	float vz = ax[_ZC];

	data[0] = cosphi + (1.0 - cosphi)*vx*vx;
	data[1] = (1.0 - cosphi)*vxvy - sinphi*vz;
	data[2] = (1.0 - cosphi)*vxvz + sinphi*vy;

	data[3] = (1.0 - cosphi)*vxvy + sinphi*vz;
	data[4] = cosphi + (1.0 - cosphi)*vy*vy;
	data[5] = (1.0 - cosphi)*vyvz - sinphi*vz;

	data[6] = (1.0 - cosphi)*vxvz - sinphi*vy;
	data[7] = (1.0 - cosphi)*vyvz + sinphi*vx;
	data[8] = cosphi + (1.0 - cosphi)*vz*vz;
}

float mtrx3_c::det() {
	return data[0]*data[4]*data[8] +
		   data[6]*data[1]*data[5] +
		   data[2]*data[3]*data[7] -
		   data[0]*data[7]*data[5] -
		   data[8]*data[3]*data[1];
}

mtrx3_c mtrx3_c::inverse() {
	return mtrx3_c();
}

void mtrx3_c::tranpose() {

}

mtrx3_c mtrx3_c::mult(const mtrx3_c &b) {
	const int _N = 3;
	mtrx3_c result;

	for (int j = 0; j < _N; j++)
		for (int i = 0; i < _N; i++)
			result.data[id_rw(i,j,_N)] = data[id_rw(0,j,_N)]*b.data[id_rw(i,0,_N)] + 
								   		 data[id_rw(1,j,_N)]*b.data[id_rw(i,1,_N)] +
								    	 data[id_rw(2,j,_N)]*b.data[id_rw(i,2,_N)];

	return result;
}

vec3_c mtrx3_c::mult_vec3(const vec3_c &b) {
	vec3_c result(data[0]*b[_XC] + data[1]*b[_YC] + data[2]*b[_ZC],
				  data[3]*b[_XC] + data[4]*b[_YC] + data[5]*b[_ZC],
				  data[6]*b[_XC] + data[7]*b[_YC] + data[8]*b[_ZC]);
	
	return result;
}

/*------------------------------------------------------------------------*/
/*																		  */
/* Quaternion definition												  */
/*																		  */
/*------------------------------------------------------------------------*/

qtnn_c::~qtnn_c() {

}

void qtnn_c::set(const qtnn_c &q) {
	this->data[_WC] = q.data[_WC];
	this->data[_XC] = q.data[_XC];
	this->data[_YC] = q.data[_YC];
	this->data[_ZC] = q.data[_ZC];
}

qtnn_c&	qtnn_c::operator=(const qtnn_c &q) {
	if (this == &q) {
		return *this;
	}

    for (int i = 0; i < 3; i++) this->data[i] = q.data[i];

	return *this;
}


float qtnn_c::lenght() const {
	return sqrt(this->data[_XC]*this->data[_XC] + 
				this->data[_YC]*this->data[_YC] + 
				this->data[_ZC]*this->data[_ZC] +
				this->data[_WC]*this->data[_WC]);
}

void qtnn_c::normalize_self() {
	float l = this->lenght();
	
	if (l != 0.0) {
		this->data[_WC] = this->data[_WC] / l;
		this->data[_XC] = this->data[_XC] / l;
		this->data[_YC] = this->data[_YC] / l;
		this->data[_ZC] = this->data[_ZC] / l;
	}
}

qtnn_c qtnn_c::invert() const {
	qtnn_c res;
	
	res.data[_WC] =  this->data[_WC];
	res.data[_XC] = -this->data[_XC];
	res.data[_YC] = -this->data[_YC];
	res.data[_ZC] = -this->data[_ZC];
	
	res.normalize_self();
	
	return res;
}

qtnn_c qtnn_c::from_mtrx3(const mtrx3_c &b) {
	return qtnn_c();
}

void qtnn_c::from_euler(const float yaw, const float pitch, const float roll) {
	qtnn_c qyaw, qpitch, qroll, tmp;
	
	qyaw.from_axis_angl(vec3_c(1.0, 0.0, 0.0), yaw);
	qpitch.from_axis_angl(vec3_c(0.0, 1.0, 0.0), pitch);
	qroll.from_axis_angl(vec3_c(0.0, 0.0, 1.0), roll);
	
	tmp = qyaw.mult(qpitch);
	
	this->set(tmp.mult(qroll));
}

void qtnn_c::from_axis_angl(const vec3_c &ax, const float phi) {
	float sinhalfphi = sin(phi * 0.5);  
	this->data[_WC] = cos(phi * 0.5);
	this->data[_XC] = ax.data[_XC] * sinhalfphi;
	this->data[_YC] = ax.data[_YC] * sinhalfphi;
	this->data[_ZC] = ax.data[_ZC] * sinhalfphi;
}

qtnn_c qtnn_c::scale(const float t) {
	qtnn_c res;
	res.data[_WC] = this->data[_WC] * t;
	res.data[_XC] = this->data[_XC] * t;
	res.data[_YC] = this->data[_YC] * t;
	res.data[_ZC] = this->data[_ZC] * t;
	return res;
}

qtnn_c qtnn_c::add(const qtnn_c &q) {
	return qtnn_c(data[_WC] + q.data[_WC], 
				  data[_XC] + q.data[_XC], 
				  data[_YC] + q.data[_YC], 
				  data[_ZC] + q.data[_ZC]);
}

qtnn_c qtnn_c::mult(const qtnn_c &b) {
	qtnn_c res;
	res.data[_WC] = data[_WC] * b.data[_WC] - data[_XC] * b.data[_XC] - data[_YC] * b.data[_YC] - data[_ZC] * b.data[_ZC];
    res.data[_XC] = data[_WC] * b.data[_XC] + data[_XC] * b.data[_WC] + data[_YC] * b.data[_ZC] - data[_ZC] * b.data[_YC];
    res.data[_YC] = data[_WC] * b.data[_YC] - data[_XC] * b.data[_ZC] + data[_YC] * b.data[_WC] + data[_ZC] * b.data[_XC];
    res.data[_ZC] = data[_WC] * b.data[_ZC] + data[_XC] * b.data[_YC] - data[_YC] * b.data[_XC] + data[_ZC] * b.data[_WC];
	return res;
}

/* function is broken */
qtnn_c qtnn_c::mult_vec3(const vec3_c &b) {
	qtnn_c res;
	res.data[_WC] = -data[_WC] * b.data[_XC] - data[_YC] * b.data[_YC] - data[_ZC] * b.data[_ZC];
    res.data[_XC] =  data[_WC] * b.data[_XC] + data[_YC] * b.data[_ZC] - data[_ZC] * b.data[_YC];
    res.data[_YC] =  data[_WC] * b.data[_YC] - data[_XC] * b.data[_ZC] + data[_ZC] * b.data[_XC];
    res.data[_ZC] =  data[_WC] * b.data[_ZC] + data[_XC] * b.data[_YC] - data[_YC] * b.data[_XC];
	return res;
}

vec3_c qtnn_c::transform_vec3(const vec3_c &b) {
	qtnn_c vq(b);
	qtnn_c tmp;// = this->mult_vec3(b);
	
	tmp = this->mult(vq);
	tmp = tmp.mult(this->invert());
	
	return tmp.to_vec3();
}

vec3_c qtnn_c::to_vec3() {
	return vec3_c(data[_XC], data[_YC], data[_ZC]);
}

/*------------------------------------------------------------------------*/
/*																		  */
/* N*N square matrix definition											  */
/*																		  */
/*------------------------------------------------------------------------

mtrxN_c::mtrxN_c(const int n) {
	int i;

	range = n;

	mtrx = new float[range*range];

	for (i = 0; i < range*range; i++)
		mtrx[i] = 0.0;

	for (i = 0; i < range; i++)
		mtrx[range*i] = 1.0;

}

mtrxN_c::~mtrxN_c() {
	delete mtrx;
}

int mtrxN_c::get_range() const {
	return range;
}

mtrxN_c mtrxN_c::mult(const mtrxN_c &b) {
	mtrxN_c ret(range);

	if (range != b.range) return ret;

	for (int j = 0; j < range; j++)
		for (int i = 0; i < range; i++) {
			ret.mtrx[INDXn(i,j,range)] = 0.0;
			for (int s = 0; s < range; s++)
				ret.mtrx[INDXn(i,j,range)] = ret.mtrx[INDXn(i,j,range)] + mtrx[INDXn(s,j,range)]*b.mtrx[INDXn(i,s,range)];
		}
}


void mtrxN_c::tranpose_self() {
	int i, j;
	float tmp;

	for (i = 0; i < range; i++) 
		for 
}*/
