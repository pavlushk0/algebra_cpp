
#include <iostream>
#include <tuple>
#include <cmath>
#include "algebra.h"

using namespace std;

mtrx3_t mtrx_set_euler(float yaw, float pitch, float roll) {
	mtrx3_t rt;
	float cosy, siny, cosp, sinp, cosr, sinr;
	
	cosy = cos(yaw);
	siny = sin(yaw);
	cosp = cos(pitch);
	sinp = sin(pitch);
	cosr = cos(roll);
	sinr = sin(roll);

	rt[0] = cosy*cosr - siny*cosp*sinr;
	rt[1] = -cosy*sinr - siny*cosp*cosr;
	rt[2] = siny * sinp;

	rt[3] = siny*cosr + cosy*cosp*sinr;
	rt[4] = -siny*sinr + cosy*cosp*cosr;
	rt[5] = -cosy * sinp;

	rt[6] = sinp * sinr;
	rt[7] = sinp * cosr;
	rt[8] = cosp;

	return rt;
}

mtrx3_t mtrx_set_axisangl(vec3_t &ax, float phi) {
	mtrx3_t rt;
	float cosphi, sinphi, vxvy, vxvz, vyvz, vx, vy, vz;

	cosphi = cos(phi);
	sinphi = sin(phi);
	vxvy = ax[_XC] * ax[_YC];
	vxvz = ax[_XC] * ax[_ZC];
	vyvz = ax[_YC] * ax[_ZC];
	vx = ax[_XC];
	vy = ax[_YC];
	vz = ax[_ZC];

	rt[0] = cosphi + (1.0-cosphi)*vx*vx;
	rt[1] = (1.0-cosphi)*vxvy - sinphi*vz;
	rt[2] = (1.0-cosphi)*vxvz + sinphi*vy;

	rt[3] = (1.0-cosphi)*vxvy + sinphi*vz;
	rt[4] = cosphi + (1.0-cosphi)*vy*vy;
	rt[5] = (1.0-cosphi)*vyvz - sinphi*vz;

	rt[6] = (1.0-cosphi)*vxvz - sinphi*vy;
	rt[7] = (1.0-cosphi)*vyvz + sinphi*vx;
	rt[8] = cosphi + (1.0-cosphi)*vz*vz;

	return rt;
}

void mtrx_show(const mtrx3_t &m) {
	printf("%5.2f %5.2f %5.2f\n", m[0], m[1], m[2]);
	printf("%5.2f %5.2f %5.2f\n", m[3], m[4], m[5]);
	printf("%5.2f %5.2f %5.2f\n", m[6], m[7], m[8]);
}

void mtrx_show(const mtrx4_t &m) {
    printf("%5.2f %5.2f %5.2f %5.2f\n", m[0],  m[1],  m[2],  m[3]);
	printf("%5.2f %5.2f %5.2f %5.2f\n", m[4],  m[5],  m[6],  m[7]);
	printf("%5.2f %5.2f %5.2f %5.2f\n", m[8],  m[9],  m[10], m[11]);
    printf("%5.2f %5.2f %5.2f %5.2f\n", m[11], m[12], m[13], m[15]);
}

float mtrx_det(const mtrx3_t &m) {
	return m[0]*m[4]*m[8] +
		   m[6]*m[1]*m[5] +
		   m[2]*m[3]*m[7] -
		   m[0]*m[7]*m[5] -
		   m[8]*m[3]*m[1];
}

template <typename mtrxT_t, int mrange>
mtrxT_t mtrx_mult(const mtrxT_t &a, const mtrxT_t &b) {
	mtrxT_t rt;
	int32_t i, j, k;
	float tmp;
	
	for (i = 0; i < mrange; i++) {
		for (j = 0; j < mrange; j++) {
			tmp = 0;
			for (k = 0; k < mrange; k++) {
				tmp = tmp + a[id_rw(k, j, mrange)]*b[id_rw(i, k, mrange)];
			}
			rt[id_rw(i, j, mrange)] = tmp;	
		}
	}

	return rt;
}
template mtrx3_t mtrx_mult<mtrx3_t, 3>(const mtrx3_t &a, const mtrx3_t &b);
template mtrx4_t mtrx_mult<mtrx4_t, 4>(const mtrx4_t &a, const mtrx4_t &b);

vec3_t mtrx_mult_vec3(const mtrx3_t &m, const vec3_t &v) {
	vec3_t rt;

	rt[_XC] = m[0]*v[_XC] + m[1]*v[_YC] + m[2]*v[_ZC];
	rt[_YC] = m[3]*v[_XC] + m[4]*v[_YC] + m[5]*v[_ZC];
	rt[_ZC] = m[6]*v[_XC] + m[7]*v[_YC] + m[8]*v[_ZC];

	return rt;
}

/*
	Где-то здесь ошибка, долго искал
	ничего не вышло и взял код из сети
*/
/*
func mtrx3_lu(m mtrx3_t) (l, u mtrx3_t) {
	var (
		i, j, k int32
		lm, um  mtrx3_t
		sum     float32
	)

	for j = 0; j < 3; j++ {
		um[id_rw(0, j, 3)] = m[id_rw(0, j, 3)]
	}

	for j = 0; j < 3; j++ {
		lm[id_rw(j, 0, 3)] = m[id_rw(j, 0, 3)] / um[id_rw(0, 0, 3)]
	}

	for i = 1; i < 3; i++ {
		for j = i; j < 3; j++ {
			sum = 0.0
			for k = 0; k < i; k++ {
				sum = sum + lm[id_rw(i, k, 3)]*um[id_rw(k, j, 3)]
			}
			um[id_rw(i, j, 3)] = m[id_rw(i, j, 3)] - sum
		}
	}

	for i = 1; i < 3; i++ {
		for j = i; j < 3; j++ {
			if i > j {
				lm[id_rw(j, i, 3)] = 0.0
			} else {
				sum = 0.0
				for k = 0; k < i; k++ {
					sum = sum + lm[id_rw(j, k, 3)]*um[id_rw(k, i, 3)]
				}
				lm[id_rw(j, i, 3)] = (1.0 / um[id_rw(i, i, 3)]) * (m[id_rw(j, i, 3)] - sum)
			}
		}
	}

	return lm, um
}
*/

/*
	Нижнетреугольная (L, lm) матрица имеет единицы по диагонали
*/

template <typename mtrxT_t, int mrange>
tuple<mtrxT_t, mtrxT_t> mtrx_lu(const mtrxT_t &m) {
	mtrxT_t lm, um;
	int32_t	i, j, k; 
	float sum;
	
	for (i = 0; i < mrange; i++) {
		for (k = i; k < mrange; k++) {
			sum = 0;
			for (j = 0; j < i; j++) {
				sum += (lm[id_rw(i, j, mrange)] * um[id_rw(j, k, mrange)]);
			}
			um[id_rw(i, k, mrange)] = m[id_rw(i, k, mrange)] - sum;
		}

		for (k = i; k < mrange; k++) {
			if (i == k) {
				lm[id_rw(i, i, mrange)] = 1.0;
			} else {
				sum = 0;
				for (j = 0; j < i; j++) {
					sum += lm[id_rw(k, j, mrange)] * um[id_rw(j, i, mrange)];
				}
				lm[id_rw(k, i, mrange)] = (m[id_rw(k, i, mrange)] - sum) / um[id_rw(i, i, mrange)];
			}
		}
	}

	return {lm, um};
}
template tuple<mtrx3_t, mtrx3_t> mtrx_lu<mtrx3_t, 3>(const mtrx3_t &m);
template tuple<mtrx4_t, mtrx4_t> mtrx_lu<mtrx4_t, 4>(const mtrx4_t &m);


template <typename mtrxT_t, typename vecT_t, int mrange>
tuple<mtrxT_t, vecT_t> mtrx_ldlt(const mtrxT_t &m) {
	mtrxT_t lm;
	vecT_t dv;
	int32_t	i, j, k; 
	float sum;
	
	for (i = 0; i < mrange; i++) {
		for (j = i; j < mrange; j++) {
			sum = m[id_rw(j, i, mrange)];
			for (k = 0; k < i; k++) {
				sum = sum - lm[id_rw(i, k, mrange)]*dv[k]*lm[id_rw(j, k, mrange)];
				if (i == j) {
					if (sum <= 0) {
						cout << "mtrx_ldlt(): matrix is not positive definite \n" ;
						return {mtrxT_t(), vecT_t()};
					}
					dv[i] = sum;
					lm[id_rw(i, i, mrange)] = 1.0;
				} else {
					lm[id_rw(j, i, mrange)] = sum / dv[i];
				}
			}
		}
	}

	return {lm, dv};
}
template tuple<mtrx3_t, vec3_t> mtrx_ldlt<mtrx3_t, vec3_t, 3>(const mtrx3_t &m);
template tuple<mtrx4_t, vec4_t> mtrx_ldlt<mtrx4_t, vec4_t, 4>(const mtrx4_t &m);

template <typename mtrxT_t, int mrange>
mtrxT_t	mtrx_get_transpose(const mtrxT_t &m) {
	mtrxT_t rt;
	int32_t i, j;
	float tmp;
	
	rt = m;

	for (i = 0; i < mrange; i++) {
		for (j = 0; j < i; j++) {
			tmp = rt[id_rw(i, i, mrange)];
			rt[id_rw(i, j, mrange)] = rt[id_rw(j, i, mrange)];
			rt[id_rw(j, i, mrange)] = tmp;
		}
	}

	return rt;
}
template mtrx3_t mtrx_get_transpose<mtrx3_t, 3>(const mtrx3_t &m);
template mtrx4_t mtrx_get_transpose<mtrx4_t, 4>(const mtrx4_t &m);

mtrx3_t	mtrx_get_inv(const mtrx3_t &m) {
	mtrx3_t inverse, rt;
	float det, invDet;

	inverse[id_rw(0, 0, 3)] = m[id_rw(1, 1, 3)] * m[id_rw(2, 2, 3)] - m[id_rw(1, 2, 3)] * m[id_rw(2, 1, 3)];
	inverse[id_rw(1, 0, 3)] = m[id_rw(1, 2, 3)] * m[id_rw(2, 0, 3)] - m[id_rw(1, 0, 3)] * m[id_rw(2, 2, 3)];
	inverse[id_rw(2, 0, 3)] = m[id_rw(1, 0, 3)] * m[id_rw(2, 1, 3)] - m[id_rw(1, 1, 3)] * m[id_rw(2, 0, 3)];

	det = m[id_rw(0, 0, 3)] * inverse[id_rw(0, 0, 3)] + m[id_rw(0, 1, 3)] * inverse[id_rw(1, 0, 3)] + 
		  m[id_rw(0, 2, 3)] * inverse[id_rw(2, 0, 3)];

	if (fabs(det) < f_eps) {
		return mtrx3_t();
	}

	invDet = 1.0f / det;

	inverse[id_rw(0, 1, 3)] = m[id_rw(0, 2, 3)] * m[id_rw(2, 1, 3)] - m[id_rw(0, 1, 3)] * m[id_rw(2, 2, 3)];
	inverse[id_rw(0, 2, 3)] = m[id_rw(0, 1, 3)] * m[id_rw(1, 2, 3)] - m[id_rw(0, 2, 3)] * m[id_rw(1, 1, 3)];
	inverse[id_rw(1, 1, 3)] = m[id_rw(0, 0, 3)] * m[id_rw(2, 2, 3)] - m[id_rw(0, 2, 3)] * m[id_rw(2, 0, 3)];
	inverse[id_rw(1, 2, 3)] = m[id_rw(0, 2, 3)] * m[id_rw(1, 0, 3)] - m[id_rw(0, 0, 3)] * m[id_rw(1, 2, 3)];
	inverse[id_rw(2, 1, 3)] = m[id_rw(0, 1, 3)] * m[id_rw(2, 0, 3)] - m[id_rw(0, 0, 3)] * m[id_rw(2, 1, 3)];
	inverse[id_rw(2, 2, 3)] = m[id_rw(0, 0, 3)] * m[id_rw(1, 1, 3)] - m[id_rw(0, 1, 3)] * m[id_rw(1, 0, 3)];

	rt[id_rw(0, 0, 3)] = inverse[id_rw(0, 0, 3)] * invDet;
	rt[id_rw(0, 1, 3)] = inverse[id_rw(0, 1, 3)] * invDet;
	rt[id_rw(0, 2, 3)] = inverse[id_rw(0, 2, 3)] * invDet;

	rt[id_rw(1, 0, 3)] = inverse[id_rw(1, 0, 3)] * invDet;
	rt[id_rw(1, 1, 3)] = inverse[id_rw(1, 1, 3)] * invDet;
	rt[id_rw(1, 2, 3)] = inverse[id_rw(1, 2, 3)] * invDet;

	rt[id_rw(2, 0, 3)] = inverse[id_rw(2, 0, 3)] * invDet;
	rt[id_rw(2, 1, 3)] = inverse[id_rw(2, 1, 3)] * invDet;
	rt[id_rw(2, 2, 3)] = inverse[id_rw(2, 2, 3)] * invDet;

	return rt;
}

mtrx4_t mtrx_get_inv(const mtrx4_t m) {
	mtrx4_t rt;
 	float det;
    int i;

    rt[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    rt[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    rt[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    rt[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    rt[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    rt[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    rt[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    rt[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    rt[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    rt[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    rt[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    rt[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    rt[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    rt[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    rt[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    rt[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * rt[0] + m[1] * rt[4] + m[2] * rt[8] + m[3] * rt[12];

    if (det == 0) {
		cout << "mtrx4_get_inv(): determinant of inv matrix is zero\n";
		return mtrx4_t();
	}
        

    det = 1.0f / det;

    for (i = 0; i < 16; i++) {
        rt[i] = rt[i] * det;
	}

	return rt;
}

mtrx4_t mtrx_get_inv_gauss(const mtrx4_t &m) {
	return mtrx4_t();	
}