#pragma once

const int _XC = 0;
const int _YC = 1;
const int _ZC = 2;
const int _WC = 3;

/*	multidimensional array mapping, array[i][j]
	row-wise (C, C++):
	(0	1)
	(2	3)

	column-wise (Fortran, Matlab):
	(0	2)
	(1	3)
*/
inline int id_rw(int i, int j, int n);	//row-wise
inline int id_cw(int i, int j, int n);	//column-wise

float lerp(float a, float b, float t);
void swap_float(float &a, float &b);

class vec3_c {
	public:
		float data[3];
		
		vec3_c(): 
			data{0.0, 0.0, 0.0} {};
		
		vec3_c(const float x, const float y, const float z): 
			data{x, y, z} {};
		
		vec3_c(const vec3_c &v): 
			data{v.data[_XC], v.data[_XC], v.data[_XC]} {};

		~vec3_c();
		
		void	set(const int index, const float x);
		void	set(const vec3_c &vec);
		float	get(const int index) const;
		
		vec3_c&	operator=(const vec3_c &v);
		float	operator[](const int index) const;
		float&  operator[](const int index);
		vec3_c	operator*(const float scalar);	//scale
		float	operator*(const vec3_c &v); 	//dot
		vec3_c	operator^(const vec3_c &v); 	//cross
		vec3_c	operator+(const vec3_c &v);
		vec3_c	operator-(const vec3_c &v);

		float	dot(const vec3_c &b) const;
		vec3_c	cross(const vec3_c &b) const;
		float	lenght() const;
		vec3_c 	normalize() const;
		void	normalize_self();
		vec3_c	invert() const;
		void	invert_self();
		
		vec3_c	scale(const float &t);
		vec3_c	add(const vec3_c &b);
		vec3_c	sub(const vec3_c &b);
};

class mtrx3_c {
	public:
		float data[9];

		mtrx3_c(): 
			data{1.0, 0.0, 0.0,
				 0.0, 1.0, 0.0,
				 0.0, 0.0, 1.0} {};
		
		mtrx3_c(const float yaw, const float pitch, const float roll) {
			from_euler(yaw, pitch, roll);
		};

		mtrx3_c(const vec3_c &ax, const float phi) {
			from_axis_angl(ax, phi);
		};

		mtrx3_c(const mtrx3_c &m):
			data{m.data[0], m.data[1], m.data[2],
				 m.data[3], m.data[4], m.data[5],
				 m.data[6], m.data[7], m.data[8]} {};
		
		~mtrx3_c();

		void	from_euler(const float yaw, const float pitch, const float roll);
		void	from_axis_angl(const vec3_c &ax, const float phi);
		void	set_idtt(); 
		
		float	det();
		mtrx3_c	inverse();
		void	tranpose();
		mtrx3_c	mult(const mtrx3_c &b);
		vec3_c	mult_vec3(const vec3_c &b);
};

class mtrx4_c {
	public:
		float mtrx[16];

		mtrx4_c(): mtrx{1.0, 0.0, 0.0, 0.0,
						0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 1.0, 0.0,
						0.0, 0.0, 0.0, 1.0} {};
		~mtrx4_c();

		void	from_euler(const float yaw, const float pitch, const float roll);
		void	from_axis_angl(const vec3_c &ax, const float phi);
		
		void	set_idtt(); 
		float	det();
		mtrx3_c	inverse();
		void	tranpose();
		mtrx3_c	mult(const mtrx3_c &b);
		vec3_c	mult_vec3(const vec3_c &b);
};

class qtnn_c {
	public:
		float data[4];

		qtnn_c(): 
			data{0.0, 0.0, 0.0, 1.0} {};
		
		qtnn_c(const float x, const float y, const float z, const float w):
			data{x, y, z, w} {};
		
		qtnn_c(const vec3_c &vq):
			data{vq.data[_XC], vq.data[_YC], vq.data[_ZC], 0.0} {};
		
		qtnn_c(const qtnn_c &q):
			data{q.data[_XC], q.data[_YC], q.data[_ZC], q.data[_WC]} {};

		qtnn_c(const mtrx3_c &b) {
			from_mtrx3(b);
		};

		qtnn_c(const float yaw, const float pitch, const float roll) {
			from_euler(yaw, pitch, roll);
		};

		qtnn_c(const vec3_c &ax, const float phi) {
			from_axis_angl(ax, phi);
		};

		~qtnn_c();
		
		void	set(const int index, const float x); /*empty*/
		void 	set(const qtnn_c &q);
		float	get(const int index) const; 		/*empty*/
		
		qtnn_c&	operator=(const qtnn_c &q);			
		float	operator[](const int index) const;	/*empty*/
		qtnn_c	operator*(const float scalar);		/*empty*/
		float	operator*(const qtnn_c &v); 		/*empty*/
		qtnn_c	operator^(const qtnn_c &v); 		/*empty*/
		qtnn_c	operator+(const qtnn_c &v);			/*empty*/
		qtnn_c	operator-(const qtnn_c &v);			/*empty*/

		float	lenght() const;
		qtnn_c	normalize() const;		/*empty*/
		void	normalize_self();
		qtnn_c	invert() const;
		void	invert_self();			/*empty*/

		qtnn_c	from_mtrx3(const mtrx3_c &b);	/*empty*/
		void	from_euler(const float yaw, const float pitch, const float roll);
		void	from_axis_angl(const vec3_c &ax, const float phi);
		qtnn_c	scale(const float t);
		qtnn_c	add(const qtnn_c &q);
		qtnn_c	mult(const qtnn_c &b);
		qtnn_c	mult_vec3(const vec3_c &b);
		vec3_c	transform_vec3(const vec3_c &b);
		
		vec3_c 	to_vec3();
};

class plane_c {
	public:
		float plane[4];
		
		plane_c();
		~plane_c();
};

class line_c {
	public:
		float line[3];
		
		line_c();
		~line_c();
};

/*
class mtrxN_c
{
	private:
		float *mtrx;
		int range;

	public:
		mtrxN_c(const int n);
		~mtrxN_c();
		
		int			get_range() const;
		mtrxN_c		mult(const mtrxN_c &b);
		void		tranpose_self();
};
*/
