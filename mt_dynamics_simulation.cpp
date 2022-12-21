#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <fstream>
#include <random>

using namespace std;


const double dc = pow(10, -10);//coordinate difference step
const double dp = pow(10, -2);//angle difference step
const double dt = 2.0 * pow(10.0, -12.0);//sec -- time difference step

//Model parameters
const double mt_rad = 25 / 2 * pow(10, -9); // nm -- microtubule radius
const double k_B = 1.38 * pow(10.0, -23.0);//J/K -- Boltzmann constant
const double T = 300;//K -- temperature
const double r_0 = 0.12 * pow(10.0, -9.0); // potential well width m(nm)
const double d = 0.25 * pow(10.0, -9.0); //tubulin-tubulin bond shape parameter, m(nm)
const double A_la = 16.9 * k_B * T;
const double b_la = 9.1 * k_B * T; //potential well depth
const double A_lo = 17.6 * k_B * T;
const double b_lo = 15.5 * k_B * T;
const double k = 517 * k_B * T * pow(10.0, 18.0); // inter-dimer bond strength 1/m^2(nm)

const double pi = acos(-1.);
const double h = 0.2;//Pa*s -- viscosity
const double R = 2.0 * pow(10.0, -9.0);//m -- particle radius

//3j6e (GTP)
const vector<double> AE_intra = { 0,0,0,8.2*pi/180,-40*pi/180,7.0*pi/180 };// GTP-GTP equilibrium angle
const vector<double> AE_inter = { 0,0,0,9.1*pi/180,-73*pi/180,7.7*pi/180 }; //GTP-GDP equilibrium angle
const vector<double> B_intra = { 0,0,0,1100 * k_B * T,39 * k_B * T,1160 * k_B * T };
const vector<double> B_inter = { 0,0,0,350 * k_B * T,13 * k_B * T,410 * k_B * T };

//3j6f (GDP)
// const vector<double> AE_intra = { 0,0,0,9.4 * pi / 180,-48 * pi / 180,5.2 * pi / 180 };// GTP-GTP equilibrium angle
// const vector<double> AE_inter = { 0,0,0,5.2 * pi / 180,-141 * pi / 180,4.4 * pi / 180 }; // GTP-GDP equilibrium angle
// const vector<double> B_intra = { 0,0,0,930 * k_B * T,37 * k_B * T,990 * k_B * T };
// const vector<double> B_inter = { 0,0,0,1290 * k_B * T,12 * k_B * T,760 * k_B * T };

const vector<double> us = { 0,0,R };//upper interaction site coord
const vector<double> bs = { 0,0,-R };//lower interaction site coord
const vector<double> rs = { -R * cos(pi / 13),R * sin(pi / 13), -12 / 13 / 2 * pow(10,-9) };//right interaction site coord
const vector<double> ls = { R * cos(pi / 13),R * sin(pi / 13), 12 / 13 / 2 * pow(10,-9) };//left interaction site coord


const vector<double> coef = { dt / (6 * pi * R * h),dt / (6 * pi * R * h),dt / (6 * pi * R * h),dt / (8 * pi * pow(R,3.) * h),dt / (8 * pi * pow(R,3.) * h),dt / (8 * pi * pow(R,3.) * h) };

vector<double> operator+(vector<double> left, vector<double> right)//checked
{
	vector<double> result(left.size());
	for (unsigned int i = 0; i < left.size(); i++)
	{
		result[i] = left[i] + right[i];
	}
	return result;
}

vector<double> operator-(vector<double> left, vector<double> right)//checked
{
	vector<double> result(right.size());
	for (unsigned int i = 0; i < right.size(); i++)
	{
		result[i] = left[i] - right[i];
	}
	return result;
}

vector<double> operator*(vector<double> left, vector<double> right)
{
	vector<double> result(right.size());
	for (unsigned int i = 0; i < right.size();i++)
	{
		result[i] = left[i] * right[i];
	}
	return result;
}

vector<double> operator*(double left, vector<double> right)//checked
{
	vector<double> result(right.size());
	for (unsigned int i = 0; i < right.size();i++)
	{
		result[i] = left * right[i];
	}
	return result;
}

vector<double> sqrt(vector<double> x)
{
	vector<double> result(x.size());
	for (unsigned int i = 0; i < x.size(); i++)
	{
		result[i] = sqrt(x[i]);
	}
	return result;
}

vector<double> random_vector(unsigned int dim)//checked
{
	//randomizer
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(0, 1.0);
	vector<double> result(dim);
	for (unsigned int i = 0; i < dim; i++)
	{
		result[i] = distribution(generator);
	}
	return result;
}

vector<double> zero(unsigned int dim)//checked
{
	vector<double> result(dim);
	for (unsigned int i = 0; i < dim; i++)
	{
		result[i] = 0;
	}
	return result;
}

vector<vector<double>> split(vector<double> coord)
{
	vector<double> C(3);
	vector<double> A(3);

	for (int i = 0; i < 3; i++)
	{
		C[i] = coord[i];
		A[i] = coord[i + 3];
	}

	return { C,A };
}

vector<double> pr(unsigned int axis, vector<double> v)
{
	switch (axis)
	{
	case 0: return { 0,v[1],v[2] };
			break;
	case 1: return { v[0],0,v[2] };
			break;
	case 2: return { v[0],v[1],0 };
	}
}

double dot(vector<double> v1, vector<double> v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double abs(vector<double> v)
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

vector<double> cross(vector<double> v1, vector<double> v2)
{
	vector<double> result(3);
	result[0] = v1[1] * v2[2] - v1[2] * v2[1];
	result[1] = v1[2] * v2[0] - v1[0] * v2[2];
	result[2] = v1[0] * v2[1] - v1[1] * v2[0];

	return result;
}

vector<double> quat_rot(vector<double> v, vector<double> angles)
{
	double a = angles[0];
	double b = angles[1];
	double c = angles[2];//углы

	vector<vector<double>> RM(3);//rotation matrix

	for (int i = 0; i < 3; i++)
	{
		RM[i].resize(3);
	}
	//1-st row
	RM[0][0] = cos(b) * cos(c);
	RM[0][1] = cos(c) * sin(a) * sin(b) - cos(a) * sin(c);
	RM[0][2] = cos(a) * cos(c) * sin(b) + sin(a) * sin(c);

	//2-d row
	RM[1][0] = cos(b) * sin(c);
	RM[1][1] = cos(a) * cos(c) + sin(a) * sin(b) * sin(c);
	RM[1][2] = -cos(c) * sin(a) + cos(a) * sin(b) * sin(c);

	//3-d row
	RM[2][0] = -sin(b);
	RM[2][1] = cos(b) * sin(a);
	RM[2][2] = cos(a) * cos(b);

	vector<double> v_res(3);

	v_res[0] = RM[0][0] * v[0] + RM[0][1] * v[1] + RM[0][2] * v[2];//row*column;
	v_res[1] = RM[1][0] * v[0] + RM[1][1] * v[1] + RM[1][2] * v[2];
	v_res[2] = RM[2][0] * v[0] + RM[2][1] * v[1] + RM[2][2] * v[2];

	if (abs(v_res[0] * pow(10., 15.)) < 1.) { v_res[0] = 0; }
	if (abs(v_res[1] * pow(10., 15.)) < 1.) { v_res[1] = 0; }
	if (abs(v_res[2] * pow(10., 15.)) < 1.) { v_res[2] = 0; }

	return v_res;
}

vector<double> axis_rot(vector<double> v, vector<double> axis, double angle)
{
	vector<double> v_res(3);

	double w = cos(angle / 2);
	double x = sin(angle / 2) * axis[0];
	double y = sin(angle / 2) * axis[1];
	double z = sin(angle / 2) * axis[2];

	vector<vector<double>> RM(3);
	for (int i = 0; i < 3; i++)
	{
		RM[i].resize(3);
	}

	RM[0][0] = 1 - 2 * y * y - 2 * z * z;
	RM[0][1] = 2 * x * y - 2 * z * w;
	RM[0][2] = 2 * x * z + 2 * y * w;

	RM[1][0] = 2 * x * y + 2 * z * w;
	RM[1][1] = 1 - 2 * x * x - 2 * z * z;
	RM[1][2] = 2 * y * z - 2 * x * w;

	RM[2][0] = 2 * x * z - 2 * y * w;
	RM[2][1] = 2 * y * z + 2 * x * w;
	RM[2][2] = 1 - 2 * x * x - 2 * y * y;

	v_res[0] = RM[0][0] * v[0] + RM[0][1] * v[1] + RM[0][2] * v[2];//row*column;
	v_res[1] = RM[1][0] * v[0] + RM[1][1] * v[1] + RM[1][2] * v[2];
	v_res[2] = RM[2][0] * v[0] + RM[2][1] * v[1] + RM[2][2] * v[2];

	if (abs(v_res[0] * pow(10., 15.)) < 1.) { v_res[0] = 0; }
	if (abs(v_res[1] * pow(10., 15.)) < 1.) { v_res[1] = 0; }
	if (abs(v_res[2] * pow(10., 15.)) < 1.) { v_res[2] = 0; }

	return v_res;
}

vector<double> UpIS(vector<double> coordinates)
{
	vector<double> C(3);
	for (int i = 0; i < 3; i++)
	{
		C[i] = coordinates[i];
	}
	vector<double> A(3);
	for (int i = 3; i < 6; i++)
	{
		A[i - 3] = coordinates[i];
	}

	return C + quat_rot(us, A);
}

vector<double> BoIS(vector<double> coordinates)
{
	vector<double> C(3);
	for (int i = 0; i < 3; i++)
	{
		C[i] = coordinates[i];
	}

	vector<double> A(3);
	for (int i = 3; i < 6; i++)
	{
		A[i - 3] = coordinates[i];
	}


	return C + quat_rot(bs, A);
}

vector<double> LeIS(vector<double> coordinates, double psi)
{
	vector<double> C(3);
	for (int i = 0; i < 3; i++)
	{
		C[i] = coordinates[i];
	}
	vector<double> A(3);
	for (int i = 3; i < 6; i++)
	{
		A[i - 3] = coordinates[i];
	}

	vector<double> is(3);
	is[0] = -ls[0]*sin(pi + psi) + ls[1]*cos(pi + psi);
	is[1] = -ls[0]*cos(pi + psi) - ls[1]*sin(pi + psi);
	is[2] = ls[2];

	return C + quat_rot(is, A);
}

vector<double> RiIS(vector<double> coordinates, double psi)
{
	vector<double> C(3);
	for (int i = 0; i < 3; i++)
	{
		C[i] = coordinates[i];
	}
	vector<double> A(3);
	for (int i = 3; i < 6; i++)
	{
		A[i - 3] = coordinates[i];
	}

	vector<double> is(3);
	is[0] = -rs[0]*sin(pi + psi) + rs[1]*cos(pi + psi);
	is[1] = -rs[0]*cos(pi + psi) - rs[1]*sin(pi + psi);
	is[2] = rs[2];

	return C + quat_rot(is, A);
}

vector<double> dh(unsigned int c)
{
	vector<double> h(6);
	h = zero(6);
	if (c < 3) { h[c] = dc; }
	else { h[c] = dp; }
	return h;
}

vector<double> angles(vector<double> coord1, vector<double> coord2, double psi) //pf - filament number
{
	vector<vector<double>> v1 = split(coord1);
	vector<vector<double>> v2 = split(coord2);

	vector<double> ex = { - R * sin(psi), R * cos (psi), 0 };
	vector<double> ey = { - R * cos(psi), - R * sin (psi), 0 };
	vector<double> ez = { 0, 0, R };

	vector<double> eX = { - R * sin(psi), R * cos (psi), 0 };
	vector<double> eY = { - R * cos(psi), - R * sin (psi), 0 };
	vector<double> eZ = { 0, 0, R };

	ex = (1 / R) * quat_rot(ex, v1[1]);
	ey = (1 / R) * quat_rot(ey, v1[1]);
	ez = (1 / R) * quat_rot(ez, v1[1]);

	eX = (1 / R) * quat_rot(eX, v2[1]);
	eY = (1 / R) * quat_rot(eY, v2[1]);
	eZ = (1 / R) * quat_rot(eZ, v2[1]);

	vector<double> P = cross(ez, eZ);
	double n = abs(P);
	P = (1 / n) * P;

	double thetta;

	if (abs(dot(ez, eZ) - 1) * pow(10, 15) < 1)
	{
		thetta = 0;
	}
	else
	{
		thetta = acos(dot(ez, eZ));
	}

	if (abs(dot(ez, eZ)) > 1)
	{
		thetta = 0;
		//cout<<"ERR1:"<<dot(ez,eZ)<<endl;
	}

	double Z_x = dot(eZ, ex);
	double Z_y = dot(eZ, ey);
	double Z_z = dot(eZ, ez);

	//cout<<"Z:"<<Z_x<< ' '<<Z_y<<' '<<Z_z<<endl;

	double phi;
	if (Z_x > 0)
	{
		phi = acos(-Z_y / sqrt(1 - Z_z * Z_z));
	}
	else
	{
		phi = -acos(-Z_y / sqrt(1 - Z_z * Z_z));
	}

	if (abs(-Z_y / sqrt(1 - Z_z * Z_z)) > 1)
	{
		phi = 0;
		//cout<<"ERR2"<<endl;
	}

	vector<double> ex_(3);
	vector<double> ey_(3);
	vector<double> ez_(3);

	if (thetta == 0)
	{
		ex_ = ex;
		ey_ = ey;
		ez_ = ez;
	}
	else
	{
		ex_ = axis_rot(ex, P, thetta);
		ey_ = axis_rot(ey, P, thetta);
		ez_ = axis_rot(ez, P, thetta);
	}

	double delta;
	if (abs(dot(ex_, eX) - 1) * pow(10, 15) < 1)
	{
		delta = 0;
	}
	else
	{
		delta = acos(dot(ex_, eX));
	}

	if (abs(dot(ex_, eX)) > 1)
	{
		delta = 0;
		//cout<<"ERR3"<<endl;
	}

	if (dot(ex_, cross(eX, eZ)) < 0) { delta = -delta; }

	vector<double> result(3);

	result = { thetta, phi, delta };

	//cout<<result[0]<<' '<<result[1]<<' '<<result[2]<<endl;

	return result;
}

double AI(vector<double> coord1, vector<double> coord2, unsigned int pf)
{
	vector<double> A(3);
	A = angles(coord1, coord2, pf);
	return 0.5* B_intra[3] * pow((A[0] - AE_intra[3]), 2) + 0.5 * B_intra[4] * pow((A[1] - AE_intra[4]), 2) + 0.5 * B_intra[5] * pow((A[2] - AE_intra[5]), 2);
}

double AE(vector<double> coord1, vector<double> coord2, double psi)
{
	vector<double> A(3);
	A = angles(coord1, coord2, psi);
	return 0.5* B_inter[3] * pow((A[0] - AE_inter[3]), 2) + 0.5 * B_inter[4] * pow((A[1] - AE_inter[4]), 2) + 0.5 * B_inter[5] * pow((A[2] - AE_inter[5]), 2);
}

double LI(double dr)
{
	return 0.5* k* pow(dr, 2);
}

double LE(double dr)
{
	return A_lo * pow(dr / r_0, 2)* exp(-dr / r_0) - b_lo * exp(-pow(dr, 2) / d / r_0);
}

double LA(double dr, unsigned int b)//b - наличие связи (1), либо её отсутствие (0)
{
    if (b == 1)
	{return A_la * pow(dr / r_0, 2)* exp(-dr / r_0) - b_la * exp(-pow(dr, 2) / d / r_0);}
    else
    {return 0;}
}

double U0(vector<double> coord1, vector<double> coord2, double psi)
{
	double dr21 = sqrt(pow((BoIS(coord2)[0] - UpIS(coord1)[0]), 2) + pow((BoIS(coord2)[1] - UpIS(coord1)[1]), 2) + pow((BoIS(coord2)[2] - UpIS(coord1)[2]), 2));

	return LI(dr21) + AI(coord1,coord2, psi);
}

double U1(vector<double> coord1, vector<double> coord2, vector<double> coord3, double psi)
{
	double dr21 = sqrt(pow((BoIS(coord2)[0] - UpIS(coord1)[0]), 2) + pow((BoIS(coord2)[1] - UpIS(coord1)[1]), 2) + pow((BoIS(coord2)[2] - UpIS(coord1)[2]), 2));
	double dr32 = sqrt(pow((BoIS(coord3)[0] - UpIS(coord2)[0]), 2) + pow((BoIS(coord3)[1] - UpIS(coord2)[1]), 2) + pow((BoIS(coord3)[2] - UpIS(coord2)[2]), 2));


	return LE(dr32) + LI(dr21) + AE(coord2, coord3, psi) + AI(coord1,coord2, psi);
}

double U2(vector<double> coord1, vector<double> coord2, vector<double> coord3, double psi)
{
	double dr21 = sqrt(pow((BoIS(coord2)[0] - UpIS(coord1)[0]), 2) + pow((BoIS(coord2)[1] - UpIS(coord1)[1]), 2) + pow((BoIS(coord2)[2] - UpIS(coord1)[2]), 2));
	double dr32 = sqrt(pow((BoIS(coord3)[0] - UpIS(coord2)[0]), 2) + pow((BoIS(coord3)[1] - UpIS(coord2)[1]), 2) + pow((BoIS(coord3)[2] - UpIS(coord2)[2]), 2));

	return LI(dr32) + LE(dr21) + AI(coord2,coord3, psi) + AE(coord1,coord2, psi);
}

double U3(vector<double> coord2, vector<double> coord3, double psi)
{
	double dr32 = sqrt(pow((BoIS(coord3)[0] - UpIS(coord2)[0]), 2) + pow((BoIS(coord3)[1] - UpIS(coord2)[1]), 2) + pow((BoIS(coord3)[2] - UpIS(coord2)[2]), 2));

	return LI(dr32) + AI(coord2,coord3, psi);
}

double UL(vector<double> coord, vector<double> coord_l, vector<double> coord_r, double psi, unsigned int b)
{
	double psi_l = psi + 2 * pi / 13;
	double psi_r = psi - 2 * pi / 13;
	double dr52 = sqrt(pow((LeIS(coord, psi)[0] - RiIS(coord_l, psi_l)[0]), 2) + pow((LeIS(coord, psi)[1] - RiIS(coord_l, psi_l)[1]), 2) + pow((LeIS(coord, psi)[2] - RiIS(coord_l, psi_l)[2]), 2));
	double dr24 = sqrt(pow((LeIS(coord_r, psi_r)[0] - RiIS(coord, psi)[0]), 2) + pow((LeIS(coord_r, psi_r)[1] - RiIS(coord, psi)[1]), 2) + pow((LeIS(coord_r, psi_r)[2] - RiIS(coord, psi)[2]), 2));
	return LA(dr24, b) + LA(dr52, b);
}

vector<vector<double>> SFL(vector<vector<double>> coord, double psi)
{
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
	}

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		//monomer classification
		unsigned int counter;
		if (i == 0) { counter = 0; }
		if ((i != 0) && (i != (coord.size() - 1)) && (i % 2 != 0)) { counter = 1; }
		if ((i != 0) && (i != (coord.size() - 1)) && (i % 2 == 0)) { counter = 2; }
		if (i == (coord.size() - 1)) { counter = 3; }
		double delta;

		switch (counter)
		{
		case 0:
            delta = dc;
			result[i][0] = (U0(coord[i] + dh(0), coord[i + 1], psi) - U0(coord[i] - dh(0), coord[i + 1], psi)) / 2 / delta;
			result[i][1] = (U0(coord[i] + dh(1), coord[i + 1], psi) - U0(coord[i] - dh(1), coord[i + 1], psi)) / 2 / delta;
			result[i][2] = (U0(coord[i] + dh(2), coord[i + 1], psi) - U0(coord[i] - dh(2), coord[i + 1], psi)) / 2 / delta;
			delta = dp;
			result[i][3] = (U0(coord[i] + dh(3), coord[i + 1], psi) - U0(coord[i] - dh(3), coord[i + 1], psi)) / 2 / delta;
			result[i][4] = (U0(coord[i] + dh(4), coord[i + 1], psi) - U0(coord[i] - dh(4), coord[i + 1], psi)) / 2 / delta;
			result[i][5] = (U0(coord[i] + dh(5), coord[i + 1], psi) - U0(coord[i] - dh(5), coord[i + 1], psi)) / 2 / delta;
			result[i] = coef * result[i];
			break;

		case 1:
            delta = dc;
			result[i][0] = (U1(coord[i - 1], coord[i] + dh(0), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(0), coord[i + 1], psi)) / 2 / delta;
			result[i][1] = (U1(coord[i - 1], coord[i] + dh(1), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(1), coord[i + 1], psi)) / 2 / delta;
			result[i][2] = (U1(coord[i - 1], coord[i] + dh(2), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(2), coord[i + 1], psi)) / 2 / delta;
			delta = dp;
			result[i][3] = (U1(coord[i - 1], coord[i] + dh(3), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(3), coord[i + 1], psi)) / 2 / delta;
			result[i][4] = (U1(coord[i - 1], coord[i] + dh(4), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(4), coord[i + 1], psi)) / 2 / delta;
			result[i][5] = (U1(coord[i - 1], coord[i] + dh(5), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(5), coord[i + 1], psi)) / 2 / delta;
			result[i] = coef * result[i];
			break;

		case 2:
            delta = dc;
			result[i][0] = (U2(coord[i - 1], coord[i] + dh(0), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(0), coord[i + 1], psi)) / 2 / delta;
			result[i][1] = (U2(coord[i - 1], coord[i] + dh(1), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(1), coord[i + 1], psi)) / 2 / delta;
			result[i][2] = (U2(coord[i - 1], coord[i] + dh(2), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(2), coord[i + 1], psi)) / 2 / delta;
			delta = dp;
			result[i][3] = (U2(coord[i - 1], coord[i] + dh(3), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(3), coord[i + 1], psi)) / 2 / delta;
			result[i][4] = (U2(coord[i - 1], coord[i] + dh(4), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(4), coord[i + 1], psi)) / 2 / delta;
			result[i][5] = (U2(coord[i - 1], coord[i] + dh(5), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(5), coord[i + 1], psi)) / 2 / delta;
			result[i] = coef * result[i];
			break;

		case 3:
            delta = dc;
			result[i][0] = (U3(coord[i - 1], coord[i] + dh(0), psi) - U3(coord[i - 1], coord[i] - dh(0), psi)) / 2 / delta;
			result[i][1] = (U3(coord[i - 1], coord[i] + dh(1), psi) - U3(coord[i - 1], coord[i] - dh(1), psi)) / 2 / delta;
			result[i][2] = (U3(coord[i - 1], coord[i] + dh(2), psi) - U3(coord[i - 1], coord[i] - dh(2), psi)) / 2 / delta;
			delta = dp;
			result[i][3] = (U3(coord[i - 1], coord[i] + dh(3), psi) - U3(coord[i - 1], coord[i] - dh(3), psi)) / 2 / delta;
			result[i][4] = (U3(coord[i - 1], coord[i] + dh(4), psi) - U3(coord[i - 1], coord[i] - dh(4), psi)) / 2 / delta;
			result[i][5] = (U3(coord[i - 1], coord[i] + dh(5), psi) - U3(coord[i - 1], coord[i] - dh(5), psi)) / 2 / delta;
			result[i] = coef * result[i];
			break;
		}
	}
	
	return result;

}

vector<vector<double>> RFL(vector<vector<double>> coord)//Longitudinal random force
{
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
	}

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i] = sqrt(2 * k_B * T) * sqrt(coef) * random_vector(6);
	}

	return(result);
}

vector<vector<double>> TFL(vector<vector<double>> coord, double psi)//Longitudinal total force
{
	vector<vector<double>> sf(coord.size());
	vector<vector<double>> rf(coord.size());
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		sf[i].resize(6);
		rf[i].resize(6);
		result[i].resize(6);
	}

	sf = SFL(coord, psi);
	rf = RFL(coord);

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i] = zero(6) - sf[i];
	    //result[i] = rf[i] - sf[i];
	}

	return result;
}

vector<vector<double>> ev(vector<vector<double>> coord, double psi)//Equations
{
	vector<vector<double>> result(coord.size());
	vector<vector<double>> force(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
		force[i].resize(6);
	}

	force = TFL(coord, psi);

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i] = coord[i] + force[i];
	}

	return result;
}

vector<vector<double>> SF_long(vector<vector<double>> coord, double psi)
{
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
	}

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		//monomer classification
		unsigned int counter;
		if (i == 0) { counter = 0; }
		if ((i != 0) && (i != (coord.size() - 1)) && (i % 2 != 0)) { counter = 1; }
		if ((i != 0) && (i != (coord.size() - 1)) && (i % 2 == 0)) { counter = 2; }
		if (i == (coord.size() - 1)) { counter = 3; }
		double delta;

		switch (counter)
		{
		case 0:
            delta = dc;
			result[i][0] = (U0(coord[i] + dh(0), coord[i + 1], psi) - U0(coord[i] - dh(0), coord[i + 1], psi)) / 2 / delta;
			result[i][1] = (U0(coord[i] + dh(1), coord[i + 1], psi) - U0(coord[i] - dh(1), coord[i + 1], psi)) / 2 / delta;
			result[i][2] = (U0(coord[i] + dh(2), coord[i + 1], psi) - U0(coord[i] - dh(2), coord[i + 1], psi)) / 2 / delta;
			delta = dp;
			result[i][3] = (U0(coord[i] + dh(3), coord[i + 1], psi) - U0(coord[i] - dh(3), coord[i + 1], psi)) / 2 / delta;
			result[i][4] = (U0(coord[i] + dh(4), coord[i + 1], psi) - U0(coord[i] - dh(4), coord[i + 1], psi)) / 2 / delta;
			result[i][5] = (U0(coord[i] + dh(5), coord[i + 1], psi) - U0(coord[i] - dh(5), coord[i + 1], psi)) / 2 / delta;
			result[i] = coef * result[i];
			break;

		case 1:
            delta = dc;
			result[i][0] = (U1(coord[i - 1], coord[i] + dh(0), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(0), coord[i + 1], psi)) / 2 / delta;
			result[i][1] = (U1(coord[i - 1], coord[i] + dh(1), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(1), coord[i + 1], psi)) / 2 / delta;
			result[i][2] = (U1(coord[i - 1], coord[i] + dh(2), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(2), coord[i + 1], psi)) / 2 / delta;
			delta = dp;
			result[i][3] = (U1(coord[i - 1], coord[i] + dh(3), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(3), coord[i + 1], psi)) / 2 / delta;
			result[i][4] = (U1(coord[i - 1], coord[i] + dh(4), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(4), coord[i + 1], psi)) / 2 / delta;
			result[i][5] = (U1(coord[i - 1], coord[i] + dh(5), coord[i + 1], psi) - U1(coord[i - 1], coord[i] - dh(5), coord[i + 1], psi)) / 2 / delta;
			result[i] = coef * result[i];
			break;

		case 2:
            delta = dc;
			result[i][0] = (U2(coord[i - 1], coord[i] + dh(0), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(0), coord[i + 1], psi)) / 2 / delta;
			result[i][1] = (U2(coord[i - 1], coord[i] + dh(1), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(1), coord[i + 1], psi)) / 2 / delta;
			result[i][2] = (U2(coord[i - 1], coord[i] + dh(2), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(2), coord[i + 1], psi)) / 2 / delta;
			delta = dp;
			result[i][3] = (U2(coord[i - 1], coord[i] + dh(3), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(3), coord[i + 1], psi)) / 2 / delta;
			result[i][4] = (U2(coord[i - 1], coord[i] + dh(4), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(4), coord[i + 1], psi)) / 2 / delta;
			result[i][5] = (U2(coord[i - 1], coord[i] + dh(5), coord[i + 1], psi) - U2(coord[i - 1], coord[i] - dh(5), coord[i + 1], psi)) / 2 / delta;
			result[i] = coef * result[i];
			break;

		case 3:
            delta = dc;
			result[i][0] = (U3(coord[i - 1], coord[i] + dh(0), psi) - U3(coord[i - 1], coord[i] - dh(0), psi)) / 2 / delta;
			result[i][1] = (U3(coord[i - 1], coord[i] + dh(1), psi) - U3(coord[i - 1], coord[i] - dh(1), psi)) / 2 / delta;
			result[i][2] = (U3(coord[i - 1], coord[i] + dh(2), psi) - U3(coord[i - 1], coord[i] - dh(2), psi)) / 2 / delta;
			delta = dp;
			result[i][3] = (U3(coord[i - 1], coord[i] + dh(3), psi) - U3(coord[i - 1], coord[i] - dh(3), psi)) / 2 / delta;
			result[i][4] = (U3(coord[i - 1], coord[i] + dh(4), psi) - U3(coord[i - 1], coord[i] - dh(4), psi)) / 2 / delta;
			result[i][5] = (U3(coord[i - 1], coord[i] + dh(5), psi) - U3(coord[i - 1], coord[i] - dh(5), psi)) / 2 / delta;
			result[i] = coef * result[i];
			break;
		}
	}

	return result;

}

vector<vector<double>> SF_lat(vector<vector<double>> coord, vector<vector<double>> coord_l, vector<vector<double>> coord_r, double psi, unsigned int b)
{
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
	}
	
	for(unsigned int i = 0; i < coord.size(); i++)
	{
		//cout << "PF: " << i << endl;

		double delta;

		delta = dc;
		result[i][0] = (UL(coord[i]+dh(0), coord_l[i], coord_r[i], psi, b) - UL(coord[i]-dh(0), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][1] = (UL(coord[i]+dh(1), coord_l[i], coord_r[i], psi, b) - UL(coord[i]-dh(1), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][2] = (UL(coord[i]+dh(2), coord_l[i], coord_r[i], psi, b) - UL(coord[i]-dh(2), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		delta = dp;
		result[i][3] = (UL(coord[i]+dh(3), coord_l[i], coord_r[i], psi, b) - UL(coord[i]-dh(3), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][4] = (UL(coord[i]+dh(4), coord_l[i], coord_r[i], psi, b) - UL(coord[i]-dh(4), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][5] = (UL(coord[i]+dh(5), coord_l[i], coord_r[i], psi, b) - UL(coord[i]-dh(5), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i] = coef * result[i];
	}

	return result;
}


double US_0(vector<double> coord, vector<double> coord_l, vector<double> coord_r, double psi, unsigned int b)
{
	double psi_l = psi + 2 * pi / 13;
	//double psi_r = psi - 2 * pi / 13;
	double dr52 = sqrt(pow((LeIS(coord, psi)[0] - RiIS(coord_l, psi_l)[0]), 2) + pow((LeIS(coord, psi)[1] - RiIS(coord_l, psi_l)[1]), 2) + pow((LeIS(coord, psi)[2] - RiIS(coord_l, psi_l)[2]), 2));
	//double dr24 = sqrt(pow((LeIS(coord_r, psi_r)[0] - RiIS(coord, psi)[0]), 2) + pow((LeIS(coord_r, psi_r)[1] - RiIS(coord, psi)[1]), 2) + pow((LeIS(coord_r, psi_r)[2] - RiIS(coord, psi)[2]), 2));
	return /*LA(dr24, b) +*/ LA(dr52, b);
}

double US_12(vector<double> coord, vector<double> coord_l, vector<double> coord_r, double psi, unsigned int b)
{
	//double psi_l = psi + 2 * pi / 13;
	double psi_r = psi - 2 * pi / 13;
	//double dr52 = sqrt(pow((LeIS(coord, psi)[0] - RiIS(coord_l, psi_l)[0]), 2) + pow((LeIS(coord, psi)[1] - RiIS(coord_l, psi_l)[1]), 2) + pow((LeIS(coord, psi)[2] - RiIS(coord_l, psi_l)[2]), 2));
	double dr24 = sqrt(pow((LeIS(coord_r, psi_r)[0] - RiIS(coord, psi)[0]), 2) + pow((LeIS(coord_r, psi_r)[1] - RiIS(coord, psi)[1]), 2) + pow((LeIS(coord_r, psi_r)[2] - RiIS(coord, psi)[2]), 2));
	return LA(dr24, b) /*+ LA(dr52, b)*/;
}


vector<vector<double>> seam_0(vector<vector<double>> coord, vector<vector<double>> coord_l, vector<vector<double>> coord_r, double psi, unsigned int b)
{
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
	}
	
	for(unsigned int i = 0; i < coord.size(); i++)
	{
		//cout << "PF: " << i << endl;

		double delta;

		delta = dc;
		result[i][0] = (US_0(coord[i]+dh(0), coord_l[i], coord_r[i], psi, b) - US_0(coord[i]-dh(0), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][1] = (US_0(coord[i]+dh(1), coord_l[i], coord_r[i], psi, b) - US_0(coord[i]-dh(1), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][2] = (US_0(coord[i]+dh(2), coord_l[i], coord_r[i], psi, b) - US_0(coord[i]-dh(2), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		delta = dp;
		result[i][3] = (US_0(coord[i]+dh(3), coord_l[i], coord_r[i], psi, b) - US_0(coord[i]-dh(3), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][4] = (US_0(coord[i]+dh(4), coord_l[i], coord_r[i], psi, b) - US_0(coord[i]-dh(4), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][5] = (US_0(coord[i]+dh(5), coord_l[i], coord_r[i], psi, b) - US_0(coord[i]-dh(5), coord_l[i], coord_r[i], psi, b)) / 2 / delta;

		if (i == (coord.size() - 1))
	{

		delta = dc;
		result[i][0] = (UL(coord[i]+dh(0), coord_l[i], coord_r[i-3], psi, b) - UL(coord[i]-dh(0), coord_l[i], coord_r[i-3], psi, b)) / 2 / delta;
		result[i][1] = (UL(coord[i]+dh(1), coord_l[i], coord_r[i-3], psi, b) - UL(coord[i]-dh(1), coord_l[i], coord_r[i-3], psi, b)) / 2 / delta;
		result[i][2] = (UL(coord[i]+dh(2), coord_l[i], coord_r[i-3], psi, b) - UL(coord[i]-dh(2), coord_l[i], coord_r[i-3], psi, b)) / 2 / delta;
		delta = dp;
		result[i][3] = (UL(coord[i]+dh(3), coord_l[i], coord_r[i-3], psi, b) - UL(coord[i]-dh(3), coord_l[i], coord_r[i-3], psi, b)) / 2 / delta;
		result[i][4] = (UL(coord[i]+dh(4), coord_l[i], coord_r[i-3], psi, b) - UL(coord[i]-dh(4), coord_l[i], coord_r[i-3], psi, b)) / 2 / delta;
		result[i][5] = (UL(coord[i]+dh(5), coord_l[i], coord_r[i-3], psi, b) - UL(coord[i]-dh(5), coord_l[i], coord_r[i-3], psi, b)) / 2 / delta;
	}


		result[i] = coef * result[i];
	}

	return result;
}


vector<vector<double>> seam_12(vector<vector<double>> coord, vector<vector<double>> coord_l, vector<vector<double>> coord_r, double psi, unsigned int b)
{
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
	}
	
	for(unsigned int i = 0; i < coord.size(); i++)
	{
		//cout << "PF: " << i << endl;

		double delta;

		delta = dc;
		result[i][0] = (US_12(coord[i]+dh(0), coord_l[i], coord_r[i], psi, b) - US_12(coord[i]-dh(0), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][1] = (US_12(coord[i]+dh(1), coord_l[i], coord_r[i], psi, b) - US_12(coord[i]-dh(1), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][2] = (US_12(coord[i]+dh(2), coord_l[i], coord_r[i], psi, b) - US_12(coord[i]-dh(2), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		delta = dp;
		result[i][3] = (US_12(coord[i]+dh(3), coord_l[i], coord_r[i], psi, b) - US_12(coord[i]-dh(3), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][4] = (US_12(coord[i]+dh(4), coord_l[i], coord_r[i], psi, b) - US_12(coord[i]-dh(4), coord_l[i], coord_r[i], psi, b)) / 2 / delta;
		result[i][5] = (US_12(coord[i]+dh(5), coord_l[i], coord_r[i], psi, b) - US_12(coord[i]-dh(5), coord_l[i], coord_r[i], psi, b)) / 2 / delta;

		if (i == 0)
	{

		delta = dc;
		result[i][0] = (UL(coord[i]+dh(0), coord_l[i+3], coord_r[i], psi, b) - UL(coord[i]-dh(0), coord_l[i+3], coord_r[i], psi, b)) / 2 / delta;
		result[i][1] = (UL(coord[i]+dh(1), coord_l[i+3], coord_r[i], psi, b) - UL(coord[i]-dh(1), coord_l[i+3], coord_r[i], psi, b)) / 2 / delta;
		result[i][2] = (UL(coord[i]+dh(2), coord_l[i+3], coord_r[i], psi, b) - UL(coord[i]-dh(2), coord_l[i+3], coord_r[i], psi, b)) / 2 / delta;
		delta = dp;
		result[i][3] = (UL(coord[i]+dh(3), coord_l[i+3], coord_r[i], psi, b) - UL(coord[i]-dh(3), coord_l[i+3], coord_r[i], psi, b)) / 2 / delta;
		result[i][4] = (UL(coord[i]+dh(4), coord_l[i+3], coord_r[i], psi, b) - UL(coord[i]-dh(4), coord_l[i+3], coord_r[i], psi, b)) / 2 / delta;
		result[i][5] = (UL(coord[i]+dh(5), coord_l[i+3], coord_r[i], psi, b) - UL(coord[i]-dh(5), coord_l[i+3], coord_r[i], psi, b)) / 2 / delta;
	}


		result[i] = coef * result[i];
	}

	return result;
}



vector<vector<double>> RF(vector<vector<double>> coord)//random force
{
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
	}

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i] = sqrt(2 * k_B * T) * sqrt(coef) * random_vector(6);
	}

	return(result);
}

vector<vector<double>> TF(vector<vector<double>> coord, vector<vector<double>> coord_l, vector<vector<double>> coord_r, double psi, unsigned int b)//Longitudinal total force
{
	vector<vector<double>> sf(coord.size());
	vector<vector<double>> rf(coord.size());
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		sf[i].resize(6);
		rf[i].resize(6);
		result[i].resize(6);
	}

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		sf[i] = SF_long(coord, psi)[i] + SF_lat(coord, coord_l, coord_r, psi, b)[i];
	}

	rf = RF(coord);

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		//result[i] = zero(6) - sf[i];
	    result[i] = rf[i] - sf[i];
	}

	return result;
}

vector<vector<double>> TF_seam_0(vector<vector<double>> coord, vector<vector<double>> coord_l, vector<vector<double>> coord_r, double psi, unsigned int b)//Longitudinal total force
{
	vector<vector<double>> sf(coord.size());
	vector<vector<double>> rf(coord.size());
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		sf[i].resize(6);
		rf[i].resize(6);
		result[i].resize(6);
	}

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		sf[i] = SF_long(coord, psi)[i] + seam_0(coord, coord_l, coord_r, psi, b)[i];
	}

	rf = RF(coord);

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		//result[i] = zero(6) - sf[i];
	    result[i] = rf[i] - sf[i];
	}

	return result;
}

vector<vector<double>> TF_seam_12(vector<vector<double>> coord, vector<vector<double>> coord_l, vector<vector<double>> coord_r, double psi, unsigned int b)//Longitudinal total force
{
	vector<vector<double>> sf(coord.size());
	vector<vector<double>> rf(coord.size());
	vector<vector<double>> result(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		sf[i].resize(6);
		rf[i].resize(6);
		result[i].resize(6);
	}

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		sf[i] = SF_long(coord, psi)[i] + seam_12(coord, coord_l, coord_r, psi, b)[i];
	}

	rf = RF(coord);

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		//result[i] = zero(6) - sf[i];
	    result[i] = rf[i] - sf[i];
	}

	return result;
}

vector<vector<double>> evaluate(vector<vector<double>> coord, vector<vector<double>> coord_l, vector<vector<double>> coord_r, double psi, unsigned int b)//Equations
{
	vector<vector<double>> result(coord.size());
	vector<vector<double>> force(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
		force[i].resize(6);
	}

	force = TF(coord, coord_l, coord_r, psi, b);

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i] = coord[i] + force[i];
	}

	return result;
}

vector<vector<double>> ev_seam_0(vector<vector<double>> coord, vector<vector<double>> coord_l, vector<vector<double>> coord_r, double psi, unsigned int b)//Equations
{
	vector<vector<double>> result(coord.size());
	vector<vector<double>> force(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
		force[i].resize(6);
	}

	force = TF_seam_0(coord, coord_l, coord_r, psi, b);

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i] = coord[i] + force[i];
	}

	return result;
}

vector<vector<double>> ev_seam_12(vector<vector<double>> coord, vector<vector<double>> coord_l, vector<vector<double>> coord_r, double psi, unsigned int b)//Equations
{
	vector<vector<double>> result(coord.size());
	vector<vector<double>> force(coord.size());

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i].resize(6);
		force[i].resize(6);
	}

	force = TF_seam_12(coord, coord_l, coord_r, psi, b);

	for (unsigned int i = 0; i < coord.size(); i++)
	{
		result[i] = coord[i] + force[i];
	}

	return result;
}



void data_file(vector<vector<vector<double>>> data, string filename)
{
	ofstream results;
	results.open(filename);

    	for (unsigned int i = 0; i < data.size(); i++)
	{
		for (unsigned int k = 0; k < data[i][0].size(); k++)
		{
			results << 'Q' << i + 1 <<  k + 1 << ',';
		}
    }
    results << endl;
    
    unsigned int cntr = 0;
	
    for (unsigned int j = 0; j < data[cntr].size(); j++)
	{
		for (unsigned int i = 0; i < data.size(); i++)
		{
			for (unsigned int k = 0; k < data[i][j].size(); k++)
			{
				if (k < 3) { results << data[i][j][k] * pow(10, 9.) << ','; }
				else { results << data[i][j][k] << ','; }
			}
		}
        cntr++;
        results << endl;
	}
	results.close();
}

void data_screen(vector<vector<vector<double>>> data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		cout << "PF #" << i + 1 << endl;
		for (unsigned int k = 0; k < data[i][0].size(); k++)
		{
			cout << 'Q' << k + 1 << ' ';
		}
		cout << endl;
		for (unsigned int j = 0; j < data[i].size(); j++)
		{
			for (unsigned int k = 0; k < data[i][j].size(); k++)
			{
				if (k < 3) { cout << data[i][j][k] * pow(10, 9.) << ' '; }
				else { cout << data[i][j][k] << ' '; }
			}
			cout << endl;
		}
	}
}


int main()
{
    unsigned int b = 1;
    unsigned int k = 0;
    unsigned int nsize = 13;
    vector<unsigned int> N(nsize);
    N = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

	//cout<<N.size()<<endl;

	vector<vector<vector<double>>> init_coord(nsize);
	vector<vector<vector<double>>> coord(nsize);
	vector<vector<vector<double>>> coord_backup(nsize);

	for (int i = 0; i < nsize; i++)
	{
		init_coord[i].resize(N[i]);
		coord[i].resize(N[i]);
		coord_backup[i].resize(N[i]);
		for (unsigned int j = 0; j < N[i]; j++)
		{
			init_coord[i][j].resize(6);
			coord[i][j].resize(6);
			coord_backup[i][j].resize(6);
		}
	}

	for (int i = 0; i < nsize; i++)
	{
		for (unsigned int j = 0; j < N[i]; j++)
		{
			init_coord[i][j] = { mt_rad * cos( i * 2. * pi / 13.), mt_rad * sin( i * 2. * pi / 13.), (2. * j + 1.) * R + i * 12.0 / 13.0 * pow(10., -9.),0,0,0 };
		}
	}

	coord = init_coord;

	int ts_amount;
	cout << "Enter the amount of steps." << endl;
	cin >> ts_amount;
	double T = ts_amount * dt;
	double psi;

	cout << "Observation time is going to be " << T << "sec." << endl;
    double cnt = 1.;
	double d1;
	double d2;

	for (int time = 0; time < ts_amount; time++)
	{
		coord_backup = coord;
		for (int i = 0; i < nsize; i++)
		{
			 if (i == 0) {d1 = 12; d2 = i+1;}
			 if (i == nsize - 1) { d1 = i-1, d2 = 0;}
			 if ((i != 0) && (i != nsize - 1)) {d1 = i - 1; d2 = i + 1;}
			psi = i * 2. * pi / 13;
			 coord[i] = evaluate(coord_backup[i], coord_backup[d1], coord_backup[d2], psi, b);
			 if (i == 0) {coord[i] = ev_seam_0(coord_backup[i], coord_backup[d1], coord_backup[d2], psi, b);}
			 if (i == nsize-1) {coord[i] = ev_seam_12(coord_backup[i], coord_backup[d1], coord_backup[d2], psi, b);}
			//coord[i] = ev(coord[i],psi);
			coord[i][0] = init_coord[i][0];

		}
		cout << cnt << endl;
        cnt = cnt + 1.;
	}

	cout << "Enter the name of the file:" << endl;
	string filename;

	cin >> filename;

	data_screen(coord);
    data_file(init_coord, "init_coord.csv");
    data_file(coord, filename);
}
