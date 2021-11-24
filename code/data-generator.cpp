#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <string.h>
#include <random>
#include<cmath>
#include <stdlib.h>
#include "func.h"


enum room_size{ KIDS=10, M=1, K_B=1};

const double LENGTH=8.6;
const double WIDTH=4.8;
static int STEPS = 25000000;
static int UNIT_TIME=1000;

const double T=0.02;

const double D=0.2;
const double EPSI=0.2;

const double A=1.1;
const double B = 1.1*0.7; // initial value of A*R0
const double R0=B/A;

const double Rc = 1/A;
const double SIGMA=0.5;
const double N=2;



class Variable
{
private:
	double val;
	const char* type;
public:
	double get_val() { return val; };
	Variable()
	{
		std::cout << "Create an empty variale!" << std::endl;
	};
	Variable(double val1, const char* type1)
	{
		val = periodic_boundary(val1, type1);
	}
	double periodic_boundary(double var_val, const char* unknow_type)
	{
		if (!strcmp(unknow_type, "x")) check1D(var_val, LENGTH, 0);
		if (!strcmp(unknow_type, "y")) check1D(var_val, WIDTH, 0);
		if (!strcmp(unknow_type, "a")) check1D(var_val, M_PI, -M_PI);
		return var_val;
	}

	void check1D(double &var, double upper_bound, double lower_bound)
	{
		var -= lower_bound;
		var = fmod(var, upper_bound - lower_bound);
		var += lower_bound;
		if (var < lower_bound) var += (upper_bound - lower_bound);
	}

	Variable operator+(const Variable &variable2)
	{
		Variable variable(this->val + variable2.val, this->type);
		return variable;
	};
};

struct User
{
	double x, y, a;
};

std::vector<User> InitUsers()
{
	std::ifstream file("./ini_status.dat");
	std::vector<User> users;
	std::string line;
	User user;
//	while ( std::getline(file, line) )
	for (int i=0; i<KIDS; ++i)
	{
		std::getline(file, line);
		std::istringstream is(line);
		if (line.size())
		{
			vector<string> tokens;
			copy(istream_iterator<string>(is),
		    istream_iterator<string>(),
		    back_inserter(tokens));
		    Variable var_x(stod(tokens[0]), "x"), var_y(stod(tokens[1]), "y"), var_a(stod(tokens[2]), "a");
		    user.x = var_x.get_val();
		    user.y = var_y.get_val();
		    user.a = var_a.get_val();

		    users.push_back(user);
		}
	}
	return users;
};

double ModelPotential(vector2& distances)
{
	return - D * (exp(-2 * A * (distances.val1 - R0)) - 2 * exp(-A * (distances.val1 - R0))) * cos(distances.val2);
};

double SoftPotential(double ds)
{
	if (ds <= SIGMA){ return EPSI * pow(SIGMA / ds, N); } else { return 0; };
};

vector2 clct_dis(User& user1, User& user2)
{
	double dx = user1.x-user2.x; double dy = user1.y-user2.y;
	// if (dx > LENGTH/2) dx = LENGTH - dx;
	// if (dy > WIDTH/2) dy = WIDTH - dy;

	vector2 distances;
	distances.val1 = sqrt(dx*dx + dy*dy);
	distances.val2 = user1.a - user2.a;
	return distances;
};

double TotalPotential(std::vector<User>& users)
{
	double potential = 0;
	double pair_potential;
	for (int i = 0; i < KIDS; i ++)
	{
		for (int j = i+1; j < KIDS; j ++)
		{
			User user_i = users[i];
			User user_j = users[j];
			vector2 distances = clct_dis(user_i, user_j);
			cout << distances.val1 << '\t' << distances.val2 << endl;
			pair_potential = ModelPotential(distances) + SoftPotential(distances.val1);
			potential += pair_potential;
			if (potential > 10)
			{
				cout << i << '\t' << j << '\t' << user_i.x << '\t' << user_i.y<< '\t' << user_j.x << '\t' << user_j.y << endl;
			}
		}
	}
	return potential;
};

void output(std::vector<User>& users)
{
	for (int i = 0; i < users.size(); ++i)
	{
		std::string index;
		std::stringstream num;
		num << i;
		num >> index;
		std::ofstream save("info_" + index + ".csv", std::ios_base::app);
		save << users[i].x << "," << users[i].y << "," << users[i].a << std::endl;
	}
};


double ClctPotential(std::vector<User>& users, int rand_kid_idx, User& random_mover)
{
	double potential = 0;
	double pair_potential;
	for (int i = 0; i < KIDS; ++i)
	{
		if (i != rand_kid_idx)
		{
			User other_kid = users[i];
			vector2 distances = clct_dis(random_mover, other_kid);
			pair_potential = ModelPotential(distances) + SoftPotential(distances.val1);
			potential += pair_potential;
		}
	}
	return potential;
};


bool UpdateStatus(double dV)
{
	if (dV < 0) { return true; }
	float p;
	p = r4_uniform_ab(0, 1, seed);
	if (exp(-dV/T) > p) { return true; }
	else { return false; }
};

void update(std::vector<User>& users)
{
	cout << "here" << 't' << STEPS << endl;
	double dx, dy, da, curr_potential, possible_potential, tot_potential;
	int rand_kid_idx;
	User random_mover;
	string T_str = to_string(T);
	ofstream save0 ("energy" + T_str + ".csv");
	int OUT = 10000000;
	for (int i=0; i<STEPS; i++)
	{
		if (i%1000000==0)
		{
			cout << i << endl;
		}
		if (i%UNIT_TIME == 0 && i > OUT)
		{
			output(users);
			// for (int i = 0; i < KIDS; ++i)
			// {
			// 	cout << users[i].x << '\t'  << users[i].y << '\t' << users[i].a << '\t' << endl;
			// }
			// tot_potential = TotalPotential(users);
			// cout << tot_potential << endl;
			// save0 << tot_potential << endl;
		}

		rand_kid_idx = i4_uniform_ab(0, KIDS-1, seed); // sample a random kid using random distribution
		curr_potential = ClctPotential(users, rand_kid_idx, users[rand_kid_idx]);

		dx = r4_uniform_ab(-0.1, 0.1, seed), dy = r4_uniform_ab(-0.1, 0.1, seed), da = r4_uniform_ab(-5*M_PI/180, 5*M_PI/180, seed);
		Variable update_x(users[rand_kid_idx].x+dx, "x"), update_y(users[rand_kid_idx].y+dy, "y"), update_a(users[rand_kid_idx].a+da, "a");
		random_mover.x = update_x.get_val();
		random_mover.y = update_y.get_val();
		random_mover.a = update_a.get_val();

		possible_potential = ClctPotential(users, rand_kid_idx, random_mover);

		if (UpdateStatus(possible_potential-curr_potential))
		{
			// double e1 = TotalPotential(users);
			users[rand_kid_idx] = random_mover;
		}
	}
};

void InitOutput()
{
	for (int i = 0; i < KIDS; ++i)
	{
		std::string index;
		std::stringstream num;
		num << i;
		num >> index;
		std::ofstream save("info_" + index + ".csv");
	}
};

void InitParameters(char* argv[])
{
	// STEPS=atoi(argv[1]);
	// UNIT_TIME=atoi(argv[2]);
};

int main(int argc, char* argv[])
{

	// InitParameters(argv);
	InitOutput();
	std::vector<User> users = InitUsers();
	// std::vector<double> potential = ClctPotential(users);
	update(users);
	std::cout << "Completed!" << std::endl;
};
