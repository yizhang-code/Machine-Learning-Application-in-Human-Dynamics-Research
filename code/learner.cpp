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


enum room_size{ KIDS=10, M=1};

const double LENGTH=8.6;
const double WIDTH=4.8;
const int STEPS = 300000;
const int UNIT_TIME=100;


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

struct Params
{
	double D, A, R0;
	double EPS, SIG, N;
};

class Obs
{
	vector<User> loc;
	User user;

public:
	void add(double x, double y, double theta)
	{
		user.x = x;
		user.y = y;
		user.a = theta;
		loc.push_back(user);
	}

	int getSize()
	{
		return loc.size();
	}

	User findObs(int step)
	{
		if (step >= loc.size() || step < 0)
		{
			cout << "err ---- " << "length is " << loc.size() << '\t' <<step << " is not out of the range!!" << endl;
		}
		return loc[step];
	}
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

void ModelPotential(vector2& distances, const Params& params, vector2& expTerm)
{
	double A = params.A; double R0 = params.R0;
	expTerm.val1 = exp(-2 * A * (distances.val1 - R0));
	expTerm.val2 = exp(-A * (distances.val1 - R0));
};

double SoftPotential(double ds, const Params& params)
{
	double EPS = params.EPS; double SIG = params.SIG; double N = params.N;
	if (ds <= SIG){ return EPS * pow(SIG / ds, N); } else { return 0; };
};

void clct_dis(const User& user1, const User& user2, vector2& distances)
{
	double dx = user1.x-user2.x; double dy = user1.y-user2.y;
	distances.val1 = sqrt(dx*dx + dy*dy);
	distances.val2 = user1.a - user2.a;

};

void TotalPotential(const vector<User>& users, const Params& params, Params& deriv)
{
	double soft_pot, cos_theta;
	vector2 distances, expTerm;
	deriv.D = 0; deriv.A = 0; deriv.R0 = 0;
	deriv.EPS = 0; deriv.SIG = 0; deriv.N = 0;
	for (int i = 0; i < KIDS; i ++)
	{
		for (int j = i+1; j < KIDS; j ++)
		{
			User user_i = users[i]; User user_j = users[j];
			clct_dis(user_i, user_j, distances);
			soft_pot  = SoftPotential(distances.val1, params);
			ModelPotential(distances, params, expTerm);
			cos_theta = cos(distances.val2);

			deriv.D -= (expTerm.val1 - 2*expTerm.val2)*cos_theta;
			deriv.A += 2*params.D*params.A*(expTerm.val1*(distances.val1 - params.R0) + expTerm.val2*(params.R0-distances.val1))*cos_theta;
			deriv.R0 -= 2*params.D*params.A*(expTerm.val1 - expTerm.val2)*cos_theta;

			deriv.EPS += soft_pot/params.EPS;
			deriv.SIG += params.N*soft_pot/(params.SIG/distances.val1);
			deriv.N += soft_pot*log(params.SIG/distances.val1);

		}
	}
};

void output(vector<User>& users)
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


double ClctPotential(const vector<User>& users, int rand_kid_idx, User& random_mover, const Params& params)
{
	double potential = 0; double pair_potential, soft_pot;
	vector2 distances, expTerm;

	for (int i = 0; i < KIDS; ++i)
	{
		if (i != rand_kid_idx)
		{
			User other_kid = users[i];
			clct_dis(random_mover, other_kid, distances);
			soft_pot  = SoftPotential(distances.val1, params);
			ModelPotential(distances, params, expTerm);
			pair_potential = soft_pot - params.D*(expTerm.val1 - 2*expTerm.val2)*cos(distances.val2);
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
	if (exp(-dV) > p) { return true; }
	else { return false; }
};



void sliceObs(vector<Obs> &observations, int obs_indx, vector<User>& snapshot)
{
	for (int kid = 0; kid < KIDS; ++kid)
	{
		snapshot[kid] = observations[kid].findObs(obs_indx);
	}
}

void update(vector<User> &users, vector<Obs> &observations)
{
	double dx, dy, da, curr_potential, possible_potential, tot_potential;
	User random_mover; int rand_kid_idx;
	int obs_num = observations[0].getSize();

	Params params, deriv_obs, deriv;
	params.D = 15; params.A = 1.3; params.R0 = 0.9;
	params.EPS = 10; params.SIG = 0.5; params.N = 2;

	double alpha; double dD, dA, dR0, dEPS, dSIG, dN;

	ofstream save0("learn.csv");
	int rand_obs_indx; vector<User> snapshot(KIDS);

	for (int i=0; i<STEPS; i++)
	{
		alpha = 0.002/(1+2*i/15000.0);

		rand_obs_indx = i4_uniform_ab(0, obs_num-1, seed);
		sliceObs(observations, rand_obs_indx, snapshot);

		TotalPotential(snapshot, params, deriv_obs);
		TotalPotential(users, params, deriv);

		dD = deriv_obs.D - deriv.D; dA = deriv_obs.A - deriv.A; dR0 = deriv_obs.R0 - deriv.R0;
		// dEPS = deriv_obs.EPS - deriv.EPS; dSIG = deriv_obs.SIG - deriv.SIG; dN = deriv_obs.N - deriv.N;
		params.D -= alpha*dD;
		params.A -= alpha*dA/100;
		params.R0 -= alpha*dR0/500;
		// params.EPS -= alpha*dEPS;

		// cout << params.D << '\t' << dD << endl;
		// cout << params.A << '\t' << dA << endl;
		// cout << params.R0 << '\t' << dR0 << endl;
		// cout << params.EPS << '\t' << dEPS << endl;

		if (i%UNIT_TIME == 0) save0 << params.D << ',' << params.A << ',' << params.R0 << endl;

		rand_kid_idx = i4_uniform_ab(0, KIDS-1, seed); // sample a random kid using random distribution
		curr_potential = ClctPotential(users, rand_kid_idx, users[rand_kid_idx], params);

		dx = r4_uniform_ab(-0.1, 0.1, seed), dy = r4_uniform_ab(-0.1, 0.1, seed), da = r4_uniform_ab(-5*M_PI/180, 5*M_PI/180, seed);
		Variable update_x(users[rand_kid_idx].x+dx, "x"), update_y(users[rand_kid_idx].y+dy, "y"), update_a(users[rand_kid_idx].a+da, "a");
		random_mover.x = update_x.get_val();
		random_mover.y = update_y.get_val();
		random_mover.a = update_a.get_val();

		possible_potential = ClctPotential(users, rand_kid_idx, random_mover, params);

		if (UpdateStatus(possible_potential-curr_potential)) users[rand_kid_idx] = random_mover;
	}
};



void LoadObservation(vector<Obs> &observations)
{
	for (int kid = 0; kid < KIDS; ++kid)
	{
		Obs kid_obs;
		string label = to_string(kid);
		ifstream infile("./info_" + label +".csv");
		double x, y, theta; char c1, c2; string line;
		while (getline(infile, line))
		{
			std::istringstream ss(line);
			if (line.size())
			{
				ss >> x >> c1 >> y >> c2 >> theta;
				kid_obs.add(x, y, theta);
			}
		}
		observations.push_back(kid_obs);
	}
}


// void InitParameters(char* argv[])
// {
	// STEPS=atoi(argv[1]);
	// UNIT_TIME=atoi(argv[2]);
// };

int main(int argc, char* argv[])
{

	// InitParameters(argv);
	vector<User> users = InitUsers();
	vector<Obs> observations;
	LoadObservation(observations);
	update(users, observations);
	std::cout << "Completed!" << std::endl;
};
