#pragma once
#include <iostream>
using namespace std;


int seed = 12345678;

struct vector3
{
  double x, y, a;
};

float r4_uniform_ab( float a, float b, int &seed )
{
  const int i4_huge = 2147483647;
  int k;
  float value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;
  // cout << k << endl;
  seed = 16807 * ( seed - k * 127773 ) - k * 2836;
  // cout << seed << endl;
  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  value = (double) (seed) * 4.656612875E-10;
  value = a + ( b - a ) * value;
  return value;
};

int i4_uniform_ab( int a, int b, int &seed )
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;
  seed = 16807 * ( seed - k * 127773 ) - k * 2836;
  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 )
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }
  return value;
};

struct vector2
{
  double val1, val2;
};

vector2 ClctDistance(vector3& loc1, vector3& loc2)
{

  double ds = sqrt( (loc1.x-loc2.x)*(loc1.x-loc2.x) + (loc1.y-loc2.y)*(loc1.y-loc2.y) );
  double da = loc1.a - loc2.a;
  vector2 distances = { ds, da };
  return distances;
};

vector2 normal_distr(double mu, double sigma)
{
    double u1 = r4_uniform_ab(0, 1, seed);
    double u2 = r4_uniform_ab(0, 1, seed);
    // double z1 = sqrt(-2*log(u1))*cos(2*M_PI*u2);
    // double z2 = sqrt(-2*log(u1))*sin(2*M_PI*u2);
    vector2 norm_vals = { mu + sigma*sqrt(-2*log(u1))*cos(2*M_PI*u2), mu + sigma*sqrt(-2*log(u1))*sin(2*M_PI*u2)};
    return norm_vals;
}

void LoadNoise(std::vector<vector2> &eta, int kid, int start, int CHUNK)
{
  std::string index_x, index_y; std::stringstream num_x, num_y;
  num_x << 2*kid;
  num_x >> index_x;
  std::ifstream infile("../eta" + index_x + ".csv");
  string line;
  double noise;
  int line_num = 0;
  while (infile.good())
  {
    std::getline(infile, line);
    // std::cout << line_num << std::endl;
    if (line_num >= start)
    {
      noise = stod(line);
      eta[line_num-start].val1 = noise;
    }
    line_num ++;
    if (line_num >= start+CHUNK) { break; }
  }

  num_y << 2*kid + 1;
  num_y >> index_y;
  // cout << num_y <<endl;
  std::ifstream infiley("../eta" + index_y + ".csv");
  // cout << "./eta" + index_y + ".csv" << endl;
  line_num = 0;
  while (infiley.good())
  {
    std::getline(infiley, line);
    // std::cout << line_num << std::endl;
    if (line_num >= start)
    {
      noise = stod(line);
      eta[line_num-start].val2 = noise;
    }
    line_num ++;
    if (line_num >= start+CHUNK) { break; }
  }
}

