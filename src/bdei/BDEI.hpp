#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <array>
#include <nlopt.hpp>
#include <algorithm>
#include <string>

#include <chrono>
#include <random>
#include "pool.hpp"

using namespace std::chrono;
using namespace std;


typedef double R;

struct Solution {
    R likelihood;
    R mu, la, psi, p;
    R mu_min, la_min, psi_min, p_min;
    R mu_max, la_max, psi_max, p_max;
    R cpu_time;
    int nb_iter;

    Solution(R l, R mu_, R la_, R psi_, R p_, R cpu_time_, int nb_iter_);
    Solution(R l, R mu_, R mu_min_, R mu_max_, R la_, R la_min_, R la_max_, R psi_, R psi_min_, R psi_max_, R p_, R p_min_, R p_max_, R cpu_time_, int nb_iter_);
};


Solution
*inferParameters(const string &treename, R *x0, const R *dub, R pie, R mu, R lambda, R psi, R p, R u, R ut, int nbdirerr, int size_pool, int debug_, int nstarts);

R calculateLikelihood(const string &treename, R mu, R lambda, R psi, R p, R pie, R u, R ut, int size_pool, int debug_);