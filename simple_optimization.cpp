#include <iostream>
#include <algorithm>
#include <vector>
#include <bitset>
#include <cstdlib>
#include <time.h>
#include <cmath>
using namespace std;

const unsigned n = 22;
const int pop_size = 500;
const double pi = acos(-1.0);
const double p_mutate = 0.01;
const double p_crossover = 0.45;

double
f(double x)
{
    return x*sin(10*pi*x) + 1.0;
}

double
transform_x(bitset<n> chrom)
{
    unsigned x_prime = 0;
    for(unsigned i = 0; i < n; ++i){
        x_prime += chrom[i]*(1 << i);
    }
    double var = x_prime * 3.0/((1 << n) - 1);
    return -1.0 + var;
}

double
eval(bitset<n> chrom)
{
    double x = transform_x(chrom);
    return f(x);
}

void
pop_init(bitset<n>* pop)
{
    srand (time(NULL));
    for(int j = 0; j < pop_size; ++j){
        for(int i = 0; i < n; ++i){
            pop[j].set(i,rand()%2);
        }
    }
}

bitset<n>
mutate(bitset<n> chrom)
{
    double vn = (double) rand()/RAND_MAX;
    if(vn <= p_mutate){
        int rp = rand() % n;
        chrom.set(rp, chrom[rp]^(1 << rp));
    }
    return chrom;
}

void
crossover(bitset<n>& chrom1, bitset<n>& chrom2)
{
    bitset<17> sub_chrom1(chrom1.to_string(),5,17);
    bitset<17> sub_chrom2(chrom2.to_string(),5,17);

    for(int i = 0; i < 17; ++i){
        chrom1.set(i,sub_chrom2[i]);
        chrom2.set(i,sub_chrom1[i]);
    }
}

bitset<n>
selection(bitset<n>* a)
{
    for(int i = 1; i < pop_size; ++i){
        if(eval(a[i - 1]) > eval(a[i])) a[i] = a[i - 1];
    }
    return a[pop_size - 1];
}

int
main(void)
{
    srand (time(NULL));

    bitset<n>* pop = new bitset<n>[pop_size];
    pop_init(pop);

    int t = 0;
    double T = 10.0;
    bitset<n> curr = selection(pop);
    double curr_eval = eval(curr);
    double thresh;
    while(t < 15.0){
        int pos = rand() % pop_size;
        while(T > 0.0){
            curr = mutate(curr);
            bitset<n> next_cand = pop[pos];
            if(eval(curr) < eval(next_cand)){
                curr = next_cand;
            } else {
                thresh = rand()/RAND_MAX;
                if(thresh < exp((eval(curr) - eval(next_cand))/T)) curr = next_cand;
                //if(thresh <= p_crossover) crossover(curr, next_cand);
            }
            T -= 0.1;
        }
        t += 1.0;
    }
    cout << eval(curr) << endl;
}
