#ifndef __COMMUNICATE_H
#define __COMMUNICATE_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include "model.h"


void print_info(const model &Sphere,int MCs, double stepsize);
void save_info(std::ofstream &file,const model &Sphere,int MCs, double stepsize);

#endif