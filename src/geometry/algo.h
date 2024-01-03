#ifndef ALGO_H
#define ALGO_H

#include <iostream>
#include <vector>

#include "mymesh.h"

// mesh collision detection
std::vector<std::pair<std::string, double>> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points);

//bb collision detection
std::vector<std::string> collision_detection_bb(std::vector<Mymesh> &organ, Mymesh &tissue);

std::vector<std::pair<int, double>> collision_detection_single_tissue_2(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points);
#endif