
#include <algorithm>
#include <vector>
#include <functional>
#include <iostream>

#include "unpickle.h"

using namespace std;

// This class is required to bind arguments for stable_partitioning()
// since stable partitioning does not accept additional arguments
template<class T1, class T2>
class JudgePart:public binary_function<T1, T2, bool>
{
public:
	bool operator()(T1 val, T2 pivot) const {
		return (val <= pivot);
	}
};

double calc_span(vector<double> *p_coord, int axis){

	// Extract specified axis data
	vector<double> axis_coord;
	for (unsigned int i=0; i< p_coord->size(); i++)
		axis_coord.push_back(p_coord[axis].at(i));

	vector<double>::iterator min = min_element(axis_coord.begin(), axis_coord.end());
	vector<double>::iterator max = max_element(axis_coord.begin(), axis_coord.end());


//	double min = min_element(p_coord[axis].begin(), p_coord[axis].end());
//	double max = max_element(p_coord[axis].begin(), p_coord[axis].end());

	return (*max-*min);
}

double find_pivot(vector<double> *p_coord, int axis){

	// Extract specified axis data
	vector<double> axis_coord;
	for (unsigned int i=0; i< p_coord->size(); i++)
		axis_coord.push_back(p_coord[axis].at(i));

	// Sort axis data
	sort(axis_coord.begin(), axis_coord.end());

	// find mid point
	int len = axis_coord.size();
	if(len % 2 == 0){
	// even
		int right = len / 2;
		int left = right - 1;
		return (axis_coord.at(right) + axis_coord.at(left))/2.0;

	}else{
	// odd
		return axis_coord.at(len/2);
	}

}

void rcb_rec(vector<double> *p_coord, int depth){

	// end of partitioning
	if(depth == 0) return;

	// calculate max distance on each axis
	double x_span = calc_span(p_coord, 0);
	double y_span = calc_span(p_coord, 1);
	double z_span = calc_span(p_coord, 2);

	// choose axis
	int axis = -1;
	if( (x_span >= y_span) && (x_span >= z_span) )
		axis = 0;
	else if( (y_span >= x_span) && (y_span >= z_span) )
		axis = 1;
	else if( (z_span >= x_span) && (z_span >= y_span) )
		axis = 2;
	else
		return;

	// find mid-point
	double pivot = find_pivot(p_coord, axis);
	cout << "pivot=" << pivot;

	// partition into two
	vector<double>::iterator bound;
	bound = stable_partition(p_coord->begin(), p_coord->end(), bind2nd(JudgePart<double, double>(), 2.0));

	// set label


	// next partitioning
//	rcb_rec(depth-1);
//	rcb_rec(depth-1);

}


int main(int argc, char **argv){

	vector<double> coord;

	// Read input coordinate data
	string filename(argv[1]);
	struct params *p = parse_file(filename);


 	// 
	rcb_rec(&coord,3);

}


