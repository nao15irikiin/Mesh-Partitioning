
#include <algorithm>
#include <vector>
#include <functional>
#include <iostream>
#include <stdio.h>
#include <pthread.h>
#include <omp.h>

#include "op_part_rcb.h"
#include "../common/simple_timer.h"

using namespace std;

// Struct Declaration
template <class T>
struct thdata{
	pthread_t th;
	vector<T> *p_coord;
	vector<int> *p_map;
	vector<int> *p_part;
	int start;
	int end;
	int dims;
	int cur_lvl;
	int start1;
	int end1;
	int start2;
	int end2;
};

struct databound{
	int start;
	int end;
};


// Global Variable Declaration

SimpleTimer tmr_span("span");
SimpleTimer tmr_pivot1("pivot1");
SimpleTimer tmr_part1("part1");
SimpleTimer tmr_part2("part2");
SimpleTimer tmr_part3("part3");
SimpleTimer tmr_out("out");


template <typename T>
double calc_span(vector<T> *p_coord, int start, int end, int dims, int axis){

	if(dims==0){
		cout << "specified dims is 0" << endl; exit(-1);
	}

	// Set initial value
	T minval, maxval;
	minval = p_coord->at(start*dims+axis); maxval = minval;

	// Search min/max
	for (int i=start; i<end; i++){
		if (minval > p_coord->at(i*dims+axis)) minval = p_coord->at(i*dims+axis);
		if (maxval < p_coord->at(i*dims+axis)) maxval = p_coord->at(i*dims+axis);
	}

	//cout << "    @axis:" << axis << " min=" << minval << " max=" << maxval << endl;

	T span = maxval-minval;
	return (span>0)?span:(T)-1*span;
}


template <typename T>
T find_pivot(vector<T> *p_coord, int start, int end, int dims, int axis){

	if(dims==0){
		cout << "specified dims is 0" << endl; exit(-1);
	}

	// Extract specified axis data
	vector<T> axis_coord(end-start);
	int cnt=0;
	for (int i=start; i< end; i++){
		axis_coord[cnt]=(p_coord->at(i*dims+axis));
		cnt++;
	}

	// Sort and find mid point
	//sort(axis_coord.begin(), axis_coord.end());
	int mid = axis_coord.size()/2;
	nth_element(axis_coord.begin(), axis_coord.begin()+mid, axis_coord.end());

	return axis_coord.at(mid);
}

template <typename T>
int partition(vector<T> *p_in, 
              vector<int> *p_map,
              vector<int> *p_part,
              int start, 
              int end, 
              int dims, 
              int axis, 
              T pivot,
              int level){

/*
tmr_part1.start();
	// Buffer
	vector<T> out1, out2;
	vector<int> map1, map2;
	int cnt1=0, cnt2=0;
	for(int i=start; i < end; i++){
		if(p_in->at(i*dims+axis)> pivot){
			for(int d_index=0; d_index<dims; d_index++)
				out1.push_back(p_in->at(i*dims+d_index));
			map1.push_back(p_map->at(i));
			cnt1+=1;
		}else{
			for(int d_index=0; d_index<dims; d_index++)
				out2.push_back(p_in->at(i*dims+d_index));
			map2.push_back(p_map->at(i));
			cnt2+=1;
		}
	}

tmr_part1.stop_and_add_to_total();
*/


tmr_part1.start();

	// Buffer
	vector<T> out1, out2;
	vector<int> map1, map2;
	int cnt1=0, cnt2=0;
	int halfcnt = (end-start)/2;
	map1.reserve(halfcnt);
	map2.reserve(halfcnt);
	out1.reserve(halfcnt*dims);
	out2.reserve(halfcnt*dims);
	T tmpcoord;
	for(int i=start; i < end; i++){
	
		tmpcoord = p_in->at(i*dims+axis);

		// - balance load when multiple node has same coordinate value
		if(tmpcoord == pivot){
			// if cnt1 is less than half size, data belong to map1
			if(halfcnt>cnt1){
				for(int d_index=0; d_index<dims; d_index++)
					out1.push_back(p_in->at(i*dims+d_index));
				map1.push_back(p_map->at(i));
				cnt1+=1;
			// if cnt1 is more than half size, data belong to map2
			}else{
				for(int d_index=0; d_index<dims; d_index++)
					out2.push_back(p_in->at(i*dims+d_index));
				map2.push_back(p_map->at(i));
				cnt2+=1;
			}
		}else if(tmpcoord > pivot){
			for(int d_index=0; d_index<dims; d_index++)
				out1.push_back(p_in->at(i*dims+d_index));
			map1.push_back(p_map->at(i));
			cnt1+=1;
		}else{
			for(int d_index=0; d_index<dims; d_index++)
				out2.push_back(p_in->at(i*dims+d_index));
			map2.push_back(p_map->at(i));
			cnt2+=1;
		}
	}
tmr_part1.stop_and_add_to_total();;

	//cout << "    out1.size =" << out1.size()/dims << " out2.size =" << out2.size()/dims << endl;
	//cout << "    level =" << level  << endl;

	// Replace to original coord data <--- can be replaced with this loop to memcpy but may not safe..

tmr_part2.start();

	for(int i=0; i < (int)out1.size(); i++){
		p_in->at(start*dims+i) = out1.at(i);
	}

	for(int i=0; i < (int)out2.size(); i++){
		p_in->at(start*dims+out1.size()+i) = out2.at(i);
	}
tmr_part2.stop_and_add_to_total();;


tmr_part3.start();
	for(int i=0; i<cnt2; i++){
		p_part->at(start+cnt1+i) |= (1 << level);
	}

	// Replace to original map/part data
	for(int i=0; i < (int)map1.size(); i++)
		p_map->at(start+i) = map1.at(i);

	for(int i=0; i < (int)map2.size(); i++)
		p_map->at(start+map1.size()+i) = map2.at(i);
tmr_part3.stop_and_add_to_total();;

	return cnt1;
}


void *rcb_rec2D_ser(void *tharg){

	// convert void* to struct
	struct thdata<float> *arg = (struct thdata<float> *)tharg;

//cout << "thread num=" << omp_get_thread_num() << endl;

	//cout << "-------------------------------------------------------------" << endl;
	//cout << "    start=" << arg->start << " end=" << arg->end << endl;

	// calculate max distance on each axis
	double x_span = calc_span(arg->p_coord, arg->start, arg->end, arg->dims, 0);
	double y_span = calc_span(arg->p_coord, arg->start, arg->end, arg->dims, 1);

	// choose axis
	int axis = -1;
	if(x_span >= y_span)
		axis = 0;
	else
		axis = 1;

	//cout << "    x_span=" << x_span << " y_span=" << y_span << " axis=" << axis << endl;

	// find mid-point
	float pivot = find_pivot(arg->p_coord, arg->start, arg->end, arg->dims, axis);
	//cout << "    pivot=" << pivot <<endl;

	// partition into two
	int level= arg->cur_lvl-1;
	int part_index = partition(arg->p_coord, arg->p_map, arg->p_part, arg->start, arg->end, arg->dims, axis, pivot, level);

	//cout << "    part_index=" << part_index << endl;

	// set next partitioning start and end boundary
	arg->start1 = arg->start;
	arg->end1   = arg->start + part_index;
	arg->start2 = arg->start + part_index;
	arg->end2   = arg->end;

	return NULL;
}


template <class T>
void rcb_rec2D(vector<T> *p_coord,
               vector<int> *p_map,
               vector<int> *p_part,
               int start,
               int end,
               int dims, 
               int cur_depth,
               int ttl_depth){

	// end of partitioning
	if(cur_depth == 0) return;

	//cout << "-------------------------------------------------------------" << endl;
	//cout << "    start=" << start << " end=" << end << " cur_depth=" << cur_depth << endl;

tmr_span.start();
	// calculate max distance on each axis
	double x_span = calc_span(p_coord, start, end, dims, 0);
	double y_span = calc_span(p_coord, start, end, dims, 1);
tmr_span.stop_and_add_to_total();;

	// choose axis
	int axis = -1;
	if(x_span >= y_span)
		axis = 0;
	else
		axis = 1;

	//cout << "    x_span=" << x_span << " y_span=" << y_span << " axis=" << axis << endl;

	// find mid-point
tmr_pivot1.start();
	T pivot = find_pivot(p_coord, start, end, dims, axis);
tmr_pivot1.stop_and_add_to_total();;
	//cout << "    pivot=" << pivot <<endl;

	// partition into two
	int level= cur_depth-1;
	int part_index = partition(p_coord, p_map, p_part, start, end, dims, axis, pivot, level);

	//cout << "    part_index=" << part_index << endl;

	// next partitioning
	rcb_rec2D(p_coord, p_map, p_part, start, start+part_index, dims, cur_depth-1, ttl_depth);
	rcb_rec2D(p_coord, p_map, p_part, start+part_index, end, dims, cur_depth-1, ttl_depth);

}


bool check_partdata(int nsets,
                    vector<int> *p_map,
                    vector<int> *p_part,
                    int depth){

	cout << "-------------------------------------------------------------" << endl;

	if( (nsets != (int)p_map->size()) | (nsets != (int)p_part->size())){
		cout << "num of edges is not the same as the size of map data" << endl;
		exit(-1);
	}

	// Check output data validity
	vector<int> tmpcheck(nsets);
	for(int i=0; i<nsets; i++){
		tmpcheck.at(i)=0;
	}
	for(int i=0; i<nsets; i++){
		tmpcheck.at(p_map->at(i))+=1;
	}

	int map_error=0;
	for(int i=0; i<nsets; i++){
		if(tmpcheck.at(i)!=1){
			cout << "index[" << i <<"]=" << tmpcheck[i] << endl;
			map_error+=1;
		}
	}
	cout << "map_error=" << map_error << endl;

	int numofpart;
	numofpart = (1 << depth);
	vector<int> cnt(numofpart);
	for(int i=0; i<nsets; i++)
		cnt[p_part->at(i)]+=1;

	int total=0;
	for(int i=0; i<numofpart; i++){
//		cout << "cnt[" << i << "]=" << cnt.at(i) << endl;
		total+=cnt.at(i);
	}
	cout << "total=" << total << endl;

	return (map_error==0)?true:false;

}


void generate_partdata(int nnode, int nedge, int nbedge, int ncell,
                       vector<int> *p_map,
                       vector<int> *p_part,
                       point *partnode,
                       point *partedge,
                       point *partbedge, 
                       point *partcell,
                       int *cell, int *ecell, int *becell){

	int numofnode=4; // For now hardcoded

	// convert map from (new -> original) to (original -> new)
	vector<int> map_OrgToNew(ncell);
	for(int i=0; i<ncell; i++)
		map_OrgToNew.at(p_map->at(i)) = i;

	// cell partition data
	for(int i=0; i<ncell; i++){
		partcell[i].part = p_part->at(map_OrgToNew.at(i));
		partcell[i].index = i;
	}

	// node partition data
	for(int idx_node=0; idx_node<nnode; idx_node++){
		partnode[idx_node].part=-1;
		partnode[idx_node].index = idx_node;
	}

	int nodeID;
	int cellPart;
	for(int idx_cell=0; idx_cell<ncell; idx_cell++){

		// get cell partition number
		cellPart = partcell[idx_cell].part;

		// apply cell partition number to its nodes
		for(int idx_node=0; idx_node<numofnode; idx_node++){

			// get node ID
			nodeID = cell[idx_cell*numofnode+idx_node];

			// set cell partition number to node
			partnode[nodeID].part=cellPart;
		}
	}

	for(int idx_node=0; idx_node<nnode; idx_node++){
		// check node partition data
		if(partnode[idx_node].part==-1)
			cout << "partnode[" << idx_node << "] does not have partition number!" << endl;
	}

	// edge partition data
	for(int i=0; i<nedge; i++){
		//int part1 = partcell[ ecell[i*2] ].part;
		int part2 = partcell[ ecell[i*2+1] ].part;
		int target_part = part2;
		
		partedge[i].part = target_part;
		partedge[i].index = i;
	}

	// bedge partition data
	for(int i=0; i<nbedge; i++){
		int part = partcell[ becell[i] ].part;
		
		partbedge[i].part = part;
		partbedge[i].index = i;
	}
}

template <class T>
void calc_cellcentre(int ncell,
                     int dims,
                     int *cell,
                     T* x,
                     vector<T> *p_coord){

	// For now hardcoded
	int numofnode=4;
	int nodeID;

	for(int idx_cell=0; idx_cell<ncell; idx_cell++){

		for(int idx_node=0; idx_node<numofnode; idx_node++){

			// get node ID
			nodeID = cell[idx_cell*numofnode+idx_node];

			// add coordinate values
			for(int idx_dim=0; idx_dim<dims; idx_dim++)
				p_coord->at(idx_cell*dims+idx_dim) += x[nodeID*dims+idx_dim];
		}

		// divide by num of node
		for(int idx_dim=0; idx_dim<dims; idx_dim++)
			p_coord->at(idx_cell*dims+idx_dim) /= (T)numofnode;
	}

	// debug print
	T minX, maxX, minY, maxY;
	minX=p_coord->at(0);
	maxX=minX;
	minY=p_coord->at(1);
	maxY=minY;
	for(int idx_cell=0; idx_cell<ncell; idx_cell++){
		if(minX>p_coord->at(2*idx_cell)) minX=p_coord->at(2*idx_cell);
		if(maxX<p_coord->at(2*idx_cell)) maxX=p_coord->at(2*idx_cell);
		if(minY>p_coord->at(2*idx_cell+1)) minY=p_coord->at(2*idx_cell+1);
		if(maxY<p_coord->at(2*idx_cell+1)) maxY=p_coord->at(2*idx_cell+1);
	}
	cout << "minX=" << minX << " maxX =" << maxX << endl;
	cout << "minY=" << minY << " maxY =" << maxY << endl;
}

template <class T>
void part_rcb_DebugPrint(int ncell, int dims,
                         vector<T> *coord_cell,
                         vector<int> *map,
                         vector<int> *part){

	FILE *fp;
	fp = fopen("map_mp.dat","w");
	for(int i=0; i<ncell; i++){
		fprintf(fp, "%d\n", map->at(i));
	}
	fclose(fp);

	fp = fopen("part_mp.dat","w");
	for(int i=0; i<ncell; i++){
		fprintf(fp, "%d\n", part->at(i));
	}
	fclose(fp);

	fp = fopen("coord_mp.dat","w");
	for(int i=0; i<ncell*dims; i++){
		fprintf(fp, "%lf\n", (float)coord_cell->at(i));
	}
	fclose(fp);
}

void part_rcb(int numoflevel, int dims,
              int nnode, int nedge, int nbedge, int ncell,
              point *partnode, point *partedge, point *partbedge, point *partcell,
              int *cell, int *ecell, int *becell, int *edge,
              float *x){

	// calc coordinate of center of gravity in each cell
	vector<float> coord_cell(ncell*dims);
	calc_cellcentre(ncell, dims, cell, x, &coord_cell);

	// initialize map and partition data
	vector<int> map, part;
	for(int i=0; i<ncell; i++){
		map.push_back(i);
		part.push_back(0);
	}

 	// call recursive coordinate bisection algorithm
	// *rewrite recursive algorithm to for loop
	int partnum = 1<<numoflevel;
	struct thdata<float> *tharg = (struct thdata<float> *) malloc(partnum*sizeof(struct thdata<float>));
	struct databound *nxtdat = (struct databound *) malloc(partnum*sizeof(struct databound));


 cout << "OpenMP : On, threads =" << omp_get_max_threads() << endl;

	for (int idx_par=0; idx_par<partnum; idx_par++){
		tharg[idx_par].p_coord = &coord_cell;
		tharg[idx_par].p_map = &map;
		tharg[idx_par].p_part = &part;
		tharg[idx_par].dims = 2;

		nxtdat[idx_par].start = 0;
		nxtdat[idx_par].end = ncell;
	}

	for (int idx_lev=0; idx_lev<numoflevel; idx_lev++){

		// Start threads
		#pragma omp parallel
		{
			#pragma omp for
			for (int idx_par=0; idx_par<(1<<idx_lev); idx_par++){
				tharg[idx_par].start = nxtdat[idx_par].start;
				tharg[idx_par].end = nxtdat[idx_par].end;
				tharg[idx_par].cur_lvl = numoflevel - idx_lev;

				//cout << "    level=" << idx_lev << " part=" << idx_par << " start=" << nxtdat[idx_par].start << " end=" << nxtdat[idx_par].end << endl;

				rcb_rec2D_ser(&tharg[idx_par]);
			
				/*
				if(pthread_create(&tharg[idx_par].th, NULL, rcb_rec2D_ser, (void*)&tharg[idx_par]) !=0){
					cout << "Error: Can not generate thread at idx_lev=" << idx_lev << " idx_par=" << idx_par << endl;
					exit(-1);
				}
				*/
				//sleep(0.2);
			}
		}

		//Join threads
		/*
		for (int idx_par=0; idx_par<(1<<idx_lev); idx_par++){
			if(pthread_join(tharg[idx_par].th, NULL)!=0)
				cout << "thread join error" << endl;
		}
		*/

		// set next partition start and end point
		for (int idx_par=0; idx_par<(1<<idx_lev); idx_par++){
			if(idx_lev<numoflevel-1){
				nxtdat[2*idx_par].start   = tharg[idx_par].start1;
				nxtdat[2*idx_par].end     = tharg[idx_par].end1;
				nxtdat[2*idx_par+1].start = tharg[idx_par].start2;
				nxtdat[2*idx_par+1].end   = tharg[idx_par].end2;
			}
		}
	}

	free(tharg);
	free(nxtdat);

	// Debug Print
	//part_rcb_DebugPrint(ncell, dims, &coord_cell, &map, &part);

/*	
	if(check_partdata(ncell, &map, &part, numoflevel))
		generate_partdata(nnode, nedge, nbedge, ncell,
                          &map, &part, 
                          partnode, partedge, partbedge, partcell, 
                          cell, ecell, becell);
	else{
		cout << "partition data is invalid" << endl;
		exit(-1);
	}
*/

printf("span  =%lf\n", tmr_span.total_time());
printf("pivot =%lf\n", tmr_pivot1.total_time());
printf("part1 =%lf\n", tmr_part1.total_time());
printf("part2 =%lf\n", tmr_part2.total_time());
printf("part3 =%lf\n", tmr_part3.total_time());
printf("out   =%lf\n", tmr_out.total_time());

}

void part_rcb_writefile(int nnode, int nedge, int nbedge, int ncell,
                        string fname_node, point *partnode,
                        string fname_edge, point *partedge,
                        string fname_bedge, point *partbedge,
                        string fname_cell, point *partcell){

	// Cell data
	cout << "writing cell partitioning data " << fname_cell << "..." << endl;

	FILE *fpcell;
 	if ( (fpcell = fopen(fname_cell.c_str(),"w")) == NULL) {
		printf("can't open file %s\n", fname_cell.c_str()); exit(-1);
	}

	for(int i=0; i<ncell; i++)
		fprintf(fpcell, "%d\n", partcell[i].part);

	fclose(fpcell);

	// Node data
	cout << "writing node partitioning data " << fname_node << "..." << endl;

	FILE *fpnode;
 	if ( (fpnode = fopen(fname_node.c_str(),"w")) == NULL) {
		printf("can't open file %s\n", fname_node.c_str()); exit(-1);
	}

	for(int idx_node=0; idx_node<nnode; idx_node++)
		fprintf(fpnode, "%d\n", partnode[idx_node].part);

	fclose(fpnode);

	// Edge data
	cout << "writing edge partitioning data " << fname_edge << "..." << endl;

	FILE *fpedge;
 	if ( (fpedge = fopen(fname_edge.c_str(),"w")) == NULL) {
		printf("can't open file %s\n", fname_edge.c_str()); exit(-1);
	}

	for(int idx_edge=0; idx_edge<nedge; idx_edge++)
		fprintf(fpedge, "%d\n", partedge[idx_edge].part);

	fclose(fpedge);


	// Edge data
	cout << "writing bedge partitioning data " << fname_bedge << "..." << endl;

	FILE *fpbedge;
 	if ( (fpbedge = fopen(fname_bedge.c_str(),"w")) == NULL) {
		printf("can't open file %s\n", fname_bedge.c_str()); exit(-1);
	}

	for(int idx_bedge=0; idx_bedge<nbedge; idx_bedge++)
		fprintf(fpbedge, "%d\n", partbedge[idx_bedge].part);

	fclose(fpbedge);

/*
	// Generate mapdata file
	FILE *fp2;
 	if ( (fp2 = fopen("map.dat","w")) == NULL) {
		printf("can't open file %s\n", fname.c_str()); exit(-1);
	}
	for(int i=0; i<nnode; i++){
		fprintf(fp2, "%d\n", p_map->at(i));
	}
	fclose(fp2);
*/

}

