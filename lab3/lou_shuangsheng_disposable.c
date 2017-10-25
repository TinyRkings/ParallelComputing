#include <stdlib.h>
#include <stdio.h>
#include <string.h>//memset
#include <time.h>
#include <omp.h>
#include "type.h"
#define each_box 7
#define NS_PER_US 1000

struct box* boxes;
static float affect_rate;
static float epsilon;
int num_box;//number of box
int num_threads;//number of threads
float*place_temp;//place temporary value of each box

//print the information of the number of threads actually created
void print_actual_thread(int actual_thread){
	printf("The requested number of OpenMP threads: %d\n",num_threads);
	printf("The actual number of threads created: %d\n",actual_thread);
}

//print the information of convergence and runing time
void print_info(float max_dsv, float min_dsv, int iteration,double diff,clock_t consum,double seconds){
	printf("\n****************************************************************\n");
	printf("dissipation converged in %d iterations,\n", iteration);
	printf("    with max DSV = %f and min DSV = %f\n", max_dsv, min_dsv);
	printf("    affect rate  = %f;    epsilon = %f\n", affect_rate, epsilon);
	printf("****************************************************************\n");
	printf("\n");
	printf("elapsed convergence loop time  (clock): %d\n", (int) consum);
	printf("elapsed convergence loop time  (time): %f\n", seconds);
	printf("elapsed convergence loop time  (realtime): %f\n", diff);
}

//calculate the average adjacent temperature
float cal_adjacent_tem(int index){
	//top
	int top_contact=0;
	float top_tem=0.0;
	int i=0;
	int left_x=boxes[index].location[1];
	int width=boxes[index].location[3];
	int temp_contact=0;
	if(boxes[index].topnum_neigh==0)
	{
		top_contact=width;
		top_tem=width*boxes[index].temperature;
	}
	else{
	struct node*first=boxes[index].toplist;
	for(i=0;i<boxes[index].topnum_neigh;i++){
		int m=first->id;
		if((left_x+width)>=(boxes[m].location[1]+boxes[m].location[3]))
		{
			if(left_x>=boxes[m].location[1])
				temp_contact=boxes[m].location[1]+boxes[m].location[3]-left_x;
			else
				temp_contact=boxes[m].location[3];
		}
		else
		{
			if(left_x>=boxes[m].location[1])
				temp_contact=width;
			else
				temp_contact=left_x+width-boxes[m].location[1];
		}
		top_contact=top_contact+temp_contact;
		top_tem=top_tem+temp_contact*boxes[m].temperature;
		first=first->next;
	}
	}//else
	//bottom
	int bottom_contact=0;
	float bottom_tem=0.0;
	if(boxes[index].bottomnum_neigh==0)
	{
		bottom_contact=width;
		bottom_tem=width*boxes[index].temperature;
	}
	else{
	struct node*first=boxes[index].bottomlist;
	for(i=0;i<boxes[index].bottomnum_neigh;i++){
		int m=first->id;
		if((left_x+width)>=(boxes[m].location[1]+boxes[m].location[3]))
		{
			if(left_x>=boxes[m].location[1])
			{	temp_contact=boxes[m].location[1]+boxes[m].location[3]-left_x;
			}
			else
			{	temp_contact=boxes[m].location[3];
			}
		}
		else
		{
			if(left_x>=boxes[m].location[1])
			{	temp_contact=width;
			}
			else
			{	temp_contact=left_x+width-boxes[m].location[1];
			}
		}
		bottom_contact=bottom_contact+temp_contact;
		bottom_tem=bottom_tem+temp_contact*boxes[m].temperature;
		first=first->next;
	}
	}//else
	//left
	int left_y=boxes[index].location[0];
	int height=boxes[index].location[2];
	int left_contact=0;
	float left_tem=0.0;
	if(boxes[index].leftnum_neigh==0)
	{
		left_contact=height;
		left_tem=height*boxes[index].temperature;
	}
	else{
	struct node*first=boxes[index].leftlist;
	for(i=0;i<boxes[index].leftnum_neigh;i++){
		int m=first->id;
		if((left_y+height)>=(boxes[m].location[0]+boxes[m].location[2]))
		{
			if(left_y>=boxes[m].location[0])
				temp_contact=boxes[m].location[0]+boxes[m].location[2]-left_y;
			else
				temp_contact=boxes[m].location[2];
		}
		else
		{
			if(left_y>=boxes[m].location[0])
				temp_contact=height;
			else
				temp_contact=left_y+height-boxes[m].location[0];
		}
		left_contact=left_contact+temp_contact;
		left_tem=left_tem+temp_contact*boxes[m].temperature;
		first=first->next;
	}
	}//else
	//right
	int right_contact=0;
	float right_tem=0.0;
	if(boxes[index].rightnum_neigh==0)
	{
		right_contact=height;
		right_tem=height*boxes[index].temperature;
	}
	else{
	struct node*first=boxes[index].rightlist;
	for(i=0;i<boxes[index].rightnum_neigh;i++){
		int m=first->id;
		if((left_y+height)>=(boxes[m].location[0]+boxes[m].location[2]))
		{
			if(left_y>=boxes[m].location[0])
				temp_contact=boxes[m].location[0]+boxes[m].location[2]-left_y;
			else
				temp_contact=boxes[m].location[2];
		}
		else
		{
			if(left_y>=boxes[m].location[0])
				temp_contact=height;
			else
				temp_contact=left_y+height-boxes[m].location[0];
		}
		right_contact=right_contact+temp_contact;
		right_tem=right_tem+temp_contact*boxes[m].temperature;
		first=first->next;
	}
	}//else
	int total_contact=top_contact+bottom_contact+left_contact+right_contact;
	float total_tem=top_tem+bottom_tem+left_tem+right_tem;
	float average_tem=total_tem/total_contact;
	return average_tem;
}


//calculate the max temperature of all boxes
float max_temp(){
	float max_t=boxes[0].temperature;
	int i=1;
	for(i=1;i<num_box;i++)
	{
		if(boxes[i].temperature>max_t)
			max_t=boxes[i].temperature;
	}
	return max_t;
}

//calculate the min temperature of all boxes
float min_temp(){
	float min_t=boxes[0].temperature;
	int i=1;
	for(i=1;i<num_box;i++)
	{
		if(boxes[i].temperature<min_t)
			min_t=boxes[i].temperature;
	}
	return min_t;
}

//calculate the new temperature of current boxes
void* disspate_box(int tn){
	int index=tn;
	int temp=0;
	if(index!=num_threads-1)//not last thread
	{
		int each_iteration=num_box/num_threads;
		int end_iteration=index*each_iteration+each_iteration;
		for(temp=index*each_iteration;temp<end_iteration;temp++)
		{
			float average_tem=cal_adjacent_tem(temp);
			float new_tem;
			if(boxes[temp].temperature<=average_tem)
			{
				new_tem=boxes[temp].temperature+(average_tem-boxes[temp].temperature)*affect_rate;
			}
			else
			{
				new_tem=boxes[temp].temperature-(boxes[temp].temperature-average_tem)*affect_rate;
			}
			place_temp[temp]=new_tem;
		}
	}
	else
	{
		int each_iteration=num_box/num_threads;
		for(temp=index*each_iteration;temp<num_box;temp++)
		{
			float average_tem=cal_adjacent_tem(temp);
			float new_tem;
			if(boxes[temp].temperature<=average_tem)
			{
				new_tem=boxes[temp].temperature+(average_tem-boxes[temp].temperature)*affect_rate;
			}
			else
			{
				new_tem=boxes[temp].temperature-(boxes[temp].temperature-average_tem)*affect_rate;
			}
			place_temp[temp]=new_tem;
		}
	}
}

//parallel the computations of temperature of each box and output the result
void compute_tem(){
	time_t timer1, timer2;
	time(&timer1);//current time
	clock_t clock_consum;
	clock_consum=clock();//clock consumption
	struct timespec start, end;
	double diff;
	clock_gettime(CLOCK_REALTIME, &start);
	float min_tem=min_temp();
	float max_tem=max_temp();
	int iteration=0;
	place_temp=malloc(sizeof(float)*num_box);
	int qq=0;
	for(qq=0;qq<num_box;qq++)//initial
	{
		place_temp[qq]=0.0;
	}
	int actual_thread=num_threads;
	while((max_tem-min_tem)>=max_tem*epsilon){
		#pragma omp parallel num_threads(num_threads)
		{
			#pragma omp single
			{
				if(actual_thread!=omp_get_num_threads())
					actual_thread=omp_get_num_threads();
			}
			disspate_box(omp_get_thread_num());
		}
		for(qq=0;qq<num_box;qq++)
			boxes[qq].temperature=place_temp[qq];
		min_tem=min_temp();
		max_tem=max_temp();
		iteration++;
	}
	free(place_temp);
	time(&timer2);//current time
	double seconds=difftime(timer2,timer1);//time
	clock_consum=clock()-clock_consum;//clock consumption
	clock_gettime(CLOCK_REALTIME, &end);
	//diff=(double)(((end.tv_sec-start.tv_sec)*CLOCKS_PER_SEC)+((end.tv_nsec-start.tv_nsec)/NS_PER_US));
	diff=(double)(((end.tv_sec-start.tv_sec))+((end.tv_nsec-start.tv_nsec)/1000000000));
	print_actual_thread(actual_thread);
	print_info(max_tem,min_tem,iteration,diff,clock_consum,seconds);
}

//read the data and recall compute_tem() function
int main(int argc, char* argv[]){
	if(argc!=4)
	{	
		printf("%d wrong input!",argc);
		return -1;
	}
  affect_rate=atof(argv[1]);
  epsilon=atof(argv[2]);
  num_threads=atoi(argv[3]);
  scanf("%d",&num_box);
  int rows;
  int cols;
  scanf("%d",&rows);
  scanf("%d",&cols);
  boxes=malloc(sizeof(struct box)*num_box);
  int j=0;
  for(j=0;j<num_box;j++){
    int i=0;
    for(i=0;i<each_box;){
	char temp;
	if(i==0){  //id
		scanf("%d",&(boxes[j].id));
	}
	else if(i==1){ //location
		int m=0;
		for(m=0;m<4;m++)
		{	scanf("%d",&(boxes[j].location[m]));
		}
	}
	else if(i==2){  //topnum_neigh + toplist
		scanf("%d",&(boxes[j].topnum_neigh));
		int ww=0;
		struct node* current=NULL;
		if(boxes[j].topnum_neigh!=0)
		{
			boxes[j].toplist=malloc(sizeof(struct node));
			current=boxes[j].toplist;
		}
		else
			boxes[j].toplist=NULL;
		for(ww=0;ww<boxes[j].topnum_neigh;ww++){
			if(ww!=0){
				current->next=malloc(sizeof(struct node));
				current=current->next;
				scanf("%d",&(current->id));
				current->next=NULL;
			}
			else{
				scanf("%d",&(current->id));
				current->next=NULL;
			}
		}
	}
	else if(i==3){  //bottomnum_neigh + bottomlist
		scanf("%d",&(boxes[j].bottomnum_neigh));
		int ww=0;
		struct node* current=NULL;
		if(boxes[j].bottomnum_neigh!=0)
		{
			boxes[j].bottomlist=malloc(sizeof(struct node));
			current=boxes[j].bottomlist;
		}
		else
			boxes[j].bottomlist=NULL;
		for(ww=0;ww<boxes[j].bottomnum_neigh;ww++){
			if(ww!=0){
				current->next=malloc(sizeof(struct node));
				current=current->next;
				scanf("%d",&(current->id));
				current->next=NULL;
			}
			else{
				scanf("%d",&(current->id));
				current->next=NULL;
			}
		}
	}
	else if(i==4){  //leftnum_neigh + leftlist
		scanf("%d",&(boxes[j].leftnum_neigh));
		int ww=0;
		struct node* current=NULL;
		if(boxes[j].leftnum_neigh!=0)
		{
			boxes[j].leftlist=malloc(sizeof(struct node));
			current=boxes[j].leftlist;
		}
		else
			boxes[j].leftlist=NULL;
		for(ww=0;ww<boxes[j].leftnum_neigh;ww++){
			if(ww!=0){
				current->next=malloc(sizeof(struct node));
				current=current->next;
				scanf("%d",&(current->id));
				current->next=NULL;
			}
			else{
				scanf("%d",&(current->id));
				current->next=NULL;
			}
		}
	}
	else if(i==5){  //rightnum_neigh + rightlist
		scanf("%d",&(boxes[j].rightnum_neigh));
		int ww=0;
		struct node* current=NULL;
		if(boxes[j].rightnum_neigh!=0)
		{
			boxes[j].rightlist=malloc(sizeof(struct node));
			current=boxes[j].rightlist;
		}
		else
			boxes[j].rightlist=NULL;
		for(ww=0;ww<boxes[j].rightnum_neigh;ww++){
			if(ww!=0){
				current->next=malloc(sizeof(struct node));
				current=current->next;
				scanf("%d",&(current->id));
				current->next=NULL;
			}
			else{
				scanf("%d",&(current->id));
				current->next=NULL;
			}
		}
	}
	else if(i==6){
		scanf("%f",&(boxes[j].temperature));
	}
	i++;	
     }   //i loop
    }   //j loop
	compute_tem();
	int www=0;
	for(www=0;www<num_box;www++){
		free(boxes[www].toplist);
		free(boxes[www].bottomlist);
		free(boxes[www].leftlist);
		free(boxes[www].rightlist);
	}
	free(boxes);
	return 0;
}
