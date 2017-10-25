#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>// sqrt function
#include <mpi.h>
#include <omp.h>
#include "read_bmp_clib.h"

FILE*file_bmp;
uint8_t *bmp_data;
#define num_omp_threads 4

//print out the information of serial program
void print_serialinfo(int threshold,double seconds){
        printf("Time taken for serial operation: %f\n", seconds);
        printf("Threshold during convergence: %d\n", threshold);
        printf("\n");
}

//print out the information of MPI parallel version
void print_mpiinfo(int threshold,double seconds){
        printf("Time taken for MPI operation: %f\n", seconds);
        printf("Threshold during convergence: %d\n", threshold);
}

//serial version
void serial_program(char*output_image){
        time_t timer1, timer2;
        time(&timer1);//current time
        //read the binary bmp file into buffer
        uint8_t *new_bmp_img;
        bmp_data=(uint8_t *)read_bmp_file(file_bmp);
        //allocate new output buffer of same size
        new_bmp_img=(uint8_t *)malloc(get_num_pixel());
        //get image attributes
        uint32_t wd=get_image_width();
        uint32_t ht=get_image_height();
        //convergence loop
        int threshold=0;
        int black_cell_count=0;
        int i,j,Gx,Gy,mag;
        while(black_cell_count<(75*wd*ht/100)){
		black_cell_count=0;
		threshold+=1;
		for(i=1;i<(ht-1);i++){
			for(j=1;j<(wd-1);j++){
				Gx=bmp_data[(i-1)*wd+(j+1)]-bmp_data[(i-1)*wd+(j-1)]+2*bmp_data[(i)*wd+(j+1)]-2*bmp_data[(i)*wd+(j-1)]+bmp_data[(i+1)*wd+(j+1)]-bmp_data[(i+1)*wd+(j-1)];
				Gy=bmp_data[(i-1)*wd+(j-1)]+2*bmp_data[(i-1)*wd+(j)]+bmp_data[(i-1)*wd+(j+1)]-bmp_data[(i+1)*wd+(j-1)]-2*bmp_data[(i+1)*wd+(j)]-bmp_data[(i+1)*wd+(j+1)];
				mag=sqrt(Gx*Gx+Gy*Gy);
				if(mag>threshold)
				{
					new_bmp_img[i*wd+j]=255;
				}
				else{
					new_bmp_img[i*wd+j]=0;
					black_cell_count++;
				}
			}
		}
	}
	//copy the boundary
	for(j=0;j<wd;j++){
		new_bmp_img[j]=bmp_data[j];//first row
		new_bmp_img[(ht-1)*wd+j]=bmp_data[(ht-1)*wd+j];//last row
	}
	for(i=0;i<ht;i++){
		new_bmp_img[i*wd]=bmp_data[i*wd];//first column
		new_bmp_img[i*wd+wd-1]=bmp_data[i*wd+wd-1];//last column
	}
	FILE*file_bmp2;
	file_bmp2=fopen(output_image,"wb");
	write_bmp_file(file_bmp2,new_bmp_img);
	time(&timer2);//current time
	double seconds=difftime(timer2,timer1);//time
	print_serialinfo(threshold,seconds);
}

void Sobel_operation(int index,int start_iteration,int end_iteration,uint32_t wd,int *place_blackcells){
	int black_cell_count=0;
	int end_index=0;
	int num_box= (end_iteration-start_iteration+1)*(wd-2);
	int each_iteration=num_box/num_omp_threads;
	int start_index=index*each_iteration;
	int start_i=start_index/(wd-2)+start_iteration;
	int start_j=start_index%(wd-2)+1;
	if(index!=num_omp_threads-1)//not last thread
	{
		end_index=start_index+each_iteration-1;
	}
	else //last thread
	{
		end_index=start_index+each_iteration+num_box%num_omp_threads-1;
	}
	int end_i=end_index/(wd-2)+start_iteration;
	int end_j=end_index%(wd-2)+1;
	int i,j,Gx,Gy,mag;
	for(i=start_i;i<=end_i;i++){
		for(j=start_j;j<=end_j;j++){
			Gx=bmp_data[(i-1)*wd+(j+1)]-bmp_data[(i-1)*wd+(j-1)]+2*bmp_data[(i)*wd+(j+1)]-2*bmp_data[(i)*wd+(j-1)]+bmp_data[(i+1)*wd+(j+1)]-bmp_data[(i+1)*wd+(j-1)];
			Gy=bmp_data[(i-1)*wd+(j-1)]+2*bmp_data[(i-1)*wd+(j)]+bmp_data[(i-1)*wd+(j+1)]-bmp_data[(i+1)*wd+(j-1)]-2*bmp_data[(i+1)*wd+(j)]-bmp_data[(i+1)*wd+(j+1)];
			mag=sqrt(Gx*Gx+Gy*Gy);
			if(mag>threshold)
			{
				//new_bmp_img[i*wd+j]=255;
			}
			else{
				//new_bmp_img[i*wd+j]=0;
				place_blackcells[index]++;
			}
		}
	}
}

int main(int argc, char* argv[]){

	file_bmp=fopen(argv[1],"rb");
	//serial program
	serial_program(argv[2]);

	//mpi program
	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	//master process
	if(rank==0){
		//initialize
		int *place_blackcells;//place temporary number of black cells
		place_blackcells=malloc(num_omp_threads*sizeof(int));
		int mm=0;
		for(mm=0;mm<num_omp_threads;mm++){
			place_blackcells[mm]=0;
		}
		//read the binary bmp file into buffer
        uint8_t *new_bmp_img;
        bmp_data=(uint8_t *)read_bmp_file(file_bmp);
        //allocate new output buffer of same size
        new_bmp_img=(uint8_t *)malloc(get_num_pixel());
        //get image attributes
        uint32_t wd=get_image_width();
        uint32_t ht=get_image_height();

		MPI_Barrier(MPI_COMM_WORLD);
		time_t timer1, timer2;
		time(&timer1);//current time

		//convergence loop
		int iteration_lines=(ht-2)/size+(ht-2)%size; //number of lines needed to execute
		int start_iteration=1;
		int end_iteration=1+iteration_lines-1;
        int threshold=0;
        int black_cell_count=0;
        int global_black_cell=0;
        while(global_black_cell<(75*wd*ht/100)){
			black_cell_count=0;
			threshold+=1;
			#pragma omp parallel num_threads(num_omp_threads)
			{
				Sobel_operation(omp_get_thread_num(),start_iteration,end_iteration,wd,place_blackcells);
			}
			for(mm=0;mm<num_omp_threads;mm++)
				black_cell_count=black_cell_count+place_blackcells[mm];
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Reduce(&black_cell_count,&global_black_cell,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		}

		//gather image pixels from all slaves

		time(&timer2);//current time
		double seconds=difftime(timer2,timer1);//time
		//copy the boundary
		for(j=0;j<wd;j++){
			new_bmp_img[j]=bmp_data[j];//first row
			new_bmp_img[(ht-1)*wd+j]=bmp_data[(ht-1)*wd+j];//last row
		}
		for(i=0;i<ht;i++){
			new_bmp_img[i*wd]=bmp_data[i*wd];//first column
			new_bmp_img[i*wd+wd-1]=bmp_data[i*wd+wd-1];//last column
		}
		FILE*file_bmp2;
		file_bmp2=fopen(output_image,"wb");
		write_bmp_file(file_bmp2,new_bmp_img);
		print_mpiinfo(threshold,seconds);
	}
	else //slave process
	{
		//initialize
		int *place_blackcells;//place temporary number of black cells
		place_blackcells=malloc(num_omp_threads*sizeof(int));
		int mm=0;
		for(mm=0;mm<num_omp_threads;mm++){
			place_blackcells[mm]=0;
		}
		//read the binary bmp file into buffer
        uint8_t *new_bmp_img;
        bmp_data=(uint8_t *)read_bmp_file(file_bmp);
        //allocate new output buffer of same size
        new_bmp_img=(uint8_t *)malloc(get_num_pixel());
        //get image attributes
        uint32_t wd=get_image_width();
        uint32_t ht=get_image_height();

		MPI_Barrier(MPI_COMM_WORLD);
		//convergence loop
		int iteration_lines=(ht-2)/size; //number of lines needed to execute
		int start_iteration=rank*((ht-2)/size)+(ht-2)%size+1;
		int end_iteration=start_iteration+iteration_lines-1;
        int threshold=0;
        int black_cell_count=0;
        int global_black_cell=0;
        while(global_black_cell<(75*wd*ht/100)){
			black_cell_count=0;
			threshold+=1;
			#pragma omp parallel num_threads(num_omp_threads)
			{
				Sobel_operation(omp_get_thread_num(),start_iteration,end_iteration,wd,place_blackcells);
			}
			for(mm=0;mm<num_omp_threads;mm++)
				black_cell_count=black_cell_count+place_blackcells[mm];
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Reduce(&black_cell_count,&global_black_cell,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		}

		//send image pixels to master

	}
	fclose(file_bmp);
	MPI_Finalize();
	return(0);
}
