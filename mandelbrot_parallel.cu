#include <stdio.h>
#include <stdlib.h>
#include "bmpfile.h"

/*Mandelbrot values*/
#define RESOLUTION 8700.0
#define XCENTER -0.55
#define YCENTER 0.4
#define MAX_ITER 1000/2

/*Colour Values*/
#define COLOUR_DEPTH 255
#define COLOUR_MAX 100.0
#define GRADIENT_COLOUR_MAX 200.0

#define FILENAME "my_mandelbrot_fractal3.bmp"
#define BLOCKSIZE 16

/**
  * Computes the color gradiant
  * color: the output vector
  * x: the gradiant (beetween 0 and 360)
  * min and max: variation of the RGB channels (Move3D 0 -> 1)
  * Check wiki for more details on the colour science: en.wikipedia.org/wiki/HSL_and_HSV
  */
void GroundColorMix(double* color, double x, double min, double max)
{
	/*
	 * Red = 0
	 * Green = 1
	 * Blue = 2
	 */
	double posSlope = (max-min)/60;
	double negSlope = (min-max)/60;

	if( x < 60 )
	{
		color[0] = max;
		color[1] = posSlope*x+min;
		color[2] = min;
		return;
	}
	else if ( x < 120 )
	{
		color[0] = negSlope*x+2.0*max+min;
		color[1] = max;
		color[2] = min;
		return;
	}
	else if ( x < 180  )
	{
		color[0] = min;
		color[1] = max;
		color[2] = posSlope*x-2.0*max+min;
		return;
	}
	else if ( x < 240  )
	{
		color[0] = min;
		color[1] = negSlope*x+4.0*max+min;
		color[2] = max;
		return;
	}
	else if ( x < 300  )
	{
		color[0] = posSlope*x-4.0*max+min;
		color[1] = min;
		color[2] = max;
		return;
	}
	else
	{
		color[0] = max;
		color[1] = min;
		color[2] = negSlope*x+6*max;
		return;
	}
}

/* Mandelbrot Set Image Demonstration
 *
 * This is a simple single-process/single thread implementation
 * that computes a mandelbrot set and produces a corresponding
 * Bitmap image. The program demonstrates the use of a colour
 * gradient
 *
 * This program uses the algorithm outlined in:
 *   "Building Parallel Programs: SMPs, Clusters And Java", Alan Kaminsky
 *
 * This program requires libbmp for all bitmap operations.
 *
 */

/*
 * Function to read in height and width values for the bmp
 * file to be produce. Exits if less than 3 arguments given.
 */
void getParameters(int argc, char** argv, int* height, int* width, size_t* size){
	if (argc < 3 || argv[1] < 0 || argv[2] < 0){
		printf("Usage: <mandelbrot_parallel> <height> <width> \n\n");
		exit(EXIT_SUCCESS);
	}

	*height = atoi(argv[1]);
	*width = atoi(argv[2]);
	*size = *height * *width;
};

/*
 * Structure to hold metadata needed for the production of a series
 * of pixels for Mandelbrot fractal.
 */
typedef struct {
	int width;
	int height;
	float xcenter;
	float ycenter;
	float resolution;
	int iterations;
} Mandelbrot;

/*
 * Kernel functions which returns an array to output of 'iter' values
 * to be used in color function.
 */
/*
	 * 	KERNEL FUNCTION TO POPULATE VALUES
	 */
//col and row values to be worked out from thread position?
//double x = XCENTER + (xoffset + col) / RESOLUTION;
//double y = YCENTER + (yoffset - row) / RESOLUTION;

/*
 * We are going to want a thread to do each of the pixels. Therefore
 * we're going to need x*y threads giving us (x*y)/1024 blocks.
 *
 * If we structure the blocks in a 2d grid, we can just have them
 * pluck out their x,y pixel based on their own x,y location.
 *
 * They can return to x+y...this causes a clash...
 *
 * To do the calculation each thread will need an array of
 * x values and an array of y values.
 *
 * To calculate x we need XCENTER, RESOLUTION, and xoffset.
 * To calculate y we need YCENTER, RESOLUTION, and yoffset.
 *
 * We want to return an iter value corresponding to that pixel that
 * it represents.
 */

__global__ void MandelbrotFractal(float* output, Mandelbrot M)
{

	//get information from 2D block/thread grid
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	//we're only interested in processing threads that fall within
	//the boundaries of the picure
	if (row < M.width && col < M.height){
		int xoffset = -(M.width - 1) /2;
		int yoffset = (M.height -1) / 2;

		//Determine where in the mandelbrot set, the pixel is referencing
		double x = M.xcenter + (xoffset + col) / M.resolution;
		double y = M.ycenter + (yoffset - row) / M.resolution;

		//Mandelbrot stuff
		double a = 0;
		double b = 0;
		double aold = 0;
		double bold = 0;
		double zmagsqr = 0;
		int iter = 0; //import one!

		//Check if the x,y coord are part of the mendelbrot set - refer to the algorithm
		while(iter < M.iterations && zmagsqr <= 4.0){
			++iter;
			a = (aold * aold) - (bold * bold) + x;
			b = 2.0 * aold*bold + y;

			zmagsqr = a*a + b*b;

			aold = a;
			bold = b;
		}
		//output is a 1D array, so we need to index using our row and
		//column number
		output[row * M.width + col] = iter;

	}
}

/*
 * Function to package globals for easier sending to
 * device.
 */
void makeMandel(MandelBrot* M){
	M.iterations = MAX_ITER;
	M.resolution = RESOLUTION;
	M.ycenter = YCENTER;
	M.xcenter = XCENTER;
}

int main(int argc, char **argv)
{
	int height, width;
	size_t size;

	cudaError_t error;

	getParameters(argc, argv, &height, &width, &size);

	bmpfile_t *bmp;
	rgb_pixel_t pixel = {0, 0, 0, 0};
	int xoffset = -(width - 1) /2;
	int yoffset = (height -1) / 2;
	bmp = bmp_create(width, height, 32);

	Mandelbrot h_mandel;
	makeMandel(&h_mandel);
	h_mandel.width = width;
	h_mandel.height = height;


	//memory to hold results
	float* h_xy = ( float*) malloc (size);

	//allocate device memory
	float* d_xy;

	error = cudaMalloc(&d_xy, size);

	if (error != cudaSuccess)
	{
		printf("cudaMalloc d_xy returned error %s (code %d), line(%d)\n", cudaGetErrorString(error), error, __LINE__);
		exit(EXIT_FAILURE);
	}

	Mandelbrot* d_Mandel;

	error = cudaMalloc(&d_Mandel, sizeof(Mandelbrot));
	if (error != cudaSuccess)
	{
		printf("cudaMalloc d_Mandel returned error %s (code %d), line(%d)\n", cudaGetErrorString(error), error, __LINE__);
		exit(EXIT_FAILURE);
	}

	//copy Mandelbrot metadata to device
	error = cudaMemcpy(d_Mandel, h_mandel, sizeof(Mandelbrot), cudaMemcpyHostToDevice);

	if (error != cudaSuccess)
	{
		printf("cudaMemcpy (d_Mandel,h_mandel) returned error %s (code %d), line(%d)\n", cudaGetErrorString(error), error, __LINE__);
		exit(EXIT_FAILURE);
	}


	//figure out blocks
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

	//figure out threads
	dim3 dimGrid(h_mandel.width / dimBlock.x + 1, h_mandel.height / dimBlock.y + 1);

	//call kernel function
	MandelbrotFractal<<<dimGrid, dimBlock>>>(d_xy, d_Mandel);

	//get data from kernel to device
	error = cudaMemcpy(h_xy, d_xy, size, cudaMemcpyDeviceToHost);
	if (error != cudaSuccess)
	{
		printf("cudaMemcpy (h_xy,d_xy) returned error %s (code %d), line(%d)\n", cudaGetErrorString(error), error, __LINE__);
		exit(EXIT_FAILURE);
	}


	/* Generate the colour of the pixel from the **iter** value */
	/* You can mess around with the colour settings to use different gradients */
	/* Colour currently maps from royal blue to red */
	/* We're interested in iter */
	int i;
	for (i = 0; i < width*height; i++){
		x_col =  (COLOUR_MAX - (( ((float) iter / ((float) MAX_ITER) * GRADIENT_COLOUR_MAX))));
		GroundColorMix(color, x_col, 1, COLOUR_DEPTH);
		pixel.red = color[0];
		pixel.green = color[1];
		pixel.blue = color[2];

		int row = i/width;
		int col = i % width;

		//adds pixel color to image
		bmp_set_pixel(bmp, col, row, pixel);
	}


	bmp_save(bmp, FILENAME);

	//free all memory used
	bmp_destroy(bmp);

	free(h_xy);
	free(h_mandel);

	cudaFree(d_xy);
	cudaFree(d_Mandel);

	return 0;
}
