#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPNGWriter.h>
#include <vtkImageData.h>
#include <cstdio>
#include <cstdlib>

#include <sys/time.h>

#define WIDTH 500
#define HEIGHT 500
#define SAMPLE_SIZE 100

int
GetPointIndex(const int x, const int y, const int z, const int *dims)
{
    return z*dims[0]*dims[1] + y*dims[0] + x;
}

double 
EvaluateFieldAtLocation(const int *dims, const float *X, const float *Y, const float *Z, const float *F, const double *pt)
{
    int x, y, z;
    x = y = z = -1;

    /*first check to see if out of bounds*/
    for(int i = 0; i < dims[0]; i++)
    {
        if(pt[0] <= X[i])
        {
            x = i-1;
            for(int j=0; j < dims[1]; j++)
            {
                if(pt[1] <= Y[j])
                {
                    y = j-1;
                    for(int k = 0; k < dims[2]; k++)
                    {
                        if(pt[2] <= Z[k])
                        {
                            z = k-1;
                            break;
                        }
                    }
                    break;
                }
            }
            break;
        }
    }
    
    if (x == -1 || y == -1 || z == -1)
        return 0.0;

    /*not out of bounds, interpolation*/
    double x_prop = (pt[0] - X[x]) / (X[x+1] - X[x]);
    double x_0 = F[GetPointIndex(x,y,z, dims)]*(1-x_prop)
        +F[GetPointIndex(x+1,y,z, dims)]*x_prop;
    double x_1 = F[GetPointIndex(x,y,z+1, dims)]*(1-x_prop)
        +F[GetPointIndex(x+1,y,z+1, dims)]*x_prop;
    double x_2 = F[GetPointIndex(x,y+1,z, dims)]*(1-x_prop)
        +F[GetPointIndex(x+1,y+1,z, dims)]*x_prop;
    double x_3 = F[GetPointIndex(x,y+1,z+1, dims)]*(1-x_prop)
        +F[GetPointIndex(x+1,y+1,z+1, dims)]*x_prop;
    
    double y_prop = (pt[1] - Y[y]) / (Y[y+1] - Y[y]);
    double y_0 = x_0*(1-y_prop)+x_2*y_prop;
    double y_1 = x_1*(1-y_prop)+x_3*y_prop;

    double z_prop = (pt[2] - Z[z]) / (Z[z+1] - Z[z]);
    double result = y_0*(1-z_prop)+y_1*z_prop;
    return result;
}

struct Camera
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};


struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    double         *opacities; // size is numBins

    int GetBin(double value){
    	double range_of_bins = (max-min)/numBins;
    	return (int)((value-min)/range_of_bins);
    }
    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
    {
        int bin = GetBin(value);
        RGB[0] = colors[3*bin+0];
        RGB[1] = colors[3*bin+1];
        RGB[2] = colors[3*bin+2];
        opacity = opacities[bin];
    }
};

TransferFunction
SetupTransferFunction(void)
{
    int  i;

    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    rv.opacities = new double[256];
    unsigned char charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = charOpacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = { 
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76 
       };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0 ; i < numControlPoints-1 ; i++)
    {
        int start = controlPointPositions[i]*rv.numBins;
        int end   = controlPointPositions[i+1]*rv.numBins+1;
		// cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << endl;
        if (end >= rv.numBins)
            end = rv.numBins-1;
        for (int j = start ; j <= end ; j++)
        {
            double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (int k = 0 ; k < 3 ; k++)
                rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                                 + controlPointColors[3*i+k];
        }
    }    

    return rv;
}



Camera
SetupCamera(void)
{
    Camera rv;
    rv.focus[0] = 0;
    rv.focus[1] = 0;
    rv.focus[2] = 0;
    rv.up[0] = 0;
    rv.up[1] = -1;
    rv.up[2] = 0;
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;

    return rv;
}

/*takes the cross product of vector A and B and outputs it in vector P*/
void 
crossProduct(const double *vect_A, const double *vect_B, double *cross_P)
{
 
    cross_P[0] = (vect_A[1] * vect_B[2]) - (vect_A[2] * vect_B[1]);
    cross_P[1] = (vect_A[2] * vect_B[0]) - (vect_A[0] * vect_B[2]);
    cross_P[2] = (vect_A[0] * vect_B[1]) - (vect_A[1] * vect_B[0]);
}

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return image;
}


int main()
{
	int i, j, k, x, y, offset;
	Camera camera = SetupCamera();
	TransferFunction transferFunction = SetupTransferFunction(); 

	vtkDataSetReader *rdr = vtkDataSetReader::New();
	rdr->SetFileName("astro64.vtk");
	rdr->Update();

	int dims[3];
	vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
	rgrid->GetDimensions(dims);

	float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
	float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
	float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
	float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

	vtkImageData *output = NewImage(WIDTH, HEIGHT);
	unsigned char *buffer = (unsigned char *) output->GetScalarPointer(0,0,0);

	for(i=0; i< 3*WIDTH*HEIGHT; i++){
		buffer[i]=0;
	}

	/***********************************
	RAY-CASTING
	For every pixel on the screen:
	1)	Find ray
	2)	Intersect volume with that ray
	3)	Calculate color from intersection
	4)	Assign color to that pixel
	************************************/
	double look[3], u[3], v[3], Dx[3], Dy[3], euclidean_distance;

	for(i=0; i<3; i++)
		look[i] = camera.focus[i]-camera.position[i];

	euclidean_distance = sqrt((look[0]*look[0]) + (look[1]*look[1]) + (look[2]*look[2]));
	for(i=0; i<3; i++)
		look[i] = look[i]/euclidean_distance;

	crossProduct( look, camera.up, u);
	euclidean_distance = sqrt((u[0]*u[0])+ (u[1]*u[1])+ (u[2]*u[2]));
	for(i=0; i<3; i++)
		u[i] = u[i]/euclidean_distance;

	crossProduct(look, u, v);
	euclidean_distance = sqrt((v[0]*v[0])+ (v[1]*v[1])+ (v[2]*v[2]));
	for(i=0; i<3; i++)
		v[i] = v[i]/euclidean_distance;

	for(i=0; i<3; i++){
		Dx[i] = (2* tan((camera.angle/2)*(M_PI/180))* u[i])/WIDTH;
		Dy[i] = (2* tan((camera.angle/2)*(M_PI/180))* v[i])/HEIGHT;
	}
	
	/*
	storing the RGBs and opacities in two different arrays for compositing
	putting pixels RGB in a double array instead of unsigned char for precision in color calculations
	*/
	unsigned char RGB[3];
	double s1_opacity, s2_opacity, pixels_RGB[3];
	

	/*****************************/
	/*The following is for timing*/

	struct timeval startTime;
	gettimeofday(&startTime, 0);
	/*****************************/

	for(x=0; x<WIDTH; x++){
		for(y=0; y<HEIGHT; y++){
			/****************************
			I) FINDING RAY FOR EACH PIXEL
			****************************/
			double ray[3], sample_position[3], eval_result;

			for(i=0; i<3; i++)
				ray[i] = look[i] + ((((2*x)+1-WIDTH)/2.0)*Dx[i]) + ((((2*y)+1-HEIGHT)/2.0)*Dy[i]);
			
			/*normalize the ray*/
			euclidean_distance = sqrt((ray[0]*ray[0])+ (ray[1]*ray[1])+ (ray[2]*ray[2]));
			for(i=0; i<3; i++)
				ray[i] = ray[i]/euclidean_distance;

			/***************************************
			II) INTERSECTING THE VOLUME WITH THE RAY
			***************************************/
			for(i=0; i<3; i++)
				sample_position[i]= camera.position[i] + (camera.near * ray[i]);
			double step_size = (camera.far - camera.near)/(SAMPLE_SIZE);

			eval_result = EvaluateFieldAtLocation(dims, X, Y, Z, F, sample_position);
			transferFunction.ApplyTransferFunction(eval_result, RGB, s1_opacity);

			/*correcting the opacity*/
			s1_opacity = 1.0 - pow((1.0-s1_opacity), (500.0/SAMPLE_SIZE));

			/*********************
			III) CALCULATING COLOR
			*********************/
			for(i=0; i<3; i++)
				pixels_RGB[i] = (RGB[i]/255.0) * s1_opacity;

			for(i=1; i<(int)SAMPLE_SIZE; i++){
				for(j=0; j<3; j++){
					sample_position[j] += (step_size*ray[j]);
				}
				eval_result = EvaluateFieldAtLocation(dims, X, Y, Z, F, sample_position);

				if(transferFunction.GetBin(eval_result) >= 0){

					/*getting the color and opacity*/
					transferFunction.ApplyTransferFunction(eval_result, RGB, s2_opacity);

					/*correcting the opacity*/
					s2_opacity = 1.0 - pow((1.0-s2_opacity), (500.0/SAMPLE_SIZE));

					/*compositing*/
					for(j=0; j<3; j++)
						pixels_RGB[j] = pixels_RGB[j] + (1-s1_opacity)*s2_opacity*(RGB[j]/255.0);
					s1_opacity = s1_opacity + (1-s1_opacity)*s2_opacity;

					
				}
			}

			/************************
			IV) ASSIGN COLOR TO PIXEL
			************************/
			offset = 3*(y*WIDTH+x);
			buffer[offset] = (unsigned char)(pixels_RGB[0]*255.0);
			buffer[offset+1] = (unsigned char)(pixels_RGB[1]*255.0);
			buffer[offset+2] = (unsigned char)(pixels_RGB[2]*255.0);
		}
	}

	/*****************************/
	/*The following is for timing*/

	struct timeval endTime;
	gettimeofday(&endTime, 0);
	double seconds = double(endTime.tv_sec - startTime.tv_sec) + double(endTime.tv_usec - startTime.tv_usec) / 1000000;
	cerr << "Time to render = " << seconds << endl;
	/*****************************/

	WriteImage(output, "volren");
}
