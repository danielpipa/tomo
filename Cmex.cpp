
#include "mexopencv.hpp"
#include "shift.hpp"
using namespace std;
using namespace cv;


void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
    // Argument vector
    vector<MxArray> rhs(prhs,prhs+nrhs);
    
    // Check arguments
    if (nrhs!=2)
        mexErrMsgIdAndTxt("myfunc:invalidArgs", "Wrong number of arguments");
    
    // Convert MxArray to cv::Mat

    cv::Mat layers(rhs[0].toMat());
	cv::Mat sizes(rhs[1].toMat());
	
	cv::Mat layer0 = layers(Rect(0,0,1,sizes.at<double>(0)*sizes.at<double>(0))).clone().reshape(0, sizes.at<double>(0));
	//Crop -> Rect(X,Y,Width,Height)
	cv::Mat layer1 = layers(Rect(0,sizes.at<double>(0)*sizes.at<double>(0),1,sizes.at<double>(1)*sizes.at<double>(1))).clone().reshape(0, sizes.at<double>(1));
	
	double data[5][5] = {{0.0047, -0.0237, 0.0685, -0.0237, 0.0047},
					{-0.0237, 0.1196, -0.3458, 0.1196, -0.0237},
					{0.0685, -0.3458,1.0000,-0.3458, 0.0685},
					{-0.0237, 0.1196, -0.3458, 0.1196, -0.0237},
					{0.0047, -0.0237, 0.0685, -0.0237, 0.0047}};
	
	cv::Mat h = Mat(5, 5, CV_64F, data);
	
	cv::Mat layer0f;
	filter2D(layer0,layer0f,-1,h,cv::Point(-1,-1),0,cv::BORDER_CONSTANT);
	//void filter2D(InputArray src, OutputArray dst, int ddepth, InputArray kernel, Point anchor=Point(-1,-1), double delta=0, int borderType=BORDER_DEFAULT )
	
	cv::Mat layer1f;
	filter2D(layer1,layer1f,-1,h,cv::Point(-1,-1),0,cv::BORDER_CONSTANT);

	
	cv::Mat out;
	hconcat(layer0f.clone().reshape(0,1),layer1f.clone().reshape(0,1),out);
	cv::transpose(out, out);
	
	plhs[0] = MxArray(out);
}