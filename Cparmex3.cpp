//The same as Cmex, but receiving an external kernel
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
    if (nrhs!=3)
        mexErrMsgIdAndTxt("myfunc:invalidArgs", "Wrong number of arguments");
    
    // Convert MxArray to cv::Mat

    cv::Mat layers(rhs[0].toMat());
	cv::Mat sizes(rhs[1].toMat());
	cv::Mat h(rhs[2].toMat());
	
    int sizes_temp = sizes.at<double>(0)*sizes.at<double>(0);
	cv::Mat layer0 = layers(Rect(0,0,1,sizes_temp)).clone().reshape(0, sizes.at<double>(0));
	//Crop -> Rect(X,Y,Width,Height)
	cv::Mat layer1 = layers(Rect(0,sizes_temp,1,sizes.at<double>(1)*sizes.at<double>(1))).clone().reshape(0, sizes.at<double>(1));
    sizes_temp += sizes.at<double>(1)*sizes.at<double>(1);
    cv::Mat layer2 = layers(Rect(0,sizes_temp,1,sizes.at<double>(2)*sizes.at<double>(2))).clone().reshape(0, sizes.at<double>(2));
	
	
	cv::Mat layer0f;
	filter2D(layer0,layer0f,-1,h,cv::Point(-1,-1),0,cv::BORDER_CONSTANT);
    //filter2D(layer0,layer0f,-1,h,cv::Point(-1,-1),0,cv::BORDER_DEFAULT);
	//void filter2D(InputArray src, OutputArray dst, int ddepth, InputArray kernel, Point anchor=Point(-1,-1), double delta=0, int borderType=BORDER_DEFAULT )
	
	cv::Mat layer1f;
	filter2D(layer1,layer1f,-1,h,cv::Point(-1,-1),0,cv::BORDER_CONSTANT);
    //filter2D(layer1,layer1f,-1,h,cv::Point(-1,-1),0,cv::BORDER_DEFAULT);
    
    cv::Mat layer2f;
	filter2D(layer2,layer2f,-1,h,cv::Point(-1,-1),0,cv::BORDER_CONSTANT);
    
	
	cv::Mat out;
	hconcat(layer0f.clone().reshape(0,1),layer1f.clone().reshape(0,1),out);
    hconcat(out.clone().reshape(0,1),layer2f.clone().reshape(0,1),out);
	cv::transpose(out, out);
	
	plhs[0] = MxArray(out);
}