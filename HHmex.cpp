//int i            = MxArray(prhs[0]).toInt();
//double d         = MxArray(prhs[0]).toDouble();
//bool b           = MxArray(prhs[0]).toBool();
//std::string s    = MxArray(prhs[0]).toString();
//cv::Mat mat      = MxArray(prhs[0]).toMat();   // For pixels
//cv::Mat ndmat    = MxArray(prhs[0]).toMatND(); // For N-D array
//cv::Point pt     = MxArray(prhs[0]).toPoint();
//cv::Size siz     = MxArray(prhs[0]).toSize();
//cv::Rect rct     = MxArray(prhs[0]).toRect();
//cv::Scalar sc    = MxArray(prhs[0]).toScalar();
//cv::SparseMat sp = MxArray(prhs[0]).toSparseMat(); // Only double to float
//
//plhs[0] = MxArray(i);
//plhs[0] = MxArray(d);
//plhs[0] = MxArray(b);
//plhs[0] = MxArray(s);
//plhs[0] = MxArray(mat);
//plhs[0] = MxArray(ndmat);
//plhs[0] = MxArray(pt);
//plhs[0] = MxArray(siz);
//plhs[0] = MxArray(rct);
//plhs[0] = MxArray(sc);
//plhs[0] = MxArray(sp); // Only 2D float to double

//double teste = scale.at<double>(10);	//Acessa valor de matriz 1D
//double teste2 = im.at<double>(0, 1, 2);	// Acessa matriz 3D -> Equivalente no Matlab a teste2(3,2,1)
//cv::Mat teste2 = im.row(0).clone();	// Copia linha
//cv::Mat teste3 = teste2.reshape(0, 3); //transforma linha em matrix 3xalgo (channels, rows)
//double teste4 = teste3.at<double>(1, 2); //equivalente no Matlab a teste3(3,2)
//cv::Mat teste5 = Mat::zeros(3, 3, CV_64F);  // Cria matriz zerada
//im.row(0).copyTo(teste5.row(0));	//Copia Linha


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
    if (nrhs!=4)
        mexErrMsgIdAndTxt("myfunc:invalidArgs", "Wrong number of arguments");
    
    // Convert MxArray to cv::Mat

    cv::Mat layers(rhs[0].toMat());
	cv::Mat x_shift(rhs[1].toMat());
	cv::Mat y_shift(rhs[2].toMat());
	cv::Mat sizes(rhs[3].toMat());
	
	cv::Mat layer0 = layers(Rect(0,0,1,sizes.at<double>(0)*sizes.at<double>(0))).clone().reshape(0, sizes.at<double>(0));
	//Crop -> Rect(X,Y,Width,Height)
	cv::Mat layer1 = layers(Rect(0,sizes.at<double>(0)*sizes.at<double>(0),1,sizes.at<double>(1)*sizes.at<double>(1))).clone().reshape(0, sizes.at<double>(1));
	
	int crop = (sizes.at<double>(1)-sizes.at<double>(0))/2;	
	
    // #WFS 0
	int i = 0;
	cv::Mat temp;
	shift(layer1,temp,Point2f(y_shift.at<double>(i),x_shift.at<double>(i)), BORDER_CONSTANT, 0);
	cv::Mat out0 = temp(Rect(crop,crop,sizes.at<double>(0),sizes.at<double>(0)));	//Crop -> Rect(X,Y,Width,Height)
	out0 += layer0;
	
    // #WFS 1
	i = 1;
	shift(layer1,temp,Point2f(y_shift.at<double>(i),x_shift.at<double>(i)), BORDER_CONSTANT, 0);
	cv::Mat out1 = temp(Rect(crop,crop,sizes.at<double>(0),sizes.at<double>(0)));	//Crop -> Rect(X,Y,Width,Height)
	out1 += layer0;
	
    // #WFS 2
	i = 2;
	shift(layer1,temp,Point2f(y_shift.at<double>(i),x_shift.at<double>(i)), BORDER_CONSTANT, 0);
	cv::Mat out2 = temp(Rect(crop,crop,sizes.at<double>(0),sizes.at<double>(0)));	//Crop -> Rect(X,Y,Width,Height)
	out2 += layer0;
	
    // #WFS 3
	i = 3;
	shift(layer1,temp,Point2f(y_shift.at<double>(i),x_shift.at<double>(i)), BORDER_CONSTANT, 0);
	cv::Mat out3 = temp(Rect(crop,crop,sizes.at<double>(0),sizes.at<double>(0)));	//Crop -> Rect(X,Y,Width,Height)
	out3 += layer0;
	
    // #WFS 4
	i = 4;
	shift(layer1,temp,Point2f(y_shift.at<double>(i),x_shift.at<double>(i)), BORDER_CONSTANT, 0);
	cv::Mat out4 = temp(Rect(crop,crop,sizes.at<double>(0),sizes.at<double>(0)));	//Crop -> Rect(X,Y,Width,Height)
	out4 += layer0;
	
	cv::Mat out;
	hconcat(out0.clone().reshape(0,1),out1.clone().reshape(0,1),out);
	hconcat(out.clone().reshape(0,1),out2.clone().reshape(0,1),out);
	hconcat(out.clone().reshape(0,1),out3.clone().reshape(0,1),out);
	hconcat(out.clone().reshape(0,1),out4.clone().reshape(0,1),out);
	
	cv::transpose(out, out);
	
	plhs[0] = MxArray(out);
}