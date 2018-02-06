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

    cv::Mat WFS(rhs[0].toMat());
	cv::Mat x_shift(rhs[1].toMat());
	cv::Mat y_shift(rhs[2].toMat());
	cv::Mat sizes(rhs[3].toMat());
	
	int w_size = sizes.at<double>(0);
	cv::Mat WFS0 = WFS(Rect(0,0,1,w_size*w_size)).clone().reshape(0, w_size);
	//Crop -> Rect(X,Y,Width,Height)
	cv::Mat WFS1 = WFS(Rect(0,w_size*w_size,1,w_size*w_size)).clone().reshape(0, w_size);
	cv::Mat WFS2 = WFS(Rect(0,2*w_size*w_size,1,w_size*w_size)).clone().reshape(0, w_size);
	cv::Mat WFS3 = WFS(Rect(0,3*w_size*w_size,1,w_size*w_size)).clone().reshape(0, w_size);
	cv::Mat WFS4 = WFS(Rect(0,4*w_size*w_size,1,w_size*w_size)).clone().reshape(0, w_size);
	
	
	
	// #Layer 0
	cv::Mat out0 = WFS0 + WFS1 + WFS2 + WFS3 + WFS4;
	
	// #Layer 1
    int crop = (sizes.at<double>(1)-sizes.at<double>(0))/2;
    int i = 0; //Layer 1
	// #WFS 0
	cv::Mat temp;
	copyMakeBorder(WFS0, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	cv::Mat WFS01;
	shift(temp,WFS01,Point2f(-y_shift.at<double>(0,i),-x_shift.at<double>(0,i)), BORDER_CONSTANT, 0);
	// #WFS 1
	copyMakeBorder(WFS1, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	cv::Mat WFS11;
	shift(temp,WFS11,Point2f(-y_shift.at<double>(1,i),-x_shift.at<double>(1,i)), BORDER_CONSTANT, 0);
	// #WFS 2
	copyMakeBorder(WFS2, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	cv::Mat WFS21;
	shift(temp,WFS21,Point2f(-y_shift.at<double>(2,i),-x_shift.at<double>(2,i)), BORDER_CONSTANT, 0);
	// #WFS 3
	copyMakeBorder(WFS3, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	cv::Mat WFS31;
	shift(temp,WFS31,Point2f(-y_shift.at<double>(3,i),-x_shift.at<double>(3,i)), BORDER_CONSTANT, 0);
	// #WFS 4
	copyMakeBorder(WFS4, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	cv::Mat WFS41;
	shift(temp,WFS41,Point2f(-y_shift.at<double>(4,i),-x_shift.at<double>(4,i)), BORDER_CONSTANT, 0);
    cv::Mat out1 = WFS01 + WFS11 + WFS21 + WFS31 + WFS41;
    
    // #Layer 2
    crop = (sizes.at<double>(2)-sizes.at<double>(0))/2;
    i = 1;  //Layer 2
	// #WFS 0
	copyMakeBorder(WFS0, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	shift(temp,WFS01,Point2f(-y_shift.at<double>(0,i),-x_shift.at<double>(0,i)), BORDER_CONSTANT, 0);
	// #WFS 1
	copyMakeBorder(WFS1, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	shift(temp,WFS11,Point2f(-y_shift.at<double>(1,i),-x_shift.at<double>(1,i)), BORDER_CONSTANT, 0);
	// #WFS 2
	copyMakeBorder(WFS2, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	shift(temp,WFS21,Point2f(-y_shift.at<double>(2,i),-x_shift.at<double>(2,i)), BORDER_CONSTANT, 0);
	// #WFS 3
	copyMakeBorder(WFS3, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	shift(temp,WFS31,Point2f(-y_shift.at<double>(3,i),-x_shift.at<double>(3,i)), BORDER_CONSTANT, 0);
	// #WFS 4
	copyMakeBorder(WFS4, temp, crop, crop, crop, crop, BORDER_CONSTANT, 0);	//Pad array
	shift(temp,WFS41,Point2f(-y_shift.at<double>(4,i),-x_shift.at<double>(4,i)), BORDER_CONSTANT, 0);
    cv::Mat out2 = WFS01 + WFS11 + WFS21 + WFS31 + WFS41;
    
	
	cv::Mat out;
	hconcat(out0.clone().reshape(0,1),out1.clone().reshape(0,1),out);
    hconcat(out.clone().reshape(0,1),out2.clone().reshape(0,1),out);
	
	cv::transpose(out, out);
	
	plhs[0] = MxArray(out);
}