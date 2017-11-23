
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
	
	double data[13][13] = {{0.000465941243316438,0.000743031610361572,0.00161656632438361,0.00472734218198899,0.0123178362186659,0.0231723083055391,0.0291640100745752,0.0231723073217654,0.0123178363481702,0.00472734142716301,0.00161656686301893,0.000743032344166441,0.000465940416718631},
		{0.000743032417026835,0.00164706415650957,0.00208658475302679,-0.00673875615654296,-0.0441569275412257,-0.113423306429761,-0.158716986455740,-0.113423305929399,-0.0441569271009932,-0.00673875520154683,0.00208658402466987,0.00164706392149224,0.000743032409467637},
		{0.00161656594697829,0.00208658463566975,-0.00177404164516209,0.0315563359198349,0.250032372921486,0.759977420348036,1.16388160714302,0.759977420938588,0.250032372085299,0.0315563354362360,-0.00177404173716947,0.00208658538701729,0.00161656596510309},
		{0.00472734142826282,-0.00673875592521003,0.0315563361296266,0.124034986955866,-0.734027019513579,-4.15531919152376,-7.86893872383399,-4.15531919259779,-0.734027018408397,0.124034987199335,0.0315563361953433,-0.00673875615072542,0.00472734182565767},
		{0.0123178359085543,-0.0441569272837529,0.250032373023320,-0.734027019112402,-0.428280821835905,22.0023125853292,61.4190684274111,22.0023125854200,-0.428280821966375,-0.734027020011068,0.250032373139896,-0.0441569280971435,0.0123178365064532},
		{0.0231723078247648,-0.113423305184675,0.759977419345479,-4.15531919208594,22.0023125859559,-9.18005327260989,-323.431054418957,-9.18005327221266,22.0023125856875,-4.15531919116582,0.759977420367020,-0.113423305899084,0.0231723080880550},
		{0.0291640106602170,-0.158716988292153,1.16388160831662,-7.86893872342928,61.4190684266374,-323.431054418805,968.988423594796,-323.431054419536,61.4190684270099,-7.86893872453772,1.16388160752838,-0.158716987365392,0.0291640103367434},
		{0.0231723079663817,-0.113423305151804,0.759977419699454,-4.15531919229161,22.0023125860733,-9.18005327286173,-323.431054419513,-9.18005327157221,22.0023125854202,-4.15531919135330,0.759977420311582,-0.113423305637743,0.0231723079825611},
		{0.0123178362900924,-0.0441569279839492,0.250032372850100,-0.734027018365024,-0.428280822510867,22.0023125858836,61.4190684277843,22.0023125853107,-0.428280822081424,-0.734027019160660,0.250032372526149,-0.0441569276901481,0.0123178360631572},
		{0.00472734111986447,-0.00673875555258005,0.0315563355633255,0.124034987547252,-0.734027019431822,-4.15531919152258,-7.86893872521336,-4.15531919216987,-0.734027018482334,0.124034986793980,0.0315563361017125,-0.00673875558224119,0.00472734176110177},
		{0.00161656663382599,0.00208658455611165,-0.00177404133497322,0.0315563349316723,0.250032373370607,0.759977419583207,1.16388160943896,0.759977419830150,0.250032373300812,0.0315563348756839,-0.00177404095059861,0.00208658346676090,0.00161656696275726},
		{0.000743031711433262,0.00164706399792697,0.00208658474046763,-0.00673875580771748,-0.0441569273761150,-0.113423305877774,-0.158716988450934,-0.113423305302550,-0.0441569283175138,-0.00673875517013869,0.00208658475022243,0.00164706453558475,0.000743032402131483},
		{0.000465941303905772,0.000743032582697057,0.00161656614572295,0.00472734167067749,0.0123178365558675,0.0231723081520764,0.0291640103393851,0.0231723080725842,0.0123178368416122,0.00472734108280510,0.00161656621560976,0.000743032205066957,0.000465940168428929}};

	
	cv::Mat h = Mat(13, 13, CV_64F, data);
	
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