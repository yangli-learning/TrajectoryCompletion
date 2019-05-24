#include <QApplication>

//#include <glog/logging.h>
//#include "caffe_wrapper.h"

#include "main_window.h"
#include "guiless.h"
#include <gflags/gflags.h>
#include <ostream>
#include <string>
using namespace std;

string ti="Turns on interactive mode",
	tb="Path to bbox file (.yaml format)",
	tp="Path to parameter file ",
	default_path="../GPStudio/parameters.ini",
	tv = "Minimum voting image threshold";
double default_v=0.08;

DEFINE_bool(interactive, true, ti.c_str());
DEFINE_string(trace_fname,  "", ti.c_str());
DEFINE_string(bbox_fname, "", tp.c_str()); //"Path to bbox file (.yaml format");
DEFINE_string(param_fname, default_path, tp.c_str()); //"../GPStudio/parameters.ini","Path to parameter file");
DEFINE_double(vmin,default_v,tv.c_str());



bool loadParameters() {
	std::string filename(FLAGS_param_fname);
	Parameters::getInstance()->setParameterFile(filename);
	if (! Parameters::getInstance()->readParameters()){
		std::cout << "Can not locate parameter file for GPStudio!" << std::endl;
		std::cout << "Expected location: \"" << filename << "\"" << std::endl;
		return false;
	}
	return true;

}

int main(int argc, char *argv[])
{
	//::google::InitGoogleLogging("Image2Scene");
	qDebug() << "------------------" << endl;
	google::ParseCommandLineFlags(&argc,&argv,true);//flags

	qDebug() << "------------------" << endl;
	
	if (FLAGS_interactive){
		// run with gui, use parameter from default file
		
		google::ShutDownCommandLineFlags();
		// Make Xlib and GLX thread safe under X11
		QApplication::setAttribute(Qt::AA_X11InitThreads);
		QApplication application(argc, argv);

		MainWindow::getInstance()->show();//Maximized();

	    application.exec();

	}else{ //run without gui, use parameters from
	  
		// read command line options to input, output trajectory,
		// bbox dimension
		// (including voting image and heading image, hence no need for
		// computing voting map)
		cout << "parameter file :" << FLAGS_param_fname << endl;
		
		cout << "trajectory file :" << FLAGS_trace_fname << endl;
		cout << "bounding box file :" << FLAGS_bbox_fname << endl;
		//		cout << "c_vmin file :" << FLAGS_vmin << endl;

		if (!loadParameters()){
			cout <<"Not loading parameters" << endl;
			google::ShutDownCommandLineFlags();
			//return false;
		}
		

		if (FLAGS_vmin){ //v_min overwrite thresh in parameters.ini
			cout << "minimum vote threshold :" << FLAGS_vmin << endl;
		}
		
		Guiless app(FLAGS_trace_fname , FLAGS_bbox_fname);
		google::ShutDownCommandLineFlags();
		 app.exec();
		
	}

}
