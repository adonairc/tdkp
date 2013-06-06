// ------------------------------------------------------------
//
// This file is part of tdkp, a simulation tool for nanostrutctures
// of optoelectronics developed at ETH Zurich
//
// (C) 2005-2009 Ratko G. Veprek, ETH Zurich, veprek@iis.ee.ethz.ch
//
// 1) As of 18.6.2009 this code is property of ETH Zurich and must not be
// transferred, modified or used by third parties without appropriate
// licenses issued by authorized agents of ETH Zurich.
//
// 2) Violation of this will result in judicial action according to civil
// and penal law.
//
// 3) Any claim of authorship other than by the author himself is
// strictly forbidden.
//
// 4) The source code must retain the copyright notice, this list
// of conditions and the following disclaimer.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS
// BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
// IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ------------------------------------------------------------


#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "tcl.h"
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Configuration.h"
#include "tdkp/common/MPIJobManager.h"
#include "bin/Version.h"


using namespace std;
using namespace tdkp;

extern "C" {
	int Tdkpshell_Init(Tcl_Interp *interp);
	
}
namespace tdkp {
void special();
}
// --------------------------------------------------------
// help and logging utilities
// --------------------------------------------------------
void why() {	
	Logger::get_instance()->emit(LOG_INFO, "weil das gut so ist");
	special();	
}

ofstream outputfile;

tdkp::Logger* get_logger_object() {
	return Logger::get_instance();	
}

void logger_output_to_file(const char* filename) {
	
    Logger::get_instance()->del_listener(&outputfile);
	outputfile.close();
	outputfile.clear();
			
	outputfile.open(filename, ios::out | ios::app);
	if(outputfile) {
		Logger::get_instance()->add_listener(&outputfile);
	} else {
		TDKP_GENERAL_EXCEPTION("can not write to file " << filename);
	}					
	
}

tdkp::Configuration* get_configuration_object() {
	return Configuration::get_instance();	
}

void help(const char* command) {
	
}

void logger_emit(const char* message, tdkp::LoggerLevel level = LOG_INFO) {
	Logger::get_instance()->emit(level, message);	
}

void logger_emit(const string& message, tdkp::LoggerLevel level = LOG_INFO) {
	Logger::get_instance()->emit(level, message);		
}

const char* string2char(const string& string) {
	return string.c_str();	
}

string char2string(const char* str) {
	return string(str); 	
}

complex<double> double2complex(double real, double imag) {
	return complex<double>(real,imag);	
}

double complex2real(const complex<double>& c) {
	return c.real();	
}
double complex2imag(const complex<double>& c) {
	return c.imag();	
}

complex<double> eval_cref(const complex<double>& c) {
	return c;	
}
double eval_dref(const double& d) {
	return d;	
}
double eval_dref(const double* d) {
	return *d;	
}

const cplx& vec_get(const vector<cplx>& vec, unsigned int idx) {
	TDKP_ASSERT(vec.size() > idx, "");
	return vec[idx];
}

void help() {
	ostringstream sout;
	sout << "tdkp shell - shell interface to 1d/2d/3d kp finite element solver.\n";
	Logger::get_instance()->emit(LOG_INFO, sout.str().c_str());	     
}

const char* explanations[] = {
	"pid  - The process id",
	"comm - The filename of the executable, in parentheses.",
	"state - One character from the string RSDZTW",
	"pid - The PID of the parent.",
	"grp - The process group ID of the process.",
	"ession - The session ID of the process.",
	"ty_nr - The tty the process uses.",
	"pgid - The process group ID of the process which currently owns the tty that the process is connected to.",
	"lags - The flags of the process.  The math bit is decimal 4, and the traced bit is decimal 10.",
	"minflt - The number of minor faults the process has made which have not required loading a memory page from disk.",
	"cminflt - The number of minor faults that the process’s waited-for children have made.",
	"majflt - The number of major faults the process has made which have required loading a memory page from disk.",
	"cmajflt - The number of major faults that the process’s waited-for children have made.",
	"utime - The number of jiffies that this process has been scheduled in user mode.",
	"stime - The number of jiffies that this process has been scheduled in kernel mode.",
	"cutime - The  number  of  jiffies  that  this  process’s  waited-for children have been scheduled in user mode.",
	"cstime - The number of jiffies that this process’ waited-for children have been scheduled in kernel mode.",
	"priority - The standard nice value, plus fifteen.  The value is never negative in the kernel.",
	"nice - The nice value ranges from 19 (nicest) to -19 (not nice to others).",
	"0 - This value is hard coded to 0 as a placeholder for a removed field.",
	"itrealvalue - The time in jiffies before the next SIGALRM is sent to the process due to an interval timer.",
	"starttime - The time in jiffies the process started after system boot.",
	"vsize - Virtual memory size in bytes.",
	"rss - Resident Set Size: number of pages the process has in real memory",
	"rlim - Current limit in bytes on the rss of the process (usually 4294967295 on i386).",
	"startcode - The address above which program text can run.",
	"endcode - The address below which program text can run.",
	"startstack - The address of the start of the stack.",
	"kstkesp - The current value of esp (stack pointer), as found in the kernel stack page for the process.",
	"kstkeip - The current EIP (instruction pointer).",
	"signal - The bitmap of pending signals (usually 0).",
	"blocked - The bitmap of blocked signals (usually 0, 2 for shells).",
	"sigignore - The bitmap of ignored signals.",
	"sigcatch - The bitmap of catched signals.",
	"wchan - This is the channel in which the process is waiting.",
	"nswap - Number of pages swapped - not maintained.",
	"cnswap - Cumulative nswap for child processes.",
	"exit_signal - Signal to be sent to parent when we die.",
	"processor - CPU number last executed on."
};

void write_proc_self_stats() {
	
	ifstream fin("/proc/self/stat");
	string tmp;
	ostringstream sout;
	// advance 9
	for(int ii = 0; ii < 39; ii++) {
		fin >> tmp;
		sout << setw(2) << ii << "  " << setw(25) << tmp << " " << explanations[ii] << "\n";	
	} 
	cout << " ------------------- /proc/self/stat ------------------\n" 
	     << sout.str()
	     << "\n-------------------------------------------------------\n";
	
}

vector<double> read_vector_from_file(const char* filename) {
	vector<double> ret;
	ifstream fin(filename);
	if(fin) {
		double d;		
		while(!fin.eof()) {
			fin >> d;
			ret.push_back(d);	
		}
		return ret;
	} else {
		std::cerr << "can not read from file " << filename << "\n";
		return ret;	
	}	
}

void write_lumi_vector_to_file(double xmin, double xmax, int dsets, const vector<double>& data, const char* filename) {
	ofstream fout(filename);	
	if(fout) {
		unsigned int dlength = data.size() / dsets;
		TDKP_ASSERT(dsets * dlength == data.size(),"");
		double dx = dlength > 1 ? (xmax - xmin) / (dlength - 1) : 0;
		for(unsigned int ii = 0; ii < dlength; ii++) {
			fout << setw(15) << xmin + dx * ii << "  ";
			for(int jj = 0; jj < dsets; jj++) {
				fout << setw(15) << data[ii * dsets + jj] << " ";
			}
			fout << "\n";	
		}
		fout.close();
	} else {
		std::cerr << "can not write to file " << filename << "\n";		
	}	
}

int Tcl_AppInit(Tcl_Interp * interp)
{
  if (Tcl_Init(interp) != TCL_OK) {
    cerr << "tcl init failed\n";
    return TCL_ERROR;
  }

  if (Tdkpshell_Init(interp) != TCL_OK) {
  	cerr << "uiuiui: hellshell init failed\n";
    return TCL_ERROR;
  }

  ostringstream sout;
  sout << "set tcl_prompt1 {puts -nonewline tdkp\\>\\  }";
  Tcl_Eval(interp,sout.str().c_str());
  return TCL_OK;
}

void quit() {
	
	// --------------------------------------------------
	// finalize (disconnect other mpi jobs)
	// --------------------------------------------------
	MPIJobManager::get_instance().finalize();
	
	TimeMeasurements::get_instance().stop("overall");	
	ostringstream sout;
	sout << TimeMeasurements::get_instance();
	Logger::get_instance()->emit(LOG_INFO, sout.str());	
	
	cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
	     << " 						t h a n k   y o u\n"
	     << " 						  for using tdkp\n"
	     << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
	     
	//write_proc_self_stats();	     

	if(outputfile) {		
		Logger::get_instance()->del_listener(&outputfile);
		outputfile.close();
	}
	     	
	Tcl_Exit(0);
	
}



int main(int argc, char* argv[]) {
	
	TimeMeasurements::get_instance().start("overall");	
	
	// ---------------------------------------------------
	// start MPI manager
	// ---------------------------------------------------
	MPIJobManager::create_instance(&argc, &argv);
	MPIJobManager::get_instance().mpi_loop();
		
	// ---------------------------------------------------
	// if this is rank 0 job, the mpi job manager will 
	// exit and we proceed as is. if this job is rank > 0,
 	// the mpi manager never leaves the constructor ...
 	// ---------------------------------------------------
	
	cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
	     << "  tdkp - a multidimensional finite element kp solver \n"
	     << "\n"
	     << "  (c) 2006-2009 ratko veprek, computational optoelectronics\n"	     
	     << "                integrated systems laboratory, eth zuerich\n"
	     << "\n"
	     << "  you are using tdkp revision " << TDKP_REVISION << ", compiled \n"
	     << "  on " << MACHINE << " at " << __DATE__ << ", " << __TIME__ << "\n"
	     << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
	     	
	Tcl_Main(argc, argv, Tcl_AppInit);

}

