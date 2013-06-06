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


#include <stdarg.h>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include "tdkp/common/Logger.h"
#include "tdkp/common/all.h"
#include <omp.h>
#include <algorithm>

namespace tdkp {

Logger::Logger(LoggerLevel level_) {
	if(level_ > LOG_INFO_DEVEL2) {
		level_ = LOG_INFO_DEVEL2;	
	}
	this->log_level = level_;
	this->stdout_registered   = false;
	this->progress_bar_active = false;
}

Logger::~Logger() {

}

Logger* Logger::singleton = 0;

Logger* Logger::get_instance() {
	
	#pragma omp critical (logmsg_generate)
	{
		if(Logger::singleton == 0) {
			Logger::singleton = new Logger(); 	
			Logger::singleton->add_listener(&std::cout);	 
		}
	}
	
	return Logger::singleton;
}
Logger* Logger::get_instance(LoggerLevel resetlevel) {
	#pragma omp critical (logmsg_generate)
	{
		if(Logger::singleton == 0) {
			Logger::singleton = new Logger(); 	
			Logger::singleton->add_listener(&std::cout);	 
		}
		Logger::singleton->set_level(resetlevel);
	} 		
	return Logger::singleton;
}


/** allows resetting output level */
void Logger::set_level(LoggerLevel level_) {
	TDKP_ASSERT(level_ >= LOG_ERROR && level_ <= LOG_INFO_DEVEL2, "invalid loglevel set");
	#pragma omp critical (logmsg_set_level)
	{
		this->log_level = level_;
	}	
}

/** attach outstream as listener */
void Logger::add_listener(std::ostream* stream) {
	#pragma omp critical (logmsg_add_listener)
	{
		if(find(this->outstreams.begin(), this->outstreams.end(), stream) == this->outstreams.end()) {
			this->outstreams.push_back(stream);		
			if(stream == &std::cout) {
				this->stdout_registered = true;				
			}
		}
	}
}
/** detach outstream as listener */
void Logger::del_listener(std::ostream* stream) {
	#pragma omp critical (logmsg_del_listener)
	{	
		this->outstreams.remove(stream);	
		if(stream == &std::cout) {
			this->stdout_registered = false;	
		}
	}
}

void Logger::emit(const std::string& str) const {
	#pragma omp critical (logmsg_emit)
	{	
		string tmp = wrapped_write(str);
		for(list<std::ostream*>::const_iterator it = this->outstreams.begin();	it != this->outstreams.end(); it++) {
			// rewind for progress bar if master thread
			if(*it == &(std::cout) && this->progress_bar_active && omp_get_thread_num() == 0) {		
				(*(*it)) << "\r";
			} else if(omp_get_num_threads() > 1 && omp_get_thread_num() != 0) {
				// inform which thread is talking
				(*(*it)) << "T" << omp_get_thread_num();
			}
			// output rest			
			(*(*it)) << tmp << endl;
		}
	}
	
}

void Logger::emit(LoggerLevel level, const std::string& str) const {
	if(level <= this->log_level) {
		ostringstream sout;
		switch(level) {
			case LOG_ERROR:
				sout << "\n\033[5;1;31m** error **\033[0;0;0m ";
				break;				
			case LOG_WARN:
				sout << "\033[0;0;1m** warning **\033[0;0;0m ";
				break;
			default:	
				sout << "-> ";
				break;			
		}
		sout << str;
		this->emit(sout.str());	
	}
}

void Logger::emit(const char* msg) const {
	string tmp(msg);
	this->emit(tmp);
}

string Logger::wrapped_write(const string& out) const {
	// first, search if there are any '\n' set
	int break_at = 80;
	if(out.find_first_of("\n") == string::npos) {
		ostringstream sout;
		int current    = 0;
		int last_ws    = 0;
		int last_write = 0;
		while(current < (signed)out.size()) {
			// find whitespace
			if(out[current] == '\t' || out[current] == ' ') {
				last_ws = current;	
			}	
			// write?
			if((current - last_write) > break_at) {
				// if there was no white space yet: break anyway
				if(last_ws > last_write) {
					current = last_ws + 1;	
				}		
				if(last_write > 0) {
					sout << "   ";	
				}
				sout << out.substr(last_write, current - last_write) << "\n";
				last_write = last_ws = current;								
			}
			current++;
		}
		sout << (last_write > 0 ? "   ":"") << out.substr(last_write);
		return sout.str();									
	} else {
		
		ostringstream sout;
		for(unsigned int ii = 0; ii < out.size(); ii++) {
			sout << out[ii];
			if(out[ii] == '\n') {
				sout << "   ";	
			}				
		}
		return sout.str();		
	}
}


/** emit exception */
/*
void Logger::emit(const Exception& e) const {
	this->emit(e.get_reason());	
}*/

void Logger::emit(const Exception* e) const {
	this->emit(e->get_reason());	
}

/** emit message with specific loglevel 
 */
void Logger::emit(LoggerLevel level, const char* msg) const {
	string tmp(msg);
	this->emit(level, tmp);
}

/** emit formatted message with specifig loglevel 
 * @param level  desired loglevel
 * @param buflen length of expected output string (MUST BE BIGGER THAN OUTPUT STRING)
 * @param fmt    format string for sprintf
 * @param ...    argument list 
 * */
void Logger::emit(LoggerLevel level, unsigned int buflen, const char* fmt, ...) const {
	if(level <= this->log_level) {				
	   va_list args;
	   va_start(args,fmt);
	   char* buf = new char[buflen];
	   vsprintf(buf, fmt, args);
	   this->emit(buf);	   	      
	   delete[] buf;
	}
}

void Logger::init_progress_bar(const char* text, int total) {
	std::string tmp(text);
	this->init_progress_bar(tmp, total);	
}

void Logger::init_progress_bar(const std::string &text, int total) {	
	if(this->stdout_registered) {
		// only for single threaded jobs or for  ...
		if(omp_get_thread_num() == 0) {
			#pragma omp critical (logmsg_emit)
			{
				std::cout << "   " << text << " " << total << "\n";	
				for(int ii = 0; ii < 70; ii++) {
					std::cout << " ";	
				}
			}		
			this->set_progress_bar(0, total);
			this->progress_bar_active = true;						
		}
	}	
}
int Logger::set_progress_bar(int current, int total) {
	
	if(omp_get_thread_num() == 0 && this->stdout_registered) {
		int prct_x = (current * 62) / total;
		int prct   = (current * 100) / total;
		#pragma omp critical (logmsg_emit)
		{		
			// delete old line
			std::cout << "\r";
	
			// write new line
			std::cout << "   [";
			for(int ii = 0; ii < prct_x; ii++) {
				std::cout << "x";	
			}
			for(int ii = prct_x; ii < 62; ii++) {
				std::cout << " ";	
			}
				
			std::cout << "] ";
			std::cout.width(3);
			std::cout << prct << "%";
		//	cout << " " << current << "/" << total;	
			std::cout.flush();
		}
		
		// calculate next
		int next = (int)floor((((float)prct + 1.0) * (float)total) / 100.0);
		return next > current ? next : (current + 1);
	}
	return total - 1;		
}

void Logger::set_line_to_stdout(const std::string& text) {
	if(omp_get_thread_num() == 0) {
		if(this->stdout_registered) {
			#pragma omp critical (logmsg_emit)
			{			
				std::cout << "\r" << text;
				std::cout.flush();
			}		
		}
	}	
}

void Logger::end_progress_bar() {
	if(omp_get_thread_num() == 0) {
		if(this->stdout_registered) {
			this->set_progress_bar(100,100);
			#pragma omp critical (logmsg_emit)
			{			
				std::cout << " ... done\n";
			}	
		}
		this->progress_bar_active = false;
	}
}

void Logger::delimiter() {
	const string delimiter("\n");
	this->emit(delimiter);	
}

} // end of namespace tdkp

