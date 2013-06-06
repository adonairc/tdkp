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

#ifndef _LOGGER_H
#define _LOGGER_H

#include <iostream>
#include <list>
#include <string>

#include "tdkp/common/all.h"
#include "tdkp/common/Exception.h"


namespace tdkp {

	// -----------------------------------------------
	// log message level definition
	// -----------------------------------------------
	enum LoggerLevel {
		LOG_ERROR 	    = 1,
		LOG_WARN		= 2,
		LOG_INFO		= 3,
		LOG_WARN_DEVEL  = 4,
		LOG_INFO_DEVEL1 = 5,
		LOG_INFO_DEVEL2 = 6
	};
		
	/** manage programm output */
	class Logger{
	public:
	
		static Logger* get_instance();
		static Logger* get_instance(LoggerLevel resetlevel);
			
				
		void add_listener(std::ostream*);
		void del_listener(std::ostream*);
		
		void set_level(LoggerLevel level_);
		LoggerLevel get_level() const { return this->log_level; }
				
//		void emit(const std::string* ps) const;		
		void emit(const std::string& str) const;
		void emit(const char* msg) const;
//		void emit(const Exception &) const;
		void emit(const Exception* ) const;
		void emit(LoggerLevel level, const char* msg) const;
		void emit(LoggerLevel level, const std::string& str) const;
//		void emit(LoggerLevel level, const std::string* str) const;		
		void emit(LoggerLevel level, unsigned int buflen, const char* fmt, ...) const;
		
		void delimiter();		
		void set_line_to_stdout(const std::string& text);
		void init_progress_bar(const char* text, int total);		
		void init_progress_bar(const std::string &text, int total);
		int  set_progress_bar(int current, int total);
		void end_progress_bar(); 
								
	protected:
		std::list<std::ostream*> outstreams;
		LoggerLevel log_level;
		bool stdout_registered;
		bool progress_bar_active;
		
		string wrapped_write(const string& out) const;
		/** singleton class (only one object may exist) */
		static Logger* singleton;	
		Logger(LoggerLevel level_ = LOG_INFO);
		~Logger();
			
	};	
		
}

#endif
