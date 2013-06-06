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

#ifndef _EXCEPTION_H
#define _EXCEPTION_H

#include <string>
#include <iostream>

using namespace std;

namespace tdkp {
	

	/** base exception class */
	class Exception {
		
	public:
		Exception();
		Exception(const string& reason, int line, const char* pfile, const char* pdate, const char* ptime);
		Exception(const string& reason, int line, const char* pfile, const char* pdate, const char* ptime, const char* pfunc);
		Exception(const char* preason,  int line, const char* pfile, const char* pdate, const char* ptime, const char* pfunc);
		Exception(const char* preason,  int line, const char* pfile, const char* pdate, const char* ptime);
		Exception(unsigned int buflen, const char* preason, ...);
		~Exception();
		
		void write_to_log() const;
						
		const string& get_reason() const;		
			
	protected:
		void append_info(int line, const char* pfile, const char* pdate, const char* ptime, const char* pfunc = 0);
		string reason;			
				
	};
	/** file exceptions */
	class FileException : public Exception {
	public:
		FileException(const char* preason, const char* filename, int line, const char* pfile, const char* pdate, const char* ptime);
			
	};
	
	class ParserException : public Exception {
	public:	
		ParserException(const char* preason, const char* filename, int line_number, string& line);
	};
			
	ostream& operator<<(std::ostream &out, const Exception& pe);
		
}
#endif 
