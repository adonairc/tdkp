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

#include "tdkp/probdefs/KP8x83D.h"
#include "tdkp/kpmatrices/KPMatrix8x8EndersForeman.h"
#include "tdkp/kpmatrices/KPMatrix8x8Wurtzite.h"
#include "tdkp/common/Configuration.h"

namespace tdkp {

KP8x83D::KP8x83D(const Geometry& geometry_, MaterialDatabase& material_database_)
: KPBase3D(geometry_, material_database_)
{
	this->init_3D();
	this->first_order_terms_exist = true;
}

KP8x83D::~KP8x83D() {

}

KPMatrixBase* KP8x83D::get_matrix() const {
	return new KPMatrix8x8EndersForeman();
}

KP8x83DWZ::KP8x83DWZ(const Geometry& geometry_, MaterialDatabase& material_database_)
: KPBase3D(geometry_, material_database_)
{
	this->init_3D();
	this->first_order_terms_exist = true;
}

KP8x83DWZ::~KP8x83DWZ() {

}

KPMatrixBase* KP8x83DWZ::get_matrix() const {
	return new KPMatrix8x8Wurtzite();
}


} // end namespace
