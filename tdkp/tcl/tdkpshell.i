 %module tdkpshell
 %{
 #include "tdkp/common/all.h"
 #include "tdkp/common/Configuration.h"
 #include "tdkp/common/Exception.h"
 #include "tdkp/common/Logger.h"
 #include "tdkp/main/PropertyContainer.h"
 #include "tdkp/geometry/Element.h"
 #include "tdkp/geometry/Node.h"
 #include "tdkp/geometry/Region.h"
 #include "tdkp/geometry/ElementBoundary.h"
 #include "tdkp/geometry/BoundaryCondition.h"
 #include "tdkp/geometry/Geometry.h"
 #include "tdkp/io/BaseGridReader.h"
 #include "tdkp/io/SebiseGridReader.h"
 #include "tdkp/io/OneDGridCreator.h"
 #include "tdkp/io/AsciiGridReader.h"
 #include "tdkp/io/MEDGridReader.h"
 #include "tdkp/io/InputParser.h"
 #include "tdkp/io/TCLDataIOWrapper.h"
 #include "tdkp/io/BaseDataIO.h"
 #include "tdkp/io/BinaryDataIO.h"
 #include "tdkp/io/AsciiDataIO.h"
 #include "tdkp/io/SebiseDataIO.h"
 #include "tdkp/io/MEDDataIO.h"
 #include "tdkp/common/Vector3D.h"
 #include "tdkp/povray/DF3Generator.h"
 #include "tdkp/povray/IsoSurfaceGenerator.h"
 #include "tdkp/main/MaterialDatabase.h"
 #include "tdkp/main/EigenSolution.h"
 #include "tdkp/common/DataTypes.h"
 #include "tdkp/common/ICurves.h"
 #include "tdkp/main/Bandstructure.h"
 #include "tdkp/probdefs/BulkBandstructureSolver.h"
 #include "tdkp/probdefs/EigenProblem.h"
 #include "tdkp/probdefs/EigenProblem3D.h" 
 #include "tdkp/probdefs/EffectiveMass.h"
 #include "tdkp/probdefs/NEGFEffectiveMass.h"
 #include "tdkp/probdefs/LinearProblem.h"
 #include "tdkp/probdefs/IntrinsicStrain.h"
 #include "tdkp/main/Fields.h"
 #include "tdkp/kpmatrices/KPMatrixBase.h"
 #include "tdkp/kpmatrices/KPMatrix4x4EndersForeman.h"
 #include "tdkp/kpmatrices/KPMatrix6x6EndersForeman.h"
 #include "tdkp/kpmatrices/KPMatrix8x8EndersForeman.h"
 #include "tdkp/kpmatrices/KPMatrix6x6Wurtzite.h"
 #include "tdkp/kpmatrices/KPMatrix8x8Wurtzite.h"
 
 #include "tdkp/probdefs/KPBase.h" 
 #include "tdkp/probdefs/KPBase3D.h" 
 #include "tdkp/probdefs/KP4x43D.h"
 #include "tdkp/probdefs/KP8x83D.h"
 
 #include "tdkp/probdefs/PoissonEquation.h"
 #include "tdkp/probdefs/ProblemDefinition.h"
 #include "tdkp/probdefs/KPBase1D2D.h"
 #include "tdkp/probdefs/KP4x41D2D.h"
 #include "tdkp/probdefs/KP6x61D2D.h"
 #include "tdkp/probdefs/KP8x81D2D.h"
 #include "tdkp/probdefs/SchroedingerPML.h"
 #include "tdkp/utilities/SchroedingerPoisson.h"
 
 #include "tdkp/probdefs/PolarizationCharge.h"
 #include "tdkp/utilities/MatrixElements.h"
 #include "tdkp/utilities/DensityOfStates.h"
 #include "tdkp/utilities/GridInterpolator.h"
 #include "tdkp/common/Domain.h"
 #include "tdkp/utilities/SLC.h"
 #include "tdkp/utilities/GraphReordering.h"
 #include "tdkp/utilities/Fermi.h"
 #include "tdkp/coulomb/CoulombMatrixElement.h"
 #include "tdkp/utilities/Overlap.h"
 #include "tdkp/clc/CLCSHF.h"
 #include "tdkp/clc/CLCScreening.h"
 
 extern void why();
 extern void help();
 using namespace tdkp;
 extern tdkp::Logger* get_logger_object();
 extern void logger_output_to_file(const char* filename);
 extern tdkp::Configuration* get_configuration_object(); 
 extern void logger_emit(const char* message, tdkp::LoggerLevel level = LOG_INFO);
 extern void logger_emit(const string& message, tdkp::LoggerLevel level = LOG_INFO);
 extern void quit();
 extern complex<double> double2complex(double real, double imag);
 extern double complex2real(const complex<double>& c);
 extern double complex2imag(const complex<double>& c); 
 extern complex<double> eval_cref(const complex<double>& c);
 extern double eval_dref(const double& d);
 extern double eval_dref(const double* d);
 extern const complex<double>& vec_get(const vector<complex<double> >& vec, unsigned int idx);
 extern const char* string2char(const string& string);
 extern string char2string(const char*);
 extern vector<double> read_vector_from_file(const char* filename);
 extern void write_lumi_vector_to_file(double xmin, double xmax, int dsets, const vector<double>& data, const char* filename);

  %}

 extern void why();
 extern void help();
 extern tdkp::Logger*        get_logger_object();
 extern tdkp::Configuration* get_configuration_object();
 extern void logger_output_to_file(const char* filename);
 extern void quit();
 extern complex<double> double2complex(double real, double imag);
 extern complex<double> double2complex(const double& real, const double& imag);
 extern double          complex2imag(const complex<double>& c);
 extern double          complex2real(const complex<double>& c);
 extern complex<double> eval_cref(const complex<double>& c);
 extern double          eval_dref(const double& d);
 extern double          eval_dref(const double* d);
 extern const complex<double>& vec_get(const vector<complex<double> >& vec, unsigned int idx);
 extern const char* string2char(const string& string);
 extern string char2string(const char*);
 extern void logger_emit(const char* message, tdkp::LoggerLevel level = LOG_INFO);
 extern void logger_emit(const string& message, tdkp::LoggerLevel level = LOG_INFO);
 extern vector<double> read_vector_from_file(const char* filename);

 extern void write_lumi_vector_to_file(double xmin, double xmax, int dsets, const vector<double>& data, const char* filename);
 
 %include "tdkp/common/all.h"

 %include "tdkp/common/Logger.h"
 %include "tdkp/common/Configuration.h"


 %include "tdkp/main/PropertyContainer.h"
 %template (PropertyContainerdouble) tdkp::PropertyContainer<double>;
 %include "tdkp/common/Vector3D.h"
 %include "tdkp/common/Domain.h"
 %include "tdkp/geometry/Node.h"
 %include "tdkp/main/MaterialDatabase.h"
 %include "tdkp/geometry/Region.h"
 %include "tdkp/geometry/Element.h" 
 %include "tdkp/geometry/ElementBoundary.h"
 %include "tdkp/geometry/BoundaryCondition.h"
 %include "tdkp/geometry/Geometry.h"
 
 %include "tdkp/common/DataTypes.h"
 %template (ElementDatadouble)    tdkp::ElementData<double>;
 %template (ElementDatacomplex)   tdkp::ElementData<complex<double> >;
 %template (StdElementDatadouble)  tdkp::StdElementData<double>;
 %template (StdElementDatacomplex) tdkp::StdElementData<complex<double>  >;
 %template (NodeDatadouble)     tdkp::NodeData<double>;
 %template (NodeDatacomplex)    tdkp::NodeData<complex<double>  >;
 %template (MergedElementDatadouble)  tdkp::MergedElementData<double>;
 %template (MergedElementDatacomplex) tdkp::MergedElementData<complex<double> >;
 %template (XYDatacomplex)        tdkp::XYData<double>;
 %template (XYDatadouble)         tdkp::XYData<complex<double>  >;
 %template (XYDataRealcomplex)    tdkp::XYDataReal<complex<double>  >; 
 %template (XYDataScaledouble)    tdkp::XYDataScale<double>;
 %template (XYDataScalecomplex)   tdkp::XYDataScale<complex<double> >; 
 %template (StdNodeDatadouble)    tdkp::StdNodeData<double>;
 %template (StdNodeDatacomplex)   tdkp::StdNodeData<complex<double>  >;
 %include "tdkp/common/ICurves.h"
 %include "tdkp/utilities/GraphReordering.h"
 %include "tdkp/io/BaseGridReader.h"
 %include "tdkp/io/SebiseGridReader.h"
 %include "tdkp/io/OneDGridCreator.h"
 %include "tdkp/io/MEDGridReader.h"
 %include "tdkp/io/InputParser.h"
 %include "tdkp/io/AsciiGridReader.h"
 %template (write_binary)  tdkp::InputParser::write_binary<double>;
 %template (write_binary)  tdkp::InputParser::write_binary<complex<double>  >;
 %template (read_binary)   tdkp::InputParser::read_binary<double>;
 %template (read_binary)   tdkp::InputParser::read_binary<complex<double>  >;
 %template (write_ascii_double)   tdkp::InputParser::write_ascii<double>;
 %template (write_ascii_cplx )     tdkp::InputParser::write_ascii<complex<double>  >; 
 %template (write_ascii_cplx2real)  tdkp::InputParser::write_ascii_real<complex<double>  >;
 
 %include "tdkp/io/AsciiDataIO.h"
 %include "tdkp/io/BinaryDataIO.h"
 %include "tdkp/io/SebiseDataIO.h"
 %include "tdkp/io/MEDDataIO.h"
 %include "tdkp/io/TCLDataIOWrapper.h"
 %template (TCLDataIOAscii)		tdkp::TCLDataIOWrapper<AsciiDataIO>;
 %template (TCLDataIOBinary)	tdkp::TCLDataIOWrapper<BinaryDataIO>;
 %template (TCLDataIOSebise)	tdkp::TCLDataIOWrapper<SebiseDataIO>;
 %template (TCLDataIOMED)		tdkp::TCLDataIOWrapper<MEDDataIO>;
 
 %template (RotationMatrix)          tdkp::RMatrix<double>;
 %include "tdkp/povray/DF3Generator.h"
 %include "tdkp/povray/IsoSurfaceGenerator.h"
 %include "tdkp/main/EigenSolution.h"
 %template (EigenSolutioncomplex)    tdkp::EigenSolution<complex<double>  >;
 %template (EigenSolutiondouble)     tdkp::EigenSolution<double>; 
 %template (EigenSolutionSetdouble)  tdkp::EigenSolutionSet<double>;
 %template (EigenSolutionSetcomplex) tdkp::EigenSolutionSet<complex<double>  >;

 %feature("notabstract") StrainField;
 %include "tdkp/main/Fields.h"

 %include "tdkp/main/Bandstructure.h"
 %template (BandstructureComplex) tdkp::Bandstructure<complex<double>  >;
 %template (BandstructureDomaincomplex) tdkp::BandstructureDomain<complex<double>  >;
  
 %include "tdkp/kpmatrices/KPMatrixBase.h"
 %include "tdkp/kpmatrices/KPMatrix4x4EndersForeman.h"
 %include "tdkp/kpmatrices/KPMatrix6x6EndersForeman.h"
 %include "tdkp/kpmatrices/KPMatrix8x8EndersForeman.h"
 %include "tdkp/kpmatrices/KPMatrix6x6Wurtzite.h"
 %include "tdkp/kpmatrices/KPMatrix8x8Wurtzite.h"
 %include "tdkp/probdefs/BulkBandstructureSolver.h"  
 %include "tdkp/probdefs/ProblemDefinition.h"
 %template (ProblemDefinitioncomplexdouble)  tdkp::ProblemDefinition<complex<double> , double, double>; 
 %template (ProblemDefinitiondoubledouble)   tdkp::ProblemDefinition<double, double, double>; 
 %template (ProblemDefinitioncomplexcomplex) tdkp::ProblemDefinition<complex<double> , complex<double>, double>;

 %include "tdkp/probdefs/LinearProblem.h"
 %template (LinearProblemDouble) tdkp::LinearProblem<double>; 
 %template (LinearProblemcomplex) tdkp::LinearProblem<complex<double> >; 
 %include "tdkp/probdefs/IntrinsicStrain.h" 
 %feature("notabstract") PoissonEquation;
 %include "tdkp/probdefs/PoissonEquation.h"
 %include "tdkp/probdefs/EigenProblem.h"
 %template (EigenProblemcomplexdouble)  tdkp::EigenProblem<complex<double> , double, double>;
 %template (EigenProblemcomplexcomplex) tdkp::EigenProblem<complex<double> , complex<double>, double >;

 %include "tdkp/probdefs/EigenProblem3D.h" 
 %template (EigenProblem3Dcomplexdouble)  tdkp::EigenProblem3D<complex<double> , double>;
 %template (EigenProblem3Dcomplexcomplex) tdkp::EigenProblem3D<complex<double> , complex<double>  >;

 %template (SchroedingerProblem3Dcomplexcomplex) tdkp::SchroedingerProblem<EigenProblem3D<complex<double>, complex<double> > >;
 %template (SchroedingerProblem1D2Dcomplexcomplex) tdkp::SchroedingerProblem<EigenProblem<complex<double>, complex<double>, double> >;
 %template (SchroedingerProblemcomplexdouble)  tdkp::SchroedingerProblem<EigenProblem3D<complex<double>, double> >;

 %include "tdkp/probdefs/EffectiveMass.h"
 %include "tdkp/probdefs/NEGFEffectiveMass.h"

 %include "tdkp/probdefs/KPBase.h" 
 %template (KPBase3DParentcomplex)    tdkp::KPBase<tdkp::SchroedingerProblem<tdkp::EigenProblem3D<complex<double>, complex<double> > > >;
 %template (KPBase1D2DParentcomplex)  tdkp::KPBase<tdkp::SchroedingerProblem<tdkp::EigenProblem<complex<double>, complex<double>, double> > >;
 %include "tdkp/probdefs/KPBase3D.h" 
 %include "tdkp/probdefs/KPBase1D2D.h"
 %feature("notabstract") KP4x43D;
 %feature("notabstract") KP6x63D;
 %feature("notabstract") KP6x63DWZ;
 %include "tdkp/probdefs/KP4x43D.h"
 %feature("notabstract") KP8x83D;
 %feature("notabstract") KP8x83DWZ;
 %include "tdkp/probdefs/KP8x83D.h"
 %feature("notabstract") KP4x41D2D;
 %include "tdkp/probdefs/KP4x41D2D.h"
 %feature("notabstract") KP6x61D2D;
 %feature("notabstract") KP6x61D2DWZ;
 %include "tdkp/probdefs/KP6x61D2D.h"
 %feature("notabstract") KP8x81D2D;
 %feature("notabstract") KP8x81D2DWZ;
 %include "tdkp/probdefs/KP8x81D2D.h"

 %include "tdkp/probdefs/SchroedingerPML.h"
 %template (EffectiveMassPMLparent) tdkp::SchroedingerPML<EffectiveMass>;
 %template (EffectiveMassPML)       tdkp::SchroedingerPML3D<EffectiveMass>;
 %template (KP4x43DPMLparent)       tdkp::SchroedingerPML<KP4x43D>;
 %template (KP4x43DPML)             tdkp::SchroedingerPML3D<KP4x43D>;
 %template (KP6x63DPMLparent)       tdkp::SchroedingerPML<KP6x63D>;
 %template (KP6x63DPML)             tdkp::SchroedingerPML3D<KP6x63D>;
 %template (KP8x83DPMLparent)       tdkp::SchroedingerPML<KP8x83D>;
 %template (KP8x83DPML)             tdkp::SchroedingerPML3D<KP8x83D>;  
 %template (KP6x63DWZPMLparent)     tdkp::SchroedingerPML<KP6x63DWZ>;
 %template (KP6x63DWZPML)           tdkp::SchroedingerPML3D<KP6x63DWZ>;
 %template (KP8x83DWZPMLparent)     tdkp::SchroedingerPML<KP8x83DWZ>;
 %template (KP8x83DWZPML)           tdkp::SchroedingerPML3D<KP8x83DWZ>;   
 %template (KP4x41D2DPMLparent)     tdkp::SchroedingerPML<KP4x41D2D>;
 %template (KP4x41D2DPML)           tdkp::SchroedingerPML1D2D<KP4x41D2D>;
 %template (KP6x61D2DPMLparent)     tdkp::SchroedingerPML<KP6x61D2D>;
 %template (KP6x61D2DPML)           tdkp::SchroedingerPML1D2D<KP6x61D2D>;
 %template (KP8x81D2DPMLparent)     tdkp::SchroedingerPML<KP8x81D2D>;
 %template (KP8x81D2DPML)           tdkp::SchroedingerPML1D2D<KP8x81D2D>;
 %template (KP6x61D2DWZPMLparent)   tdkp::SchroedingerPML<KP6x61D2DWZ>;
 %template (KP6x61D2DWZPML)         tdkp::SchroedingerPML1D2D<KP6x61D2DWZ>;
 %template (KP8x81D2DWZPMLparent)   tdkp::SchroedingerPML<KP8x81D2DWZ>;
 %template (KP8x81D2DWZPML)         tdkp::SchroedingerPML1D2D<KP8x81D2DWZ>;
 
 %include "tdkp/utilities/SchroedingerPoisson.h"
 %template (SPKP6x61D2DWZ)			tdkp::SchroedingerPoisson<EffectiveMass,KP6x61D2DWZ>;
 %template (SPKP8x81D2D)			tdkp::SchroedingerPoisson<KP8x81D2D,KP8x81D2D>;
 %template (SPKP8x83D)				tdkp::SchroedingerPoisson<KP8x83D,KP8x83D>;

 %feature("notabstract") MatrixElements;
 %feature("notabstract") MomentumOperator4x4;
 %feature("notabstract") MomentumOperator6x6;
 %feature("notabstract") MomentumOperator8x8;
 %feature("notabstract") MomentumOperator6x6WZ;
 %feature("notabstract") MomentumOperator8x8WZ;
 %include "tdkp/utilities/MatrixElements.h"


 %include "tdkp/utilities/DensityOfStates.h"
 %include "tdkp/utilities/SLC.h"
 %include "tdkp/utilities/Fermi.h"
 %include "tdkp/coulomb/CoulombMatrixElement.h" 
 %feature("notabstract") Overlap;
 %include "tdkp/utilities/Overlap.h"
 %include "tdkp/utilities/PiezoElectricTensor.h"
 %include "tdkp/probdefs/PolarizationCharge.h"
 
 %include "tdkp/utilities/GridInterpolator.h"
 %extend tdkp::GridInterpolator {
	 %template (interpolate_double) interpolate<double>;
 }

 %include "tdkp/clc/CLCIngredients.h"
 %include "tdkp/clc/CLCOpticalResults.h"
 %include "tdkp/clc/CLCSHF.h"
 %include "tdkp/clc/CLCScreening.h"
 
 
 
 