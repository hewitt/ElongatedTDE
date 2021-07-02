/// IN PROGRESS
///
/// switches: -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps
///
#include <DenseVector.h>
#include <OneD_Node_Mesh.h>
#include <TwoD_Node_Mesh.h>
#include <SparseLinearSystem.h>
#include <Types.h>
#include <PetscSession.h>
#include <Utility.h>
#include <Timer.h>
#include <TrackerFile.h>
#include <cmath>

//#define GSDC3     // the GSDC3 roughness shape
#define BASIC     // a simple test roughness
//#define SYKES     // a 2D test compact roughness

//#define TWOD      // 2D calc only, ignores m-th harmonics -- best set NZ=1

#define NONLINEAR // otherwise linearised solution found with no harmonic interaction
//#define AFROMP    // treatment of the Hilbert integral either A'=int of P or P=int of A'


enum {V,U,W}; 
 
namespace CppNoddy
{ 
  namespace Example
  {

#ifdef SYKES
    double h(3.0);

    const double L(1.0); // doesn't matter
    
    double F( double X ) {
      if (abs(X)<1) {
        return pow((1-X*X),2);
      }
      return 0.0;
    }

    // 2D
    double Fm( double X, int m ) {
      if ( m == 0 ) {        
        return 2*0.5*F(X);
      } 
      return 0.0;
    }
#endif    

    
#ifdef BASIC
    double h(3.0);
    const double L(1/(2*M_PI));
    const double d0(2./3.);
    
    double F( double X ) {
      return exp(-X*X/d0);  
    }

    double Fm( double X, int m ) {
      if ( m == 0 ) {        
        return 0.5*F(X);
      } else {
        if ( m == 1 ) {
          return 0.25*F(X);        
        }
      }
      return 0.0;
    }
#endif    

#ifdef GSDC3
    double h( 2.0 );
    const double L( 2./M_PI );
    DenseVector<D_complex> Fms;

    double F( double X ) {
      const double d0(0.5);
      return 1.0/(1+pow(X/d0,2));
    }

    double Fhat( double Z ) {
      const double d(1.0);
      const double kappa(2.0);
      if ( abs(Z) < 1 ) {
        return exp( -(Z*Z/(d*d))*(kappa*kappa + d*d/(d*d-Z*Z)) );
      } else {
        return 0.0;
      }        
    }

    double Fm( double X, int m ) {
      return Fms[m].real()*F(X);
    }

#endif   
    
  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using std::cout;

int main(int argc, char *argv[]) {
  PetscSession::getInstance(argc,argv);

  const D_complex eye(0.0,1.0);
  double relax(0.3);

  double Xup = -40.0;
  double Xdown = 40.0;
  double Ymax = 30.0;
  // number of (spatial) nodal points in y
  unsigned NY = 201;
  unsigned NX = 801;
  unsigned NZ = 8;  // [-NZ:+NZ] total Fourier modes in Z
  //
  unsigned NE0 = 2; // number of extra d.o.f per X node for the mean flow: P & A
  unsigned NV0 = 2; // number of d.o.f at each (X,Y) node for the mean flow: V_0,U_0
  unsigned NEm = 2; // number of extra d.o.f per X node for the m-th harmonic: Q_m,Q^\dagger_m
  unsigned NVm = 3; // number of d.o.f at each (X,Y) node for the m-th harmonic: V_m,U_m,W_m

  // GSDC roughness
  // Run A : -30:30 x 0:30 with 601x201 and NZ=5
  // Run B : -30:30 x 0:30 with 801x201 and NZ=5
  // Run C : -30:30 x 0:45 with 601x301 and NZ=5
  // Run D : -30:50 x 0:30 with 801x201 and NZ=5
  // Run E : -30:20 x 0:30 with 601x201 and NZ=10
  // Run F : -30:20 x 0:30 with 801x201 and NZ=15
  // Run G : -30:80 x 0:30 with 2001x101 and NZ=10
  // Run H : -30:80 x 0:30 with 4001x101 and NZ=15
  // Run I : -30:80 x 0:60 with 2001x201 and NZ=10
  
  // BASIC1 : h=2, -30:30 x 0:30 with 601x201 and NZ=4
  // BASIC2 : h=2, -30:30 x 0:30 with 601x201 and NZ=8
  // BASIC3 : h=2, -30:100 x 0:30 with 1201x201 and NZ=4
  // BASIC4 : h=2, -30:100 x 0:30 with 1201x201 and NZ=8
  // BASIC5 : h=2, -30:100 x 0:30 with 1201x201 and NZ=12
  // BASIC6 : h=2, -30:30 x 0:45 with 601x301 and NZ=8
  // BASIC7 : h=2, -30:100 x 0:45 with 1201x301 and NZ=8  
  
  cout << "\n";
  cout << "=== GSDC elongated triple deck problem for hump =====\n";
  cout << "\n";
  cout << "|| Domain: " << Xup << " to " << Xdown << " and 0 to " << Ymax << "\n";
  cout << "|| Nodes: " << NX << " by " << NY << " nodes.\n";
  cout << "|| Fourier in Z: " << NZ << "\n";

  // streamwise noodes
  DenseVector<double> Xnodes = Utility::uniform_node_vector( Xup, Xdown, NX );
  const double dX = Xnodes[1]-Xnodes[0];
  cout << "|| -----------------\n";
  cout << "|| dX=" << dX << "\n";
  // transverse nodes
  DenseVector<double> Ynodes = Utility::uniform_node_vector( 0, Ymax, NY );
  const double dY = Ynodes[1]-Ynodes[0];
  cout << "|| dY=" << dY << "\n";
  const double dY2 = dY*dY;

  cout << "|| -----------------\n";
#ifdef NONLINEAR
  cout << "|| NONLINEAR computation *including* R_m terms\n";
#else
  cout << "|| This is only a LINEAR computation, ignoring the R_m forcing\n";
#endif
  cout << "|| -----------------\n";
  
  // mean flow
  TwoD_Node_Mesh<D_complex> meanSoln( Xnodes, Ynodes, NV0 ); // V_0,U_0 at each X,Y
  OneD_Node_Mesh<D_complex> meanPA( Xnodes, 3 );             // P,A (and tau_0) at each X
  // harmonics are stored as a STL vector of pointers to 2D meshes
  std::vector< TwoD_Node_Mesh<D_complex>* > harmonicsSoln;   // storage for Fourier (Z) modes (vels)
  std::vector< OneD_Node_Mesh<D_complex>* > harmonicsQ;      // storage for Fourier (Z) modes (pressure)
  // create the containers in the STL vector
  for ( unsigned m=1; m<=NZ; ++m ){
    TwoD_Node_Mesh<D_complex>* p_mthSoln = new TwoD_Node_Mesh<D_complex>( Xnodes, Ynodes, NVm ); // V_m,U_m,W_m at each X,Y
    OneD_Node_Mesh<D_complex>*p_mthQ = new OneD_Node_Mesh<D_complex>( Xnodes, 3 );               // Q_m, Q_m^\dagger (and tau_m) at each X
    harmonicsSoln.push_back(p_mthSoln);
    harmonicsQ.push_back(p_mthQ);
  }

  // metrics/timers
  Timer timer("Assemble the R0 terms");
  double lastIntA(0.0);
  double counter(0.0);
  TrackerFile conv("./DATA/conv.dat");
  conv.push_ptr( &counter, "counter" );
  conv.push_ptr( &lastIntA, "A integrated" );
  conv.header();
  
  // the +NE0*NX below is the P and A at each X station
  // could do this as "double" but need the Complex PETSc solver anyway for the harmonics.
  SparseMatrix<D_complex> A0( NX*(NY*NV0+NE0), NX*(NY*NV0+NE0) );
  DenseVector<D_complex> B0( NX*(NY*NV0+NE0), 0.0 );
  SparseLinearSystem<D_complex> meanSparseSystem( &A0, &B0, "petsc" );

  // the +NEm*NX below is the Q_m at each X station
  // Complex solution for the m-th harmonic
  SparseMatrix<D_complex> Am( NX*(NY*NVm+NEm), NX*(NY*NVm+NEm) );
  DenseVector<D_complex> Bm( NX*(NY*NVm+NEm), 0.0 );
  SparseLinearSystem<D_complex> mthSparseSystem( &Am, &Bm, "petsc" );

  int loopCounter(1);  
 top:
  cout << "\n|| Global loop counter = " << loopCounter << "\n";
  
  /////////////////////////////////////
  // Compute the mean flow solution  //
  /////////////////////////////////////
  double absResidual(1.0);
  cout << "|| -----------------\n";
  cout << "|| h = " << Example::h << "\n";
  cout << "|| -----------------\n\n";
  cout << "Mean flow, m=0 calc.\n";

 meanflow:
  //
  {
    const unsigned nextX = NV0*NY + NE0; // +NE0 for the P & A at each X station
    A0.blank();
    {
      // i = 0 : X = Xup
      unsigned i(0);
      unsigned O( (NV0*NY+NE0)*i ); // =0 obvs.
      for ( unsigned j=0; j<NY; ++j ) {
        // each interior transverse node
        unsigned k( NV0*j ); 
        // V = 0
        A0( O + k + 0, O + k + V ) = 1.0;
        B0[ O + k + 0 ] = -(meanSoln(i,j,V)-0.0);
        // U = 0
        A0( O + k + 1, O + k + U ) = 1.0;
        B0[ O + k + 1 ] = -(meanSoln(i,j,U)-0.0);
      }
      {
        unsigned j( NY-1 );
        unsigned k( NV0*j );
        // P = 0
        A0( O + k + 2, NY*NV0 + i*nextX ) = 1.0;
        B0[ O + k + 2 ] = -(meanPA(i,0)-0.0);
        // A = 0
        A0( O + k + 3, NY*NV0 + 1 + i*nextX ) = 1.0;
        B0[ O + k + 3 ] = -(meanPA(i,1)-0.0);
      }
    }
    // INTERIOR X
    for ( unsigned i=1; i<NX-1; ++i ) {
      // each interior downstream node
      unsigned O( (NV0*NY+NE0)*i ); // =0 obvs.
      {
        // boundary Y=0
        unsigned j(0);
        unsigned k( NV0*j ); // =0 obvs.
        // V = 0
        A0( O + k + 0, O + k + V ) = 1.0;
        B0[ O + k + 0 ] = -(meanSoln(i,j,V)-0.0);
        // U = 0
        A0( O + k + 1, O + k + U ) = 1.0;
        B0[ O + k + 1 ] = -(meanSoln(i,j,U)-0.0);
      }
      for ( unsigned j=1; j<NY-1; ++j ) {
        // each interior transverse node
        unsigned k( NV0*j ); 

        ////////////////////////////
        // ctty at mid point in Y //
        ////////////////////////////
        // V_Y at X_i, Y_{j-1/2}
        A0( O + k + 0, O + k + V ) = 1.0/dY;
        A0( O + k + 0, O + k + V - NV0 ) = -1.0/dY;
        // U_X at X_i, Y_{j-1/2}
        A0( O + k + 0, O + k + U + nextX ) = 0.5/(2*dX);
        A0( O + k + 0, O + k + U - NV0 + nextX ) = 0.5/(2*dX);
        A0( O + k + 0, O + k + U - nextX ) = -0.5/(2*dX);
        A0( O + k + 0, O + k + U - NV0 - nextX ) = -0.5/(2*dX);      
        // residual for the ctty equation
        B0[ O + k + 0 ] = -( (0.5*(meanSoln(i+1,j,U)+meanSoln(i+1,j-1,U))
                              - 0.5*(meanSoln(i-1,j,U)+meanSoln(i-1,j-1,U)))/(2*dX)
                            + (meanSoln(i,j,V)-meanSoln(i,j-1,V))/dY );

        //////////////////
        //   U mmt eqn  //
        //////////////////

        // R_0^u terms
        D_complex R0(0.0);
        // "Lin=0" ignores nonlinear interaction of mean flow component with itself
        int Lin(0);
#ifdef NONLINEAR
        Lin = 1;
        for ( unsigned n = 1; n <= NZ; ++n ) {
          TwoD_Node_Mesh<D_complex>* p_nthSoln = harmonicsSoln[n-1];
          R0 += p_nthSoln->operator()(i,j,U)*(p_nthSoln->operator()(i+1,j,U)
                                              -p_nthSoln->operator()(i-1,j,U))/(2*dX);
          R0 += p_nthSoln->operator()(i,j,V)*(p_nthSoln->operator()(i,j+1,U)
                                              -p_nthSoln->operator()(i,j-1,U))/(2*dY);
          R0 -= p_nthSoln->operator()(i,j,W)*eye*(n/Example::L)*p_nthSoln->operator()(i,j,U);
        }
        R0 += std::conj(R0);
#endif
        
        // U_{i,j+1}
        A0( O + k + 1, O + k + U + NV0 ) = 1.0/dY2 - Lin*meanSoln(i,j,V)/(2*dY);
        // U_{i,j}
        A0( O + k + 1, O + k + U ) = -2.0/dY2 - Lin*(meanSoln(i+1,j,U)-meanSoln(i-1,j,U))/(2*dX);
        // U_{i,j-1}
        A0( O + k + 1, O + k + U - NV0 ) = 1.0/dY2 + Lin*meanSoln(i,j,V)/(2*dY);
        // U_{i+1,j}
        A0( O + k + 1, O + k + U + nextX ) = -Ynodes[j]/(2*dX) - Lin*meanSoln(i,j,U)/(2*dX);
        // U_{i-1,j}
        A0( O + k + 1, O + k + U - nextX ) = +Ynodes[j]/(2*dX) + Lin*meanSoln(i,j,U)/(2*dX);
        // V_{i,j}
        A0( O + k + 1, O + k + V ) = -1.0 - Lin*(meanSoln(i,j+1,U)-meanSoln(i,j-1,U))/(2*dY);
        // CENTRAL
        // P_{i+1} 
        A0( O + k + 1, NV0*NY + (i+1)*nextX ) = -1.0/(2*dX);
        // P_{i-1}
        A0( O + k + 1, NV0*NY + (i-1)*nextX ) = 1.0/(2*dX);
        // BACKWARD
        // P_{i} 
        //A0( O + k + 1, NV0*NY + (i)*nextX ) = -1.0/(dX);
        // P_{i-1}
        //A0( O + k + 1, NV0*NY + (i-1)*nextX ) = 1.0/(dX);

        // residual for the (nonlinear) U-momentum equation
        B0[ O + k + 1 ] = -( (meanSoln(i,j+1,U)-2.0*meanSoln(i,j,U)+meanSoln(i,j-1,U))/dY2
                             -(meanPA(i+1,0)-meanPA(i-1,0))/(2.*dX)
                            -Ynodes[j]*(meanSoln(i+1,j,U)-meanSoln(i-1,j,U))/(2*dX)
                            -meanSoln(i,j,V)
                            -Lin*meanSoln(i,j,U)*(meanSoln(i+1,j,U)-meanSoln(i-1,j,U))/(2*dX)
                            -Lin*meanSoln(i,j,V)*(meanSoln(i,j+1,U)-meanSoln(i,j-1,U))/(2*dY) - R0);
      }
      {
        unsigned j(NY-1);
        unsigned k( NV0*j );
        // ctty at mid point in Y
        A0( O + k + 0, O + k + V ) = 1.0/dY;
        A0( O + k + 0, O + k + V - NV0 ) = -1.0/dY;
        A0( O + k + 0, O + k + U + nextX ) = 0.5/(2*dX);
        A0( O + k + 0, O + k + U - NV0 + nextX ) = 0.5/(2*dX);
        A0( O + k + 0, O + k + U - nextX ) = -0.5/(2*dX);
        A0( O + k + 0, O + k + U - NV0 - nextX ) = -0.5/(2*dX);      
        B0[ O + k + 0 ] = -( (0.5*(meanSoln(i+1,j,U)+meanSoln(i+1,j-1,U)) - 0.5*(meanSoln(i-1,j,U)+meanSoln(i-1,j-1,U)))/(2*dX)
                            + (meanSoln(i,j,V)-meanSoln(i,j-1,V))/dY );

        // Ud=0
        A0( O + k + 1, O + k + U ) = 3.0/(2*dY);
        A0( O + k + 1, O + k + U - NV0 ) = -4.0/(2*dY);
        A0( O + k + 1, O + k + U - 2*NV0 ) = +1.0/(2*dY);
        B0[ O + k + 1 ] = -( (3*meanSoln(i,j,U)-4*meanSoln(i,j-1,U)+meanSoln(i,j-2,U))/(2*dY) );

        // displacement condition
        // U
        A0( O + k + 2, O + k + U ) = 1.0;
        // -A
        A0( O + k + 2, NV0*NY + 1 + i*nextX ) = -1.0;
        // = h*F(x) + A - U 
        B0[ O + k + 2 ] = Example::h*Example::Fm(Xnodes[i],0) + meanPA(i,1) - meanSoln(i,j,U);

        // interaction condition from the Hilbert integral

        // RHS: pi*P - int (A'(s)-A'(X))/(X-s) ds - A'(X)int(1/(X-s))ds = 0
        B0[ O + k + 3 ] = M_PI*meanPA(i,0);
        D_complex AXm = (meanPA(i+1,1)-meanPA(i-1,1))/2.;
        D_complex sum1(0.0);
        D_complex sum2(0.0);
                for ( unsigned n=1; n<NX; ++n ) {
          // mid-point integration variable
          double sm = 0.5*(Xnodes[n]+Xnodes[n-1]);
          // A'(sm)*ds : mid point
          D_complex Asm = (meanPA(n,1)-meanPA(n-1,1));
          // contribution to the RHS from the integral: (A'(s)-A'(X))/(X-s) @ mid-point
          sum1 += Asm/(Xnodes[i]-sm);
          sum2 += AXm/(Xnodes[i]-sm);
          //B0[ O + k + 3 ] -= (Asm-AXm)/(Xnodes[i]-sm);
        }
        B0[ O + k + 3 ] -= sum1;
        B0[ O + k + 3 ] += sum2;
        //if ( i == 500 ) {
        //  cout << "X = " << Xnodes[i] << " sum2 = " << sum2
        //       << " analytic = " << -(AXm/dX)*log( abs(Xnodes[i]-Xup)/abs(Xnodes[i]-Xdown) ) <<"\n";
        //}
        B0[ O + k + 3 ] -= (AXm/dX)*log( abs(Xnodes[i]-Xup)/abs(Xnodes[i]-Xdown) );

        // Jacobian terms NOT from the integral
        // P_i
        A0( O + k + 3, NV0*NY + 0 + i*nextX ) = -M_PI;
        // A_i
        A0( O + k + 3, NV0*NY + 1 + (i+1)*nextX ) = +log( abs(Xnodes[i]-Xup)/abs(Xnodes[i]-Xdown) )/(2.*dX);
        // A_{i-1}
        A0( O + k + 3, NV0*NY + 1 + (i-1)*nextX ) = -log( abs(Xnodes[i]-Xup)/abs(Xnodes[i]-Xdown) )/(2*dX);
        // Jacobian terms from the integral
        for ( unsigned n=1; n<NX; ++n ) {
          // mid point in the integration (dummy) variable
          double sm = 0.5*(Xnodes[n]+Xnodes[n-1]);
          // mid-point evaluation of the A'(s) function
          // A_n
          A0( O + k + 3, NV0*NY + 1 + n*nextX ) += 1.0/(Xnodes[i]-sm);
          // A_{n-1}
          A0( O + k + 3, NV0*NY + 1 + (n-1)*nextX ) += -1.0/(Xnodes[i]-sm);
          // A_i
          A0( O + k + 3, NV0*NY + 1 + (i+1)*nextX ) += -0.5/(Xnodes[i]-sm);
          // A_{i-1}
          A0( O + k + 3, NV0*NY + 1 + (i-1)*nextX ) += 0.5/(Xnodes[i]-sm);
        }
        
        // part of integral from Xdown to infty
        {
          double Ahat = meanPA(NX-1,1).real()/pow(Xnodes[NX-1],1./3.); // constant in far field
          double sum(0.0);
          const double ds(dX);
          double sm( Xnodes[NX-1] + .5*ds );
          double term(0.0);
          do {
            term = ds*pow(sm,-2./3.)/(Xnodes[i]-sm);
            sum += term;
            sm += ds;
          } while ( abs(term) > 1.e-8 );
          sum *= (1./3.)*Ahat;
          //if ( i == 500 ) {           
          //  cout << "X= " << Xnodes[i] << " Smax = " << sm << " sum = " << sum
          //       << " Ahat*10^3 = " << 1.e3*Ahat <<"\n";
	  // }
          B0[ O + k + 3 ] -= sum;
        }        
      }
    }
    
    {
      // i = NX-1 : X = Xdown
      unsigned i(NX-1);
      unsigned O( (NV0*NY+NE0)*i ); // last X station
      for ( unsigned j=0; j<NY; ++j ) {
        // each interior transverse node
        unsigned k( NV0*j ); 
        // XV_0X + V_0/3 = 0 => V0 \sim X^{-1/3}
        A0( O + k + 0, O + k + V ) = Xnodes[i]*3./(2*dX) + 1.0/3.0;
        A0( O + k + 0, O + k + V - nextX ) = -Xnodes[i]*4./(2*dX);
        A0( O + k + 0, O + k + V - 2*nextX ) = Xnodes[i]*1./(2*dX);
        B0[ O + k + 0 ] = -( Xnodes[i]*(3.*meanSoln(i,j,V)-4.*meanSoln(i-1,j,V)+meanSoln(i-2,j,V))/(2*dX)
                             + meanSoln(i,j,V)/3.);
        // XU_0X - U_0/3 = 0 => U0 \sim X^{1/3} 
        A0( O + k + 1, O + k + U ) = Xnodes[i]*3./(2*dX) - 1.0/3.0;
        A0( O + k + 1, O + k + U - nextX ) = -Xnodes[i]*4./(2*dX);
        A0( O + k + 1, O + k + U - 2*nextX ) = Xnodes[i]*1./(2*dX);
        B0[ O + k + 1 ] = -( Xnodes[i]*(3.*meanSoln(i,j,U)-4.*meanSoln(i-1,j,U)+meanSoln(i-2,j,U))/(2*dX)
                             - meanSoln(i,j,U)/3.);
      }
      {
        unsigned j( NY-1 );
        unsigned k( NV0*j );
        // XP_X + 2P/3 = 0 => P \sim X{-2/3}
        A0( O + k + 2, NY*NV0 + 0 + i*nextX ) = Xnodes[i]*3.0/(2*dX) + 2./3.;
        A0( O + k + 2, NY*NV0 + 0 + (i-1)*nextX ) = -Xnodes[i]*4.0/(2*dX);
        A0( O + k + 2, NY*NV0 + 0 + (i-2)*nextX ) = Xnodes[i]*1.0/(2*dX);
        B0[ O + k + 2 ] = -( Xnodes[i]*(3*meanPA(i,0)-4*meanPA(i-1,0)+meanPA(i-2,0))/(2*dX)
                            + 2.*meanPA(i,0)/3. );
        // XA_X - A/3 = 0 => A \sim X^{1/3} 
        A0( O + k + 3, NY*NV0 + 1 + i*nextX ) = Xnodes[i]*3.0/(2*dX) - 1./3.;
        A0( O + k + 3, NY*NV0 + 1 + (i-1)*nextX ) = -Xnodes[i]*4.0/(2*dX);
        A0( O + k + 3, NY*NV0 + 1 + (i-2)*nextX ) = Xnodes[i]*1.0/(2*dX);
        B0[ O + k + 3 ] = -( Xnodes[i]*(3*meanPA(i,1)-4*meanPA(i-1,1)+meanPA(i-2,1))/(2*dX)
                            - meanPA(i,1)/3. );
      }
    }

    meanSparseSystem.solve();

    for ( unsigned i=0; i<NX; ++i ) {
      for ( unsigned j=0; j<NY; ++j ) {      
        // store V then U
        meanSoln(i,j,V) += B0[ V + j*NV0 + i*nextX ]; // V
        meanSoln(i,j,U) += B0[ U + j*NV0 + i*nextX ]; // U
      }
      meanPA(i,0) += B0[ NV0*NY + i*nextX ]; // P
      meanPA(i,1) += B0[ NV0*NY+1 + i*nextX ]; // A
      meanPA(i,2) = (-3.*meanSoln(i,0,U)+4.*meanSoln(i,1,U)-meanSoln(i,2,U)).real()/(2*dY); // tau_{Y=0}
    }
    absResidual = B0.inf_norm();
    cout << "residual = " << absResidual << "\n";
    cout << "intA difference = " << std::abs( std::abs(meanPA.integral2(1)) - lastIntA ) << "\n";

  } // one iteration of the mean flow m=0 equations has been done

  if (absResidual > 1.e-8) {
    goto meanflow;
  }

  //////////////////
  // IF CONVERGED //
  //////////////////
  if ( std::abs( std::abs(meanPA.integral2(1)) - lastIntA ) < 1.e-8 ) {
  //if ( absResidual < 1.e-8 ) {
    cout << "Converged.\n";

    meanSoln.dump_gnu("./DATA/meanSoln_"+Utility::stringify(Example::h,3)+".dat");
    meanPA.dump_gnu("./DATA/meanPA_"+Utility::stringify(Example::h,3)+".dat");
    for ( unsigned n=1; n<=NZ; ++n ) {
      // [0] => n=1 (first non-mean Fourier mode)
      harmonicsSoln[n-1] -> dump_gnu("./DATA/"+Utility::stringify(n,2)
                                     + "Soln_"+Utility::stringify(Example::h,3)
                                     + ".dat");
      harmonicsQ[n-1] -> dump_gnu("./DATA/"+Utility::stringify(n,2)
                                  + "Q_"+Utility::stringify(Example::h,3)
                                  + ".dat");
    }

    {
      // sum shear and higher-order pressure at the centreline
      OneD_Node_Mesh<double> centrelineData( Xnodes, 2 );
      for (unsigned i=0; i<NX; ++i ) {
        centrelineData( i, 0 ) = (-3.*meanSoln(i,0,U)+4.*meanSoln(i,1,U)-meanSoln(i,2,U)).real()/(2*dY); // tau_{Y=0}(Z=0) 
        centrelineData( i, 1 ) = 0.0; // Q(Z=0)
        for (unsigned n=1; n<=NZ; ++n ) {
          // z = 0 so exp(i*n*z/l) = 1
          TwoD_Node_Mesh<D_complex>* mthSoln = harmonicsSoln[n-1];
          centrelineData(i,0) += 2*(-3.*(*mthSoln)(i,0,U)+4.*(*mthSoln)(i,1,U)-(*mthSoln)(i,2,U)).real()/(2*dY); // tau_{Y=0}(Z=0)          
          centrelineData(i,1) += 2.*(harmonicsQ[n-1] -> operator()(i,0).real());    // Q(Z=0)       
        } // Z step
      } // X step
      centrelineData.dump_gnu("./DATA/centreline_"+Utility::stringify(Example::h,3)+".dat");
    }

    {
      // using M points in the Z direction
      unsigned M( 101 );
      DenseVector<double> Znodes( Utility::uniform_node_vector(-M_PI*Example::L, +M_PI*Example::L, M) );
      // contours of U_Y and roughness on the plane Y=0
      TwoD_Node_Mesh<double> contours( Xnodes, Znodes, 3 );
      for ( unsigned i=0; i<NX; ++i ) {
        for ( unsigned j=0; j<M; ++j ) {
          contours( i,j, 0 ) = Example::h*Example::Fm(Xnodes[i],0); // mean shear contribution to roughness
          contours( i,j, 1 ) = meanPA(i,2).real(); // mean shear contribution to shear
          contours( i,j, 2 ) = 0.0; // second-order pressure distribution
          for ( unsigned m=1; m<=NZ; ++m ) {
            contours( i, j, 0 ) += 2.*Example::h*Example::Fm(Xnodes[i],m)*exp(eye*m*Znodes[j]/Example::L).real(); // +/- m-th harmonic contribution to roughness
            contours( i, j, 1 ) += 2.*(
                                       (harmonicsQ[m-1] -> operator()(i,2))*exp(eye*m*Znodes[j]/Example::L)
                                       ).real(); // +/- m-th harmonic contribution to shear
            contours( i, j, 2 ) += 2.*(
                                       (harmonicsQ[m-1] -> operator()(i,0))*exp(eye*m*Znodes[j]/Example::L)
                                       ).real(); // +/- m-th harmonic contribution to 2nd order pressure Q
          } // Fourier modes
        } // Z nodes
      } // X nodes
      contours.dump_gnu("./DATA/contours_"+Utility::stringify(Example::h,3)+".dat");

      // TwoD_Node_Mesh<double> linearContours( "./DATA/contours_1.dat", NX, M, 3 );
      // for ( unsigned i=0; i<NX; ++i ) {
      //   for ( unsigned j=0; j<M; ++j ) {
      //     for ( unsigned k=0; k<3; ++k ) {
      //       contours(i,j,k)-= linearContours(i,j,k)*Example::h;
      //     }
      //   }
      // }
      // contours.dump_gnu("./DATA/contours_sub_"+Utility::stringify(Example::h,3)+".dat");

      OneD_Node_Mesh<double> errorMesh( Xnodes, 2 );
      // error check using skew-reciprocal relation between A'(X) and P(X)
      // so pi*A'(x) = - p.v. int P(s)/(X-s) ds
      errorMesh(0,0)=0.0; // A' from code
      errorMesh(0,1)=0.0; // A' from integral of P
      for (unsigned i=1; i<NX-1; ++i) {
        errorMesh(i,0)=(meanPA(i+1,1).real()-meanPA(i-1,1).real())/(2*dX);
        for (unsigned j=1; j<NX; ++j) {
          // mid-point in integration variable s
          double sm = (Xnodes[j]+Xnodes[j-1])*.5;
          // P(sm); ie. mid-point for 2nd order
          double Pm = (meanPA(j,0).real()+meanPA(j-1,0).real())*.5 - meanPA(i,0).real();
          // step in the integration variable
          double ds = (Xnodes[j]-Xnodes[j-1]);
          errorMesh(i,1) += ds*Pm/(Xnodes[i]-sm);
        }
        // extracting the singularity below makes little difference on a uniform mesh
        errorMesh(i,1) += meanPA(i,0).real()*log( abs((Xnodes[i]-Xup)/(Xnodes[i]-Xdown)) );
        errorMesh(i,1) *= (-1.0/M_PI);
      }
      errorMesh.dump_gnu("./DATA/error.dat");

    }
    
    assert(false); // exit
    
    //gamma = std::min( 1.0, gamma+0.25 );
    Example::h += 0.5; // increment roughness height
    loopCounter = 1;
    goto top;
  } else {
    cout << "Not yet converged.\n";
  }
  /// finisehd IF converged section
   
  lastIntA = std::abs(meanPA.integral2(1));
  counter = loopCounter;
  conv.update();

  ////////////////////////////////////////
  // Compute the m-th harmonic solution //
  ////////////////////////////////////////
  for ( unsigned mNum = 1; mNum <= NZ; ++mNum ) {
    cout << "Mode number m = " << mNum << " calc.\n";
    absResidual = 1.0;
    TwoD_Node_Mesh<D_complex>* mthSoln;
    OneD_Node_Mesh<D_complex>* mthQ;
    mthSoln = harmonicsSoln[ mNum - 1 ]; // indexed from zero for m=1
    mthQ = harmonicsQ[ mNum - 1 ]; // indexed from zero for m=1
    //
    { // tackle the m-th harmonic
      Am.blank();
      const unsigned nextX = NVm*NY + NEm; // +NEm for the Q_m, Q_m^\dagger at each X station
      {
        // i = 0 : X = Xup
        unsigned i(0);
        unsigned O( (NVm*NY+NEm)*i ); // =0 obvs.
        for ( unsigned j=0; j<NY; ++j ) {
          // each interior transverse node
          unsigned k( NVm*j ); 
          // V = 0
          Am( O + k + 0, O + k + V ) = 1.0;
          Bm[ O + k + 0 ] = -((*mthSoln)(i,j,V)-0.0);
          // U = 0
          Am( O + k + 1, O + k + U ) = 1.0;
          Bm[ O + k + 1 ] = -((*mthSoln)(i,j,U)-0.0);
          // W = 0
          Am( O + k + 2, O + k + W ) = 1.0;
          Bm[ O + k + 2 ] = -((*mthSoln)(i,j,W)-0.0);        
        }
        {
          unsigned j( NY-1 );
          unsigned k( NVm*j );
          // Q = 0
          Am( O + k + 3, NY*NVm + i*nextX ) = 1.0;
          Bm[ O + k + 3 ] = -((*mthQ)(i,0)-0.0);
          // Q^\dagger = 0
          Am( O + k + 4, NY*NVm + 1 + i*nextX ) = 1.0;
          Bm[ O + k + 4 ] = -((*mthQ)(i,1)-0.0);          
        }
      }
      // INTERIOR X
      for ( unsigned i=1; i<NX-1; ++i ) {
        // each interior downstream node
        unsigned O( (NVm*NY+NEm)*i ); // =0 obvs.
        {
          // boundary Y=0
          unsigned j(0);
          unsigned k( NVm*j ); // =0 obvs.
          // V = 0
          Am( O + k + 0, O + k + V ) = 1.0;
          Bm[ O + k + 0 ] = -((*mthSoln)(i,j,V)-0.0);
          // U = 0
          Am( O + k + 1, O + k + U ) = 1.0;        
          Bm[ O + k + 1 ] = -((*mthSoln)(i,j,U)-0.0);
          // W = 0
          Am( O + k + 2, O + k + W ) = 1.0;
          Bm[ O + k + 2 ] = -((*mthSoln)(i,j,W)-0.0);
        }
        for ( unsigned j=1; j<NY-1; ++j ) {
          // each interior transverse node
          unsigned k( NVm*j ); 

          ////////////////////////////
          // ctty at mid point in Y //
          ////////////////////////////
          // Vm_Y at X_i, Y_{j-1/2}
          Am( O + k + 0, O + k + V ) = 1.0/dY;
          Am( O + k + 0, O + k + V - NVm ) = -1.0/dY;
          // Um_X at X_i, Y_{j-1/2}
          Am( O + k + 0, O + k + U + nextX ) = 0.5/(2*dX);
          Am( O + k + 0, O + k + U - NVm + nextX ) = 0.5/(2*dX);
          Am( O + k + 0, O + k + U - nextX ) = -0.5/(2*dX);
          Am( O + k + 0, O + k + U - NVm - nextX ) = -0.5/(2*dX);
          // Wm_Z
          Am( O + k + 0, O + k + W ) = 0.5*(eye*mNum/Example::L);
          Am( O + k + 0, O + k + W - NVm ) = 0.5*(eye*mNum/Example::L);
          // residual for the ctty equation
          Bm[ O + k + 0 ] = -( (0.5*((*mthSoln)(i+1,j,U)+(*mthSoln)(i+1,j-1,U))
                                - 0.5*((*mthSoln)(i-1,j,U)+(*mthSoln)(i-1,j-1,U)))/(2*dX)
                               + ((*mthSoln)(i,j,V)-(*mthSoln)(i,j-1,V))/dY
                               + (eye*mNum/Example::L)*0.5*((*mthSoln)(i,j,W)+(*mthSoln)(i,j-1,W))
                               );

          // R_m^u terms
          D_complex Rm(0.0);
          // "Lin=0" turns off the mean-flow*mth-harmonic terms under the linear assumption
          int Lin(0);
#ifdef NONLINEAR
          Lin = 1;
          if ( mNum >= 2) { // only arise if m=2,3,...
            for ( unsigned n = 1; n <= mNum-1; ++n ) {
              TwoD_Node_Mesh<D_complex>* p_nthSoln = harmonicsSoln[n -1];
              TwoD_Node_Mesh<D_complex>* p_mmnthSoln = harmonicsSoln[mNum-n -1];
              Rm += p_nthSoln->operator()(i,j,U)*(p_mmnthSoln->operator()(i+1,j,U)
                                                  -p_mmnthSoln->operator()(i-1,j,U))/(2*dX);            
              Rm += p_nthSoln->operator()(i,j,V)*(p_mmnthSoln->operator()(i,j+1,U)
                                                  -p_mmnthSoln->operator()(i,j-1,U))/(2*dY);
              Rm += p_nthSoln->operator()(i,j,W)*eye*((mNum-n)/Example::L)*p_mmnthSoln->operator()(i,j,U);
            }
          }
          if ( mNum < NZ ) { // mNum=1,2,..,NZ so this only applies if n!=NZ 
              for ( unsigned n = 1; n <= NZ-mNum; ++n ) {
                TwoD_Node_Mesh<D_complex>* p_nthSoln = harmonicsSoln[n -1];
                TwoD_Node_Mesh<D_complex>* p_mpnthSoln = harmonicsSoln[mNum+n -1];
                Rm += p_mpnthSoln->operator()(i,j,U)*conj((p_nthSoln->operator()(i+1,j,U)
                                                           -p_nthSoln->operator()(i-1,j,U))/(2*dX));
                Rm += p_mpnthSoln->operator()(i,j,V)*conj((p_nthSoln->operator()(i,j+1,U)
                                                           -p_nthSoln->operator()(i,j-1,U))/(2*dY));
                Rm += p_mpnthSoln->operator()(i,j,W)*(-eye)*(n/Example::L)*conj(p_nthSoln->operator()(i,j,U));
              }
              for ( unsigned n = 1; n <= NZ-mNum; ++n ) {
                TwoD_Node_Mesh<D_complex>* p_nthSoln = harmonicsSoln[n -1];
                TwoD_Node_Mesh<D_complex>* p_mpnthSoln = harmonicsSoln[mNum+n -1];
                Rm += conj(p_nthSoln->operator()(i,j,U))*((p_mpnthSoln->operator()(i+1,j,U)
                                                           -p_mpnthSoln->operator()(i-1,j,U))/(2*dX));
                Rm += conj(p_nthSoln->operator()(i,j,V))*((p_mpnthSoln->operator()(i,j+1,U)
                                                           -p_mpnthSoln->operator()(i,j-1,U))/(2*dY));
                Rm += conj(p_nthSoln->operator()(i,j,W))*eye*((mNum+n)/Example::L)*(p_mpnthSoln->operator()(i,j,U));
              }              
          }
#endif
          
          //////////////////
          //   U mmt eqn  //
          //////////////////
          // U_{i,j+1}
          Am( O + k + 1, O + k + U + NVm ) = 1.0/dY2 - Lin*meanSoln(i,j,V)/(2*dY);
          // U_{i,j}
          Am( O + k + 1, O + k + U ) = -2.0/dY2 - Lin*(meanSoln(i+1,j,U)-meanSoln(i-1,j,U))/(2*dX);
          // U_{i,j-1}
          Am( O + k + 1, O + k + U - NVm ) = 1.0/dY2 + Lin*meanSoln(i,j,V)/(2*dY);
          // U_{i+1,j}
          Am( O + k + 1, O + k + U + nextX ) = -(Ynodes[j]+Lin*meanSoln(i,j,U))/(2*dX);
          // U_{i-1,j}
          Am( O + k + 1, O + k + U - nextX ) = +(Ynodes[j]+Lin*meanSoln(i,j,U))/(2*dX);
          // V_{i,j}
          Am( O + k + 1, O + k + V ) = -(1.0+Lin*(meanSoln(i,j+1,U)-meanSoln(i,j-1,U))/(2*dY));
          // No U-mmt pressure for the m-th harmonic        
          // residual for the (nonlinear) U-momentum equation
          Bm[ O + k + 1 ] = -( ((*mthSoln)(i,j+1,U)-2.0*(*mthSoln)(i,j,U)+(*mthSoln)(i,j-1,U))/dY2
                               -(Ynodes[j]+Lin*meanSoln(i,j,U))*((*mthSoln)(i+1,j,U)-(*mthSoln)(i-1,j,U))/(2*dX)
                               -(*mthSoln)(i,j,V)*(1.0+Lin*(meanSoln(i,j+1,U)-meanSoln(i,j-1,U))/(2*dY))
                               -(*mthSoln)(i,j,U)*Lin*(meanSoln(i+1,j,U)-meanSoln(i-1,j,U))/(2*dX)
                               -Lin*meanSoln(i,j,V)*((*mthSoln)(i,j+1,U)-(*mthSoln)(i,j-1,U))/(2*dY)
                               -Rm
                               );


          // R_m^w
          Rm = 0.0;
          // "Lin=0" turns off the mean-flow*mth-harmonic terms under the linear assumption
          Lin = 0;
#ifdef NONLINEAR
          Lin = 1;
          if ( mNum >= 2) { // only arise if m=2,3,...
            for ( unsigned n = 1; n <= mNum-1; ++n ) {
              TwoD_Node_Mesh<D_complex>* p_nthSoln = harmonicsSoln[n -1];
              TwoD_Node_Mesh<D_complex>* p_mmnthSoln = harmonicsSoln[mNum-n -1];
              Rm += p_nthSoln->operator()(i,j,U)*(p_mmnthSoln->operator()(i+1,j,W)
                                                  -p_mmnthSoln->operator()(i-1,j,W))/(2*dX);            
              Rm += p_nthSoln->operator()(i,j,V)*(p_mmnthSoln->operator()(i,j+1,W)
                                                  -p_mmnthSoln->operator()(i,j-1,W))/(2*dY);
              Rm += p_nthSoln->operator()(i,j,W)*eye*((mNum-n)/Example::L)*p_mmnthSoln->operator()(i,j,W);
            }
          }
          if ( mNum < NZ ) { // mNum=1,2,..,NZ so this only applies if n!=NZ 
              for ( unsigned n = 1; n <= NZ-mNum; ++n ) {
                TwoD_Node_Mesh<D_complex>* p_nthSoln = harmonicsSoln[n -1];
                TwoD_Node_Mesh<D_complex>* p_mpnthSoln = harmonicsSoln[mNum+n -1];
                Rm += p_mpnthSoln->operator()(i,j,U)*conj((p_nthSoln->operator()(i+1,j,W)
                                                           -p_nthSoln->operator()(i-1,j,W))/(2*dX));
                Rm += p_mpnthSoln->operator()(i,j,V)*conj((p_nthSoln->operator()(i,j+1,W)
                                                           -p_nthSoln->operator()(i,j-1,W))/(2*dY));
                Rm += p_mpnthSoln->operator()(i,j,W)*(-eye)*(n/Example::L)*conj(p_nthSoln->operator()(i,j,W));
              }
              for ( unsigned n = 1; n <= NZ-mNum; ++n ) {
                TwoD_Node_Mesh<D_complex>* p_nthSoln = harmonicsSoln[n -1];
                TwoD_Node_Mesh<D_complex>* p_mpnthSoln = harmonicsSoln[mNum+n -1];
                Rm += conj(p_nthSoln->operator()(i,j,U))*((p_mpnthSoln->operator()(i+1,j,W)
                                                           -p_mpnthSoln->operator()(i-1,j,W))/(2*dX));
                Rm += conj(p_nthSoln->operator()(i,j,V))*((p_mpnthSoln->operator()(i,j+1,W)
                                                           -p_mpnthSoln->operator()(i,j-1,W))/(2*dY));
                Rm += conj(p_nthSoln->operator()(i,j,W))*eye*((mNum+n)/Example::L)*(p_mpnthSoln->operator()(i,j,W));
              }
          }
#endif
          
          //////////////////
          //   W mmt eqn  //
          //////////////////
          // Wm_{i,j+1}
          Am( O + k + 2, O + k + W + NVm ) = 1.0/dY2 - Lin*meanSoln(i,j,V)/(2*dY);
          // Wm_{i,j}
          Am( O + k + 2, O + k + W ) = -2.0/dY2;
          // Wm_{i,j-1}
          Am( O + k + 2, O + k + W - NVm ) = 1.0/dY2 + Lin*meanSoln(i,j,V)/(2*dY);
          // Wm_{i+1,j}
          Am( O + k + 2, O + k + W + nextX ) = -(Ynodes[j]+Lin*meanSoln(i,j,U))/(2*dX);
          // Wm_{i-1,j}
          Am( O + k + 2, O + k + W - nextX ) = +(Ynodes[j]+Lin*meanSoln(i,j,U))/(2*dX);
          // Qm_i
          Am( O + k + 2, NVm*NY + i*nextX  ) = -eye*mNum/Example::L;
          // residual for the (nonlinear) W-momentum equation
          Bm[ O + k + 2 ] = -( ((*mthSoln)(i,j+1,W)-2.0*(*mthSoln)(i,j,W)+(*mthSoln)(i,j-1,W))/dY2
                               -(Ynodes[j]+Lin*meanSoln(i,j,U))*((*mthSoln)(i+1,j,W)-(*mthSoln)(i-1,j,W))/(2*dX)
                               -Lin*meanSoln(i,j,V)*((*mthSoln)(i,j+1,W)-(*mthSoln)(i,j-1,W))/(2*dY)
                               -(eye*mNum/Example::L)*(*mthQ)(i,0)
                               -Rm );
        }
        {
          unsigned j(NY-1);
          unsigned k( NVm*j );
          // ctty at mid point in Y
          Am( O + k + 0, O + k + V ) = 1.0/dY;
          Am( O + k + 0, O + k + V - NVm ) = -1.0/dY;
          Am( O + k + 0, O + k + U + nextX ) = 0.5/(2*dX);
          Am( O + k + 0, O + k + U - NVm + nextX ) = 0.5/(2*dX);
          Am( O + k + 0, O + k + U - nextX ) = -0.5/(2*dX);
          Am( O + k + 0, O + k + U - NVm - nextX ) = -0.5/(2*dX);
          Am( O + k + 0, O + k + W ) = 0.5*(eye*mNum/Example::L);
          Am( O + k + 0, O + k + W - NVm ) = 0.5*(eye*mNum/Example::L);        
          Bm[ O + k + 0 ] = -( (0.5*((*mthSoln)(i+1,j,U)+(*mthSoln)(i+1,j-1,U))
                                - 0.5*((*mthSoln)(i-1,j,U)+(*mthSoln)(i-1,j-1,U)))/(2*dX)
                               + ((*mthSoln)(i,j,V)-(*mthSoln)(i,j-1,V))/dY
                               + (eye*mNum/Example::L)*0.5*((*mthSoln)(i,j,W)+(*mthSoln)(i,j-1,W)) );

          // Um' = (i*mNum/L)^2 *Q^\dagger_m / Y^2
          Am( O + k + 1, O + k + U ) = 3.0/(2*dY);
          Am( O + k + 1, O + k + U - NVm ) = -4.0/(2*dY);
          Am( O + k + 1, O + k + U - 2*NVm ) = 1.0/(2*dY);
          Am( O + k + 1, NY*NVm + 1 + i*nextX ) = -std::pow(eye*mNum/Example::L,2)/pow(Ynodes[j],2);
          Bm[ O + k + 1 ] = -( (3*(*mthSoln)(i,j,U)-4*(*mthSoln)(i,j-1,U)+(*mthSoln)(i,j-2,U))/(2*dY)
                               - std::pow(eye*mNum/Example::L,2)*((*mthQ)(i,1))/pow(Ynodes[j],2) ); 
          // displacement condition
          // Um = h*Fm(x) - (i*mNum/L)^2 *Q^\dagger_m / Y
          Am( O + k + 2, O + k + U ) = 1.0;
          Am( O + k + 2, NY*NVm + 1 + i*nextX ) = +std::pow(eye*mNum/Example::L,2)/Ynodes[j];
          Bm[ O + k + 2 ] = Example::h*Example::Fm(Xnodes[i],mNum) - (*mthSoln)(i,j,U)
            - std::pow(eye*mNum/Example::L,2)*((*mthQ)(i,1))/Ynodes[j];
          // Wm = 0
          Am( O + k + 3, O + k + W ) = Ynodes[j]*3./(2*dY) + 1.0;
          Am( O + k + 3, O + k + W - NVm ) = -Ynodes[j]*4./(2*dY);        
          Am( O + k + 3, O + k + W - 2*NVm ) = Ynodes[j]*1./(2*dY);
          Bm[ O + k + 3 ] = -( Ynodes[j]*(3*(*mthSoln)(i,j,W)-4*(*mthSoln)(i,j-1,W)+(*mthSoln)(i,j-2,W))/(2*dY)
                               + (*mthSoln)(i,j,W) );
          // Q^\dagger equation
          Am( O + k + 4,  NY*NVm + 1 + (i-1)*nextX  ) = 1.0/(dX*dX);
          Am( O + k + 4,  NY*NVm + 1 + i*nextX ) = -2.0/(dX*dX);
          Am( O + k + 4,  NY*NVm + 0 + i*nextX ) = +1.0;
          Am( O + k + 4,  NY*NVm + 1 + (i+1)*nextX  ) = 1.0/(dX*dX);
          Bm[ O + k + 4 ] = -( ((*mthQ)(i-1,1)-2.*(*mthQ)(i,1)+(*mthQ)(i+1,1))/(dX*dX)
                               + (*mthQ)(i,0) );          
        }
      }
      {
        // i = NX-1 : X = Xdown
        unsigned i(NX-1);
        unsigned O( (NVm*NY+NEm)*i ); // =0 obvs.
        for ( unsigned j=0; j<NY; ++j ) {
          // each interior transverse node
          unsigned k( NVm*j ); 
          // XVm_X + V_m/3 = 0 => Vm \sim X^{-1/3} 
          Am( O + k + 0, O + k + V ) = Xnodes[i]*3./(2*dX) + 1.0/3.0;
          Am( O + k + 0, O + k + V - nextX ) = -Xnodes[i]*4./(2*dX);
          Am( O + k + 0, O + k + V - 2*nextX ) = Xnodes[i]*1./(2*dX);
          Bm[ O + k + 0 ] = -( Xnodes[i]*(3.*(*mthSoln)(i,j,V)-4.*(*mthSoln)(i-1,j,V)+(*mthSoln)(i-2,j,V))/(2*dX)
                               + (*mthSoln)(i,j,V)/3.);          
          // XUm_X - U_m/3 = 0 => Um \sim X^{1/3} 
          Am( O + k + 1, O + k + U ) = Xnodes[i]*3./(2*dX) - 1.0/3.0;
          Am( O + k + 1, O + k + U - nextX ) = -Xnodes[i]*4./(2*dX);
          Am( O + k + 1, O + k + U - 2*nextX ) = Xnodes[i]*1./(2*dX);
          Bm[ O + k + 1 ] = -( Xnodes[i]*(3.*(*mthSoln)(i,j,U)-4.*(*mthSoln)(i-1,j,U)+(*mthSoln)(i-2,j,U))/(2*dX)
                               - (*mthSoln)(i,j,U)/3.);
          // XWm_X + 2*W_mX/3 = 0 => Wm \sim X^{-2/3} 
          Am( O + k + 2, O + k + W ) = Xnodes[i]*3./(2*dX) + 2.0/3.0;
          Am( O + k + 2, O + k + W - nextX ) = -Xnodes[i]*4./(2*dX);
          Am( O + k + 2, O + k + W - 2*nextX ) = Xnodes[i]*1./(2*dX);
          Bm[ O + k + 2 ] = -( Xnodes[i]*(3.*(*mthSoln)(i,j,W)-4.*(*mthSoln)(i-1,j,W)+(*mthSoln)(i-2,j,W))/(2*dX)
                               + 2.*(*mthSoln)(i,j,W)/3.);
        }
        {
          unsigned j( NY-1 );
          unsigned k( NVm*j );
          // XQ_mX + 4*Q_mX/3 = 0 => Q_m \sim X^{-4/3} 
          Am( O + k + 3,  NY*NVm + i*nextX  ) = Xnodes[i]*3./(2*dX) + 4.0/3.0;
          Am( O + k + 3,  NY*NVm + (i-1)*nextX ) = -Xnodes[i]*4./(2*dX);
          Am( O + k + 3,  NY*NVm + (i-2)*nextX  ) = Xnodes[i]*1./(2*dX);
          Bm[ O + k + 3 ] = -( Xnodes[i]*(3.*(*mthQ)(i,0)-4.*(*mthQ)(i-1,0)+(*mthQ)(i-2,0))/(2*dX)
                               + 4.*(*mthQ)(i,0)/3.);
          // Q^\dagger_mX = 0 at X=Xupstream!
          Am( O + k + 4,  NY*NVm + 1 + 0*nextX  ) = -3.0/(2*dX);
          Am( O + k + 4,  NY*NVm + 1 + 1*nextX ) = +4.0/(2*dX);
          Am( O + k + 4,  NY*NVm + 1 + 2*nextX ) = -1.0/(2*dX);
          Bm[ O + k + 4 ] = -( (-3.*(*mthQ)(0,1)+4.*(*mthQ)(1,1)-(*mthQ)(2,1))/(2*dX) );          
        }
      }

      mthSparseSystem.solve();

#ifndef TWOD  // the THREE-DIMENSIONAL problem
      for ( unsigned i=0; i<NX; ++i ) {
        for ( unsigned j=0; j<NY; ++j ) {      
          (*mthSoln)(i,j,V) += relax*Bm[ V + j*NVm + i*nextX ]; // V
          (*mthSoln)(i,j,U) += relax*Bm[ U + j*NVm + i*nextX ]; // U
          (*mthSoln)(i,j,W) += relax*Bm[ W + j*NVm + i*nextX ]; // W
        }
        (*mthQ)(i,0) += relax*Bm[ NVm*NY + 0 + i*nextX ]; // Q_m
        (*mthQ)(i,1) += relax*Bm[ NVm*NY + 1 + i*nextX ]; // Q^\dagger_m
        (*mthQ)(i,2) = (-3.*(*mthSoln)(i,0,U)+4.*(*mthSoln)(i,1,U)-(*mthSoln)(i,2,U))/(2*dY); // tau_{Y=0}
      }
#else
      cout << "|| Not updating the harmonics since the TWOD flag is set\n";
#endif      
      absResidual = Bm.inf_norm();
      cout << "residual = " << absResidual << "\n";

    } // now we've completed an iteration of the m-th Fourier mode (in Z)
  }
  
  // go back to the start (m=0) for another round  
  loopCounter += 1;
  goto top;
  

}
