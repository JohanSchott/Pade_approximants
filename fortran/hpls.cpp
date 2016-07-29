#include <gmp.h>
#include <gmpxx.h>
#include <mblas_gmp.h>
#include <mlapack_gmp.h>
//#include <mutils_qd.h> 
#include <stdio.h>
#include <iostream>
#include <string>
#include <algorithm>
//#include <mblas_qd.h> 
//#include <mlapack_qd.h>
//#include <qd_complex.h>
#include <fstream>
#include <sstream>

// types for gmp:
#define INTEGER mpackint
#define COMPLEX mpc_class
#define REAL mpf_class
//#define cgesv17 cgesv17_




struct cmplx {
        double re;
        double im;
};

void ersetb(void);


extern "C"
{
void padec_(char *proj,int &contin,int &N,int &M,int &mtr,mp_bitcnt_t &nbit,cmplx * cof, char *stringname, int *nindex,double &err);
}

/*************************************************************************************
 *
 *      Subroutine to be called from Fortran to find pade coefficients.
 *	The subroutine will read input data from file *stringname and use the Matsubara
 *	frequencies specified by the vector *nindex. Precision will be set to nbit.
 *
 *      Input:
 *          integer N	        - number of Matsubara frequencies (should be even)
 *          integer M           - number of Padé coefficients (should be even)
 *          integer mtr    - number determining if inverting transpose multiplied matrix. Only use mtr=0 if N=M
 *      	integer*8  nbit     - precision with which to solve matrix system to get pade coefficients, WARNING-Integer*8!
 *          char * stringname   - file containing frequencies and corresponding selfenergies
 *		    nindex		        - vector specifying which lines to read in stringname
 *
 *      Output:
 *               cmplx*8 cof      - coifficients of Pade approximant
 *
 *      NOTE input data should be stored in three columns:
 *      	col 1: Matsubara frequencies
 *      	col 2: real part of self energy at matsubara points
 *      	col 3: imaginary part of self energy ar matsubara points
 *     
 *
 * **********************************************************************************/



void *__gxx_personality_v0;

void ersetb(void){
    setbuf(stdout,NULL); /* set output to unbuffered */
}

void padec_(char *proj,int &contin,int &N, int &M,int &mtr,mp_bitcnt_t &nbit, cmplx *cof,char *stringname, int *nindex,double &err){

    ersetb();

//Convert zmat and sigma from double complex to mpc_class
	if ((float)(N)/2 != (float)((N)/2)){
          std::cout << "Warning: not an even number of Matsubara points";
          return;
        }
	if ((float)(M)/2 != (float)((M)/2)){
          std::cout << "Warning: not an even number of pade coefficients";
          return;
        }
    if ((N<M)){
        std::cout << "Warning: N<M";
        return;
    }
	int r = M/2;

	mpf_t RE;  //Temporary variables for conversion from double to mpc_class            
    mpf_t IM;
	mpf_t rr;
	mpf_t rz;
	mpf_t mats;

    mpf_set_default_prec(nbit);  //Sets default precision and initializes floats
    mpf_init(RE);
    mpf_init(IM);
	mpf_init(rr);
	mpf_init(rz);
	mpf_init(mats);
   

	mpc_class *zmat_mp,*sigma_mp; //Mp variable arrays
	zmat_mp = new mpc_class[N];
	sigma_mp = new mpc_class[N];

   int base = 10;		     

	char uno [] = "1";   //To be used in matrix
	char zer [] = "0";
	mpf_set_str(rr,uno,base);
	mpf_set_str(rz,zer,base);
  




//--------------------------------------------------------------------




/*********************************************
 * RESDING in C++
 * ******************************************/
	char buff1[100];  
	char buff2[100];
	char buff3[100];


	std::string dfile = stringname;

	//std::ofstream outf;
	//outf.open("output/extra/buffcoff");
	int nind =0;
	int kindex=0;
   // temporary doubles
	double preal;
   double pimag; 

   if(nindex[0]<=0){ // Have negative Matsubara points
      int nneg=abs(nindex[0])+1;
      int npos=N-nneg;
      int mcount=0;
      nind=nneg;
      //read positive Matsubara (negative Matsubara are just copies, using G(-i*w)=G(i*w)^*
      std::ifstream in;
      in.open(dfile.c_str());
      if(in.is_open()){
         while(in && mcount < npos){
            in >> buff1;  //Matsubara frequency
            in >> buff2;  //Real part of function to be continued (self energy, greens function, hybridization,...)
            in >> buff3;  //Imaginary part of function to be continued
            //outf << buff1 << buff2 << buff3 << std::endl;
            kindex = kindex +1;
            if(kindex==nindex[nind]){
               mpf_set_str(mats,buff1,base);
               mpf_set_str(RE,buff2,base);
               mpf_set_str(IM,buff3,base);
               zmat_mp[nind]=mpc_class(rz,mats);  //complex matsubara values stored here
               sigma_mp[nind] = mpc_class(RE,IM); //complex selfenergy values stored here
               nind = nind +1;
               mcount=mcount+1;
            }
         }
      }
      in.close();
      //copy data in zmat_mp and sigma_mp with indicies nneg to 2*nneg-1
      for (mcount=0;mcount<nneg;mcount++) {
         zmat_mp[mcount]=-zmat_mp[2*nneg-1-mcount];
         sigma_mp[mcount]=2*sigma_mp[2*nneg-1-mcount].real()-sigma_mp[2*nneg-1-mcount]; //Is complex conjugate ok?
      }
   }
   else{
      std::ifstream in;
      in.open(dfile.c_str());
      if(in.is_open()){
         while(in && nind < N){
            in >> buff1;  //Matsubara frequency
            in >> buff2;  //Real part of function to be continued (self energy, greens function, hybridization,...)
            in >> buff3;  //Imaginary part of function to be continued
            //outf << buff1 << buff2 << buff3 << std::endl;
            kindex = kindex +1;
            if(kindex==nindex[nind]){
               mpf_set_str(mats,buff1,base);
               mpf_set_str(RE,buff2,base);
               mpf_set_str(IM,buff3,base);
               zmat_mp[nind]=mpc_class(rz,mats);  //complex matsubara values stored here
               sigma_mp[nind] = mpc_class(RE,IM); //complex selfenergy values stored here
               nind = nind +1;
            }
         }
      }
      in.close();
   }
   
   ////Print out the N values stored in zmat_mp and sigma_mp
   //std::ofstream kristoffer;
   //kristoffer.open("kristoffer.dat",std::ios_base::app);
   //kristoffer << "w_re   w_im   G_re   G_im\n";
   //kristoffer.flush();
   //for (int i=0;i<N;i++) {
   //   preal=zmat_mp[i].real().get_d();
   //   pimag=zmat_mp[i].imag().get_d();
   //   kristoffer << preal;
   //   kristoffer << "   ";
   //   kristoffer << pimag;
   //   kristoffer << "   ";
   //   preal=sigma_mp[i].real().get_d();
   //   pimag=sigma_mp[i].imag().get_d();
   //   kristoffer << preal;
   //   kristoffer << "   ";
   //   kristoffer << pimag;
   //   kristoffer << "\n";
   //}
   //kristoffer.flush();
   //kristoffer.close();


	//Create matrix to get coefficients (according to Beach matrix inversion algorithm)
	COMPLEX *X = new mpc_class[N*M];	
	COMPLEX *iw = new mpc_class[N];
    //what difference does complex do vs mpc_class ???
	
	for(int j=0; j<N;j++){
            X[j]=mpc_class(rr,rz);
		    iw[j]=mpc_class(rr,rz);
        }
	for(int k=1; k<r;k++){
                for (int j=0;j<N;j++){
                        iw[j] = iw[j]*zmat_mp[j];
                        X[k*N+j]=iw[j];
                        X[(r+k)*N+j]=-sigma_mp[j]*iw[j];
                }
        }
	for(int j=0; j<N;j++){
                X[r*N+j]= -sigma_mp[j];
        }
	
    

	//Construct righthandside of system X*x=p
	COMPLEX *p = new mpc_class[N];
	for(int j=0;j<N;j++){
      iw[j]=iw[j]*zmat_mp[j];
      p[j]=sigma_mp[j]*iw[j];
   }

    //reformulate LS-problem into a square-matrix problem (by multiplying with X^T)
    COMPLEX *A = new mpc_class[M*M];
    COMPLEX *b = new mpc_class[M];

    for(int i=0;i<M;i++){
            b[i] = mpc_class(rz,rz);
        for(int j=0;j<N;j++){
            b[i] = b[i] + X[i*N+j]*p[j];
        }
    }

    for(int i=0;i<M;i++){
        for(int j=0;j<M;j++){
            A[j*M+i] = mpc_class(rz,rz);
            for(int k=0;k<N;k++){
                A[j*M+i] = A[j*M+i] + X[i*N+k]*X[j*N+k];
            } 
        }
    }
    
      
   // Temporary copies
   COMPLEX *Xo = new mpc_class[N*M];
   COMPLEX *po = new mpc_class[N];
   COMPLEX *sol = new mpc_class[M];
   COMPLEX *errv = new mpc_class[N];
   for(int i=0;i<N;i++){
      for(int j=0;j<M;j++){
         Xo[i+N*j]=X[i+N*j];
      }
      po[i]=p[i];
   }

   if(mtr==1){ 
      // Use multiplied transposed matrix 
      //Call mpack routine to solve (M,M) problem

      int nrhs2 = 1;           //Solves for a single right handside
      mpackint info2;          //info vairable for call to Cgesv
      mpackint ipiv2[M];       //vector of pivot elements for call to Cgesv
      Cgesv(M,nrhs2,A,M,ipiv2,b,M,&info2);
      //--------------------------------------------------------------------------
      //Print coefficients to file
      std::ofstream myfile;
      std::stringstream ss;
      ss << contin;
      std::string str = ss.str();
	   std::string projstr = proj;
      str=projstr+"_cof_"+str;
      myfile.open(str.c_str());
      mp_exp_t exp;
      REAL *tmpr = new mpf_class[1];
      for(int i=0;i<M;i++){
         if(b[i].real().get_d()<0){
            tmpr[0]= -b[i].real();  
            myfile << "-0." << tmpr[0].get_str(exp) << "E" << exp;
         }else{
            myfile << "0." << b[i].real().get_str(exp) << "E" << exp;
         }
         if(b[i].imag().get_d()<0){
            tmpr[0]= -b[i].imag();  
            myfile << "      -0." << tmpr[0].get_str(exp) << "E" << exp << "\n";
         }else{
            myfile << "       0." << b[i].imag().get_str(exp) << "E" << exp << "\n";
         }
         //myfile << "0." << b[i].real().get_str(exp1) << "e" << exp1  << "     0." << b[i].imag().get_str(exp2) << "e" << exp2 << "\n";
      }
      myfile.close();
      
      //Translate coefficients to doubles to be sent back to Fortran
      for(int i=0;i<M;i++){
         sol[i]=b[i];
         preal = b[i].real().get_d();
         pimag = b[i].imag().get_d();
         cof[i].re = preal;
         cof[i].im = pimag;
      }

   }
   else if(mtr==0){
      if(N != M){ 
         std::cout << "Error: Can't use mpack and mtr=0 when N!=M";
         std::exit(1);
      }
      // Don't use multiplied transpose matrix
      //Call mpack routine to solve (N,M=N) problem

      int nrhs = 1;           //Solves for a single right handside
      mpackint info;          //info vairable for call to Cgesv
      mpackint ipiv[N];       //vector of pivot elements for call to Cgesv
      Cgesv(N,nrhs,X,N,ipiv,p,N,&info);

      //----------------------------------------------------------------------
      //Test how many digits are used in the mantissa.
         //mpf_t test5;
         //mpf_t test6;
         //mpf_t test7;
         //mpf_init_set_str(test5,"0",base);
         //mpf_init_set_str(test6,"10",base);
         //mpf_init_set_str(test7,"3",base);
         //mpf_div(test5,test6,test7);
         //zmat_mp[0]=mpc_class(rz,test5);  //complex matsubara values stored here
         //REAL *niklas = new mpf_class[1];
         //niklas[0]=zmat_mp[0].imag();
         //mp_exp_t exp;
         //
         //std::ofstream testf;
         //testf.open("coeff.dat");
         //testf << "0." << zmat_mp[0].imag().get_str(exp) << "e" << exp;
         //std::cout << "\n";
         //std::cout << "0.";
         //std::cout << niklas[0].get_str(exp);
         //std::cout << "e" << exp;
         //testf.close();
      //--------------------------------------------------------------------------
      //Print coefficients to file
      std::ofstream myfile;
      std::stringstream ss;
      ss << contin;
      std::string str = ss.str();
	   std::string projstr = proj;
      str=projstr+"_cof_"+str;
      myfile.open(str.c_str());
      mp_exp_t exp;
      REAL *tmpr = new mpf_class[1];
      for(int i=0;i<M;i++){
         if(p[i].real().get_d()<0){
            tmpr[0]= -p[i].real();  
            myfile << "-0." << tmpr[0].get_str(exp) << "E" << exp;
         }else{
            myfile << "0." << p[i].real().get_str(exp) << "E" << exp;
         }
         if(p[i].imag().get_d()<0){
            tmpr[0]= -p[i].imag();  
            myfile << "      -0." << tmpr[0].get_str(exp) << "E" << exp << "\n";
         }else{
            myfile << "       0." << p[i].imag().get_str(exp) << "E" << exp << "\n";
         }
         //myfile << "0."      << p[i].real().get_str(exp) << "e" << exp ;
         //myfile << "     0." << p[i].imag().get_str(exp) << "e" << exp << "\n";
      }
      myfile.close();

      //---------------------------------------------------------------------
      
      //Translate coefficients to doubles to be sent back to Fortran
      for(int i=0;i<M;i++){
         sol[i]=p[i];
         preal = p[i].real().get_d();
         pimag = p[i].imag().get_d();
         cof[i].re = preal;
         cof[i].im = pimag;
      }


   }
   else{
      std::cout << "Error: invrtqud neither 0 nor 1";
      std::exit(1);
   }
	
   // Calculate ||Xx-p||^2/||p||^2
   REAL *errp = new mpf_class[1];
   REAL *errnorm = new mpf_class[1];
   //double err;
   //errp[0]=mpf_class(rz); // same as line below 
   errp[0]=0;
   errnorm[0]=0;
   
   for(int i=0;i<N;i++){
      errv[i]=mpc_class(rz,rz);
      for(int j=0;j<M;j++){
         errv[i]=errv[i]+Xo[i+N*j]*sol[j];
      }
      errv[i]=errv[i]-po[i];
      errp[0]=errp[0]+errv[i].real()*errv[i].real()+errv[i].imag()*errv[i].imag();
   }
   for(int i=0;i<N;i++){
      errnorm[0]=errnorm[0]+po[i].real()*po[i].real()+po[i].imag()*po[i].imag();
   }

   // errp[0]=errp[0]/errnorm[0];
   // Want to take square root of errp[0], but don't know routine name so do it outside routine instead.
   err=errp[0].get_d();

   delete [] Xo;
   delete [] po;
   delete [] errv;
   delete [] sol;
   delete [] errp;
   delete [] errnorm;
        
	delete [] zmat_mp;
   delete [] sigma_mp;
   delete [] X;
   delete [] iw;
   delete [] p;
   delete [] A;
   delete [] b;


    
	return;
}	
