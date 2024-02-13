#include"udf.h"
/*Rotation velocity of the wall*/
#define OMEGA 2.222222

DEFINE_PROFILE(inlet_profile_velocity, thread, position)
 {
	 
	
	#if RP_NODE
	
		/*local vbariables to be used*/
		int n;
		face_t f;
		/* x[ND_ND] is the array where the centroids of the boundary are stored*/
		real x[ND_ND],RR, L, Etha;
		Node *ptrN;
		real Rmax=0.0;
		real Rmin=1.0e19;
		real RN;
		
	 begin_f_loop(f,thread)
		{
			f_node_loop(f,thread,n)
			{
				ptrN=F_NODE(f,thread,n);
				/*Radial coordinates of the nodes*/
				RN=sqrt(pow(NODE_X(ptrN),2.0)+pow(NODE_Y(ptrN),2.0));
				/*Find the maximum and minimum radial coordinates of the boundary*/
				Rmax = MAX(Rmax,RN);
				Rmin = MIN(Rmin,RN);
			}
	 }end_f_loop(f, thread)

	 /*Store the max and min radial coordinate of the boundary using parallelisation macros*/
	 Rmax=PRF_GRHIGH1(Rmax);
	 Rmin=PRF_GRLOW1(Rmin);

	 
	 /*Total length of the boundary*/
	 L=Rmax-Rmin;
	 
	 begin_f_loop(f,thread)
	 {
		F_CENTROID(x, f, thread);
		/*Radial coordinates of the centroids*/
		RR=sqrt(pow(x[0],2.0)+pow(x[1],2.0));
		Etha=(RR-Rmin)/L;
		/*Newtonian profile calculation*/
		F_PROFILE(f, thread, position) = Rmin*OMEGA*(1.0-(3.0*Etha))*(1.0-Etha);
	 }end_f_loop(f, thread)
	 
	 #endif
 }