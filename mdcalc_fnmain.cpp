
#include "mdcalc_header.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>


using namespace std;

//-------------------------------------------------------------------------------------------------------------------------------
//function to copy vectors 
void MD::copy_vector(vector <string> test){

	for(int i=0;i<test.size();i++){
		tokens[i] = test[i];
	}
}


//-------------------------------------------------------------------------------------------------------------------------------
void MD::make_vector(){
	for(int i=0;i<5;i++)
		tokens.push_back("dummy");
}

//-------------------------------------------------------------------------------------------------------------------------------
void MD::calc_lat(){
	d_alat = atof(tokens[2].c_str());	
}

//-------------------------------------------------------------------------------------------------------------------------------
void MD::make_potential(){

	
	ofstream fout("potential.file");
	fout<<"pair_style eam"<<endl;
	fout<<"pair_coeff 1 1 ../"<<tokens[3]<<endl;
	fout.close();
}

//-------------------------------------------------------------------------------------------------------------------------------



    
   
    //checkin if gsfe 2 is selected------------------------------------------------------------------------------------






    //checkin if gsfe 3 is selected-------------------------------------------------------------------------------------



    

//-------------------------------------------------------------------------------------------------------------------------------
//make crystal structure
void MD::make_xtal(){

	if ((tokens[0]=="cohesiveeng")||(tokens[0]=="bulkmod"))
    {
        if(tokens[1]=="bcc")
        {   xtalf_rname = "gen.bcc.atoms";
            xtalf_wname = "generic.bcc.atoms";
            mockxtal = "BCC";

            enx = 1; eny = 1; enz = 1;
            no_atoms = 2*enx*eny*enz;
            cord_no = 2;

            xfact = 1;
            yfact = 1;
            zfact = 1;

            unitcellX[1]=0;
            unitcellY[1]=0;
            unitcellZ[1]=0;

            unitcellX[2]=0.5*d_alat;
            unitcellY[2]=0.5*d_alat;
            unitcellZ[2]=0.5*d_alat;

            pre_make_xtal(0);
        }

        if(tokens[1]=="fcc")
        {   xtalf_rname = "gen.fcc.atoms";
            xtalf_wname = "generic.fcc.atoms";
            mockxtal = "FCC";

            enx = 1; eny = 1; enz = 1;
            no_atoms = 4*enx*eny*enz;
            cord_no = 4;

            xfact = 1;
            yfact = 1;
            zfact = 1;

            unitcellX[1]=0;
            unitcellY[1]=0;
            unitcellZ[1]=0;

            unitcellX[2]=0.5*d_alat;
            unitcellY[2]=0;
            unitcellZ[2]=0.5*d_alat;

            unitcellX[3]=0;
            unitcellY[3]=0.5*d_alat;
            unitcellZ[3]=0.5*d_alat;

            unitcellX[4]=0.5*d_alat;
            unitcellY[4]=0.5*d_alat;
            unitcellZ[4]=0;

            pre_make_xtal(0);
        }


        if(tokens[1]=="hcp")
        {   xtalf_rname = "gen.hcp.atoms";
            xtalf_wname = "generic.hcp.atoms";
            mockxtal = "HCP";

            enx = 1; eny = 1; enz = 1;
            no_atoms = 4*enx*eny*enz;
            cord_no = 4;

            xfact = 1;
            yfact = sqrt(3);
            zfact = 1.633;

            unitcellX[1]=0;
            unitcellY[1]=0;
            unitcellZ[1]=0;

            unitcellX[2]=0.5*d_alat;
            unitcellY[2]=0.5*d_alat*sqrt(3);
            unitcellZ[2]=0;

            unitcellX[3]=0.5*d_alat;
            unitcellY[3]=sqrt(3)*d_alat/6;
            unitcellZ[3]=0.5*zfact*d_alat;

            unitcellX[4]=0;
            unitcellY[4]=2*d_alat/sqrt(3);
            unitcellZ[4]=0.5*zfact*d_alat;

            pre_make_xtal(0);
        }
    }

// elastic constants----*************************************************************

    if(tokens[0]=="elasticconst")
    {
        if(tokens[1]=="bcc")
        {   xtalf_rname = "gen.bcc.atoms";
            xtalf_wname = "generic2.bcc.atoms";
            mockxtal = "BCC";

            enx = 1; eny = 1; enz = 1;
            no_atoms = 2*enx*eny*enz;
            cord_no = 2;

            xfact = 1;
            yfact = 1;
            zfact = 1;

            unitcellX[1]=0;
            unitcellY[1]=0;
            unitcellZ[1]=0;

            unitcellX[2]=0.5*d_alat;
            unitcellY[2]=0.5*d_alat;
            unitcellZ[2]=0.5*d_alat;

            pre_make_xtal(1);

        }

        if(tokens[1]=="fcc")
        {   xtalf_rname = "gen.fcc.atoms";
            xtalf_wname = "generic2.fcc.atoms";
            mockxtal = "FCC";

            enx = 1; eny = 1; enz = 1;
            no_atoms = 4*enx*eny*enz;
            cord_no = 4;

            xfact = 1;
            yfact = 1;
            zfact = 1;

            unitcellX[1]=0;
            unitcellY[1]=0;
            unitcellZ[1]=0;

            unitcellX[2]=0.5*d_alat;
            unitcellY[2]=0;
            unitcellZ[2]=0.5*d_alat;

            unitcellX[3]=0;
            unitcellY[3]=0.5*d_alat;
            unitcellZ[3]=0.5*d_alat;

            unitcellX[4]=0.5*d_alat;
            unitcellY[4]=0.5*d_alat;
            unitcellZ[4]=0;

            pre_make_xtal(1);

        }

        if(tokens[1]=="hcp")
        {   xtalf_rname = "gen.hcp.atoms";
            xtalf_wname = "generic2.hcp.atoms";
            mockxtal = "HCP";

            enx = 1; eny = 1; enz = 1;
            no_atoms = 4*enx*eny*enz;
            cord_no = 4;

            xfact = 1;
            yfact = sqrt(3);
            zfact = 1.633;

            unitcellX[1]=0;
            unitcellY[1]=0;
            unitcellZ[1]=0;

            unitcellX[2]=0.5*d_alat;
            unitcellY[2]=0.5*d_alat*sqrt(3);
            unitcellZ[2]=0;

            unitcellX[3]=0.5*d_alat;
            unitcellY[3]=sqrt(3)*d_alat/6;
            unitcellZ[3]=0.5*zfact*d_alat;

            unitcellX[4]=0;
            unitcellY[4]=2*d_alat/sqrt(3);
            unitcellZ[4]=0.5*zfact*d_alat;


            pre_make_xtal(1);

        }
    }

    if(tokens[0]=="gsfe1")
    {
        if(tokens[1]=="bcc")
        {   xtalf_rname = "gen.bccgsfe1.atoms";
            xtalf_wname = "generic.bccgsfe1.atoms";
            mockxtal = "BCC";

            enx = 2; eny = 10; enz = 2;
            no_atoms = 12*enx*eny*enz;
            cord_no = 12;

            xfact = 1.732;
            yfact = 1.414;
            zfact = 2.4494801;

            unitcellX[1]=1.4433751*d_alat;
            unitcellY[1]=0*d_alat;
            unitcellZ[1]=1.632995*d_alat;

            unitcellX[2]=0.288675*d_alat;
            unitcellY[2]=0*d_alat;
            unitcellZ[2]=0.816495*d_alat;

            unitcellX[3]=0.866025*d_alat;
            unitcellY[3]=0*d_alat;
            unitcellZ[3]=5.0067902e-06*d_alat;

            unitcellX[4]=0.288675*d_alat;
            unitcellY[4]=0.7071001*d_alat;
            unitcellZ[4]=2.041235*d_alat;

            unitcellX[5]=0.57735*d_alat;
            unitcellY[5]=0*d_alat;
            unitcellZ[5]=1.63299*d_alat;

            unitcellX[6]=0.866025*d_alat;
            unitcellY[6]=0.7071001*d_alat;
            unitcellZ[6]=1.224745*d_alat;

            unitcellX[7]=1.1547*d_alat;
            unitcellY[7]=0*d_alat;
            unitcellZ[7]=0.8165*d_alat;

            unitcellX[8]=1.4433751*d_alat;
            unitcellY[8]=0.7071001*d_alat;
            unitcellZ[8]=0.408255*d_alat;

            unitcellX[9]=1.1547*d_alat;
            unitcellY[9]=0.7071*d_alat;
            unitcellZ[9]=2.04124*d_alat;

            unitcellX[10]=0*d_alat;
            unitcellY[10]=0*d_alat;
            unitcellZ[10]=0*d_alat;

            unitcellX[11]=0*d_alat;
            unitcellY[11]=0.7071*d_alat;
            unitcellZ[11]=1.22474*d_alat;

            unitcellX[12]=0.57735*d_alat;
            unitcellY[12]=0.7071*d_alat;
            unitcellZ[12]=0.40825*d_alat;


            pre_make_xtal(0);
        }

        if(tokens[1]=="fcc")
        {   xtalf_rname = "gen.fccgsfe1.atoms";
            xtalf_wname = "generic.fccgsfe1.atoms";
            mockxtal = "FCC";

            enx = 2; eny = 8; enz = 2;
            no_atoms = 24*enx*eny*enz;
            cord_no = 24;

            xfact = 1.732;
            yfact = 1.414;
            zfact = 2.4494801;

            unitcellX[1]=0*d_alat;
            unitcellY[1]=0*d_alat;
            unitcellZ[1]=0*d_alat;

            unitcellX[2]=0*d_alat;
            unitcellY[2]=0.7070999742*d_alat;
            unitcellZ[2]=0*d_alat;

            unitcellX[3]=1.154700041*d_alat;
            unitcellY[3]=0*d_alat;
            unitcellZ[3]=2.041239977*d_alat;

            unitcellX[4]=0*d_alat;
            unitcellY[4]=0*d_alat;
            unitcellZ[4]=1.224740028*d_alat;

            unitcellX[5]=0.5773500204*d_alat;
            unitcellY[5]=0*d_alat;
            unitcellZ[5]=0.4082500041*d_alat;

            unitcellX[6]=0.57735008*d_alat;
            unitcellY[6]=0.3535500169*d_alat;
            unitcellZ[6]=2.245359898*d_alat;

            unitcellX[7]=0*d_alat;
            unitcellY[7]=0.3535499871*d_alat;
            unitcellZ[7]=1.837110043*d_alat;

            unitcellX[8]=0.5773500204*d_alat;
            unitcellY[8]=0*d_alat;
            unitcellZ[8]=1.632990003*d_alat;

            unitcellX[9]=1.154700041*d_alat;
            unitcellY[9]=0.3535500169*d_alat;
            unitcellZ[9]=1.428869963*d_alat;

            unitcellX[10]=0.5773500204*d_alat;
            unitcellY[10]=0.7071000338*d_alat;
            unitcellZ[10]=1.632990003*d_alat;

            unitcellX[11]=0.5773500204*d_alat;
            unitcellY[11]=0.3535499871*d_alat;
            unitcellZ[11]=1.020619988*d_alat;

            unitcellX[12]=1.154700041*d_alat;
            unitcellY[12]=0*d_alat;
            unitcellZ[12]=0.8165000081*d_alat;

            unitcellX[13]=1.154700041*d_alat;
            unitcellY[13]=0.7071000338*d_alat;
            unitcellZ[13]=0.8165000081*d_alat;

            unitcellX[14]=1.154700041*d_alat;
            unitcellY[14]=0.3535499871*d_alat;
            unitcellZ[14]=0.2041300237*d_alat;

            unitcellX[15]=0.57735008*d_alat;
            unitcellY[15]=1.060649991*d_alat;
            unitcellZ[15]=2.245359898*d_alat;

            unitcellX[16]=1.154700041*d_alat;
            unitcellY[16]=0.7070999742*d_alat;
            unitcellZ[16]=2.041239977*d_alat;

            unitcellX[17]=1.154700041*d_alat;
            unitcellY[17]=1.060649991*d_alat;
            unitcellZ[17]=1.428869963*d_alat;

            unitcellX[18]=0*d_alat;
            unitcellY[18]=0.3535499871*d_alat;
            unitcellZ[18]=0.6123700142*d_alat;

            unitcellX[19]=0*d_alat;
            unitcellY[19]=1.060649991*d_alat;
            unitcellZ[19]=1.837110043*d_alat;

            unitcellX[20]=0*d_alat;
            unitcellY[20]=0.7070999742*d_alat;
            unitcellZ[20]=1.224740028*d_alat;

            unitcellX[21]=0.5773500204*d_alat;
            unitcellY[21]=1.060649991*d_alat;
            unitcellZ[21]=1.020619988*d_alat;

            unitcellX[22]=0*d_alat;
            unitcellY[22]=1.060649991*d_alat;
            unitcellZ[22]=0.6123700142*d_alat;

            unitcellX[23]=0.5773500204*d_alat;
            unitcellY[23]=0.7070999742*d_alat;
            unitcellZ[23]=0.4082500041*d_alat;

            unitcellX[24]=1.154700041*d_alat;
            unitcellY[24]=1.060649991*d_alat;
            unitcellZ[24]=0.2041300237*d_alat;

            pre_make_xtal(0);
        }

        if(tokens[1]=="hcp")
        {   xtalf_rname = "gen.hcpgsfe1.atoms";
            xtalf_wname = "generic.hcpgsfe1.atoms";
            mockxtal = "HCP";

            xmin = 0; xmax = 2;
            ymin = -6.532000065; ymax = 9.798000336;
            zmin = 0; zmax = 3.463999987;
            no_atoms = 160;

            pre_make_xtal(0);
        }
    }

    if(tokens[0]=="gsfe2")
    {
        if(tokens[1]=="bcc")
        {   xtalf_rname = "gen.bccgsfe2.atoms";
            xtalf_wname = "generic.bccgsfe2.atoms";
            mockxtal = "BCC";

            enx = 2; eny = 10; enz = 2;
            no_atoms = 12*enx*eny*enz;
            cord_no = 12;

            xfact = 1.732;
            yfact = 2.4494801;
            zfact = 1.414;

            unitcellX[1]=0*d_alat;
            unitcellY[1]=0*d_alat;
            unitcellZ[1]=0*d_alat;

            unitcellX[2]=1.443375111*d_alat;
            unitcellY[2]=1.632995009*d_alat;
            unitcellZ[2]=0*d_alat;

            unitcellX[3]=0.2886750102*d_alat;
            unitcellY[3]=0.8164950013*d_alat;
            unitcellZ[3]=0*d_alat;

            unitcellX[4]=0.8660250306*d_alat;
            unitcellY[4]=5.006790161e-06*d_alat;
            unitcellZ[4]=0*d_alat;

            unitcellX[5]=0.2886750102*d_alat;
            unitcellY[5]=2.04123497*d_alat;
            unitcellZ[5]=0.7071000338*d_alat;

            unitcellX[6]=0.5773500204*d_alat;
            unitcellY[6]=1.632990003*d_alat;
            unitcellZ[6]=0*d_alat;

            unitcellX[7]=0.8660250306*d_alat;
            unitcellY[7]=1.224745035*d_alat;
            unitcellZ[7]=0.7071000338*d_alat;

            unitcellX[8]=1.154700041*d_alat;
            unitcellY[8]=0.8165000081*d_alat;
            unitcellZ[8]=0*d_alat;

            unitcellX[9]=1.443375111*d_alat;
            unitcellY[9]=0.4082550108*d_alat;
            unitcellZ[9]=0.7071000338*d_alat;

            unitcellX[10]=1.154700041*d_alat;
            unitcellY[10]=2.041239977*d_alat;
            unitcellZ[10]=0.7070999742*d_alat;

            unitcellX[11]=0*d_alat;
            unitcellY[11]=1.224740028*d_alat;
            unitcellZ[11]=0.7070999742*d_alat;

            unitcellX[12]=0.5773500204*d_alat;
            unitcellY[12]=0.4082500041*d_alat;
            unitcellZ[12]=0.7070999742*d_alat;

            pre_make_xtal(0);
        }

        if(tokens[1]=="fcc")
        {   xtalf_rname = "gen.fccgsfe2.atoms";
            xtalf_wname = "generic.fccgsfe2.atoms";
            mockxtal = "FCC";

            enx = 2; eny = 8; enz = 2;
            no_atoms = 24*enx*eny*enz;
            cord_no = 24;

            xfact = 1.732;
            yfact = 2.4494801;
            zfact = 1.414;

            unitcellX[1]=0*d_alat;
            unitcellY[1]=0*d_alat;
            unitcellZ[1]=0*d_alat;

            unitcellX[2]=0*d_alat;
            unitcellY[2]=0*d_alat;
            unitcellZ[2]=0.7070999742*d_alat;

            unitcellX[3]=1.154700041*d_alat;
            unitcellY[3]=2.041239977*d_alat;
            unitcellZ[3]=0*d_alat;

            unitcellX[4]=0*d_alat;
            unitcellY[4]=1.224740028*d_alat;
            unitcellZ[4]=0*d_alat;

            unitcellX[5]=0.5773500204*d_alat;
            unitcellY[5]=0.4082500041*d_alat;
            unitcellZ[5]=0*d_alat;

            unitcellX[6]=0.57735008*d_alat;
            unitcellY[6]=2.245359898*d_alat;
            unitcellZ[6]=0.3535500169*d_alat;

            unitcellX[7]=0*d_alat;
            unitcellY[7]=1.837110043*d_alat;
            unitcellZ[7]=0.3535499871*d_alat;

            unitcellX[8]=0.5773500204*d_alat;
            unitcellY[8]=1.632990003*d_alat;
            unitcellZ[8]=0*d_alat;

            unitcellX[9]=1.154700041*d_alat;
            unitcellY[9]=1.428869963*d_alat;
            unitcellZ[9]=0.3535500169*d_alat;

            unitcellX[10]=0.5773500204*d_alat;
            unitcellY[10]=1.632990003*d_alat;
            unitcellZ[10]=0.7071000338*d_alat;

            unitcellX[11]=0.5773500204*d_alat;
            unitcellY[11]=1.020619988*d_alat;
            unitcellZ[11]=0.3535499871*d_alat;

            unitcellX[12]=1.154700041*d_alat;
            unitcellY[12]=0.8165000081*d_alat;
            unitcellZ[12]=0*d_alat;

            unitcellX[13]=1.154700041*d_alat;
            unitcellY[13]=0.8165000081*d_alat;
            unitcellZ[13]=0.7071000338*d_alat;

            unitcellX[14]=1.154700041*d_alat;
            unitcellY[14]=0.2041300237*d_alat;
            unitcellZ[14]=0.3535499871*d_alat;

            unitcellX[15]=0.57735008*d_alat;
            unitcellY[15]=2.245359898*d_alat;
            unitcellZ[15]=1.060649991*d_alat;

            unitcellX[16]=1.154700041*d_alat;
            unitcellY[16]=2.041239977*d_alat;
            unitcellZ[16]=0.7070999742*d_alat;

            unitcellX[17]=1.154700041*d_alat;
            unitcellY[17]=1.428869963*d_alat;
            unitcellZ[17]=1.060649991*d_alat;

            unitcellX[18]=0*d_alat;
            unitcellY[18]=0.6123700142*d_alat;
            unitcellZ[18]=0.3535499871*d_alat;

            unitcellX[19]=0*d_alat;
            unitcellY[19]=1.837110043*d_alat;
            unitcellZ[19]=1.060649991*d_alat;

            unitcellX[20]=0*d_alat;
            unitcellY[20]=1.224740028*d_alat;
            unitcellZ[20]=0.7070999742*d_alat;

            unitcellX[21]=0.5773500204*d_alat;
            unitcellY[21]=1.020619988*d_alat;
            unitcellZ[21]=1.060649991*d_alat;

            unitcellX[22]=0*d_alat;
            unitcellY[22]=0.6123700142*d_alat;
            unitcellZ[22]=1.060649991*d_alat;

            unitcellX[23]=0.5773500204*d_alat;
            unitcellY[23]=0.4082500041*d_alat;
            unitcellZ[23]=0.7070999742*d_alat;

            unitcellX[24]=1.154700041*d_alat;
            unitcellY[24]=0.2041300237*d_alat;
            unitcellZ[24]=1.060649991*d_alat;

            pre_make_xtal(0);
        }

        if(tokens[1]=="hcp")
        {   xtalf_rname = "gen.hcpgsfe2.atoms";
            xtalf_wname = "generic.hcpgsfe2.atoms";
            mockxtal = "HCP";

            xmin = -4; xmax = 6;
            ymin = 0; ymax = 3.266000032;
            zmin = 0; zmax = 3.463999987;
            no_atoms = 160;

            pre_make_xtal(0);
        }
    }

    if(tokens[0]=="gsfe3")
    {
        if(tokens[1]=="bcc")
        {   xtalf_rname = "gen.bccgsfe3.atoms";
            xtalf_wname = "generic.bccgsfe3.atoms";
            mockxtal = "BCC";

            enx = 2; eny = 6; enz = 2;
            no_atoms = 84*enx*eny*enz;
            cord_no = 84;

            xfact = 1.731999993;
            yfact = 3.741650105;
            zfact = 6.48074007;

            unitcellX[1]=0*d_alat;
            unitcellY[1]=0*d_alat;
            unitcellZ[1]=0*d_alat;
            unitcellX[2]=1.44338*d_alat;
            unitcellY[2]=2.93986*d_alat;
            unitcellZ[2]=6.3263*d_alat;
            unitcellX[3]=0.288675*d_alat;
            unitcellY[3]=2.13808*d_alat;
            unitcellZ[3]=6.172*d_alat;
            unitcellX[4]=0.866025*d_alat;
            unitcellY[4]=1.3363*d_alat;
            unitcellZ[4]=6.0177*d_alat;
            unitcellX[5]=1.44337*d_alat;
            unitcellY[5]=0.53452*d_alat;
            unitcellZ[5]=5.8634*d_alat;
            unitcellX[6]=0.288675*d_alat;
            unitcellY[6]=3.47438*d_alat;
            unitcellZ[6]=5.7091*d_alat;
            unitcellX[7]=0.57735*d_alat;
            unitcellY[7]=2.93986*d_alat;
            unitcellZ[7]=6.3263*d_alat;
            unitcellX[8]=0.866025*d_alat;
            unitcellY[8]=2.6726*d_alat;
            unitcellZ[8]=5.5548*d_alat;
            unitcellX[9]=1.1547*d_alat;
            unitcellY[9]=2.13808*d_alat;
            unitcellZ[9]=6.172*d_alat;
            unitcellX[10]=1.44338*d_alat;
            unitcellY[10]=1.87082*d_alat;
            unitcellZ[10]=5.4005*d_alat;
            unitcellX[11]=1.1547*d_alat;
            unitcellY[11]=3.47438*d_alat;
            unitcellZ[11]=5.7091*d_alat;
            unitcellX[12]=1.44337*d_alat;
            unitcellY[12]=3.20712*d_alat;
            unitcellZ[12]=4.9376*d_alat;
            unitcellX[13]=1.19209e-07*d_alat;
            unitcellY[13]=1.3363*d_alat;
            unitcellZ[13]=6.0177*d_alat;
            unitcellX[14]=0.288675*d_alat;
            unitcellY[14]=1.06904*d_alat;
            unitcellZ[14]=5.2462*d_alat;
            unitcellX[15]=0.57735*d_alat;
            unitcellY[15]=0.53452*d_alat;
            unitcellZ[15]=5.8634*d_alat;
            unitcellX[16]=0.866025*d_alat;
            unitcellY[16]=0.26726*d_alat;
            unitcellZ[16]=5.0919*d_alat;
            unitcellX[17]=1.19209e-07*d_alat;
            unitcellY[17]=2.6726*d_alat;
            unitcellZ[17]=5.5548*d_alat;
            unitcellX[18]=0.288675*d_alat;
            unitcellY[18]=2.40534*d_alat;
            unitcellZ[18]=4.7833*d_alat;
            unitcellX[19]=0.57735*d_alat;
            unitcellY[19]=1.87082*d_alat;
            unitcellZ[19]=5.4005*d_alat;
            unitcellX[20]=0.866025*d_alat;
            unitcellY[20]=1.60356*d_alat;
            unitcellZ[20]=4.629*d_alat;
            unitcellX[21]=1.1547*d_alat;
            unitcellY[21]=1.06904*d_alat;
            unitcellZ[21]=5.2462*d_alat;
            unitcellX[22]=1.44338*d_alat;
            unitcellY[22]=0.80178*d_alat;
            unitcellZ[22]=4.4747*d_alat;
            unitcellX[23]=0.57735*d_alat;
            unitcellY[23]=3.20712*d_alat;
            unitcellZ[23]=4.9376*d_alat;
            unitcellX[24]=0.866025*d_alat;
            unitcellY[24]=2.93986*d_alat;
            unitcellZ[24]=4.1661*d_alat;
            unitcellX[25]=1.1547*d_alat;
            unitcellY[25]=2.40534*d_alat;
            unitcellZ[25]=4.7833*d_alat;
            unitcellX[26]=1.44338*d_alat;
            unitcellY[26]=2.13808*d_alat;
            unitcellZ[26]=4.0118*d_alat;
            unitcellX[27]=1.44338*d_alat;
            unitcellY[27]=3.47438*d_alat;
            unitcellZ[27]=3.5489*d_alat;
            unitcellX[28]=0*d_alat;
            unitcellY[28]=0.26726*d_alat;
            unitcellZ[28]=5.0919*d_alat;
            unitcellX[29]=0.288675*d_alat;
            unitcellY[29]=5.96046e-08*d_alat;
            unitcellZ[29]=4.3204*d_alat;
            unitcellX[30]=0*d_alat;
            unitcellY[30]=1.60356*d_alat;
            unitcellZ[30]=4.629*d_alat;
            unitcellX[31]=0.288675*d_alat;
            unitcellY[31]=1.3363*d_alat;
            unitcellZ[31]=3.8575*d_alat;
            unitcellX[32]=0.57735*d_alat;
            unitcellY[32]=0.80178*d_alat;
            unitcellZ[32]=4.4747*d_alat;
            unitcellX[33]=0.866025*d_alat;
            unitcellY[33]=0.53452*d_alat;
            unitcellZ[33]=3.7032*d_alat;
            unitcellX[34]=1.1547*d_alat;
            unitcellY[34]=0*d_alat;
            unitcellZ[34]=4.3204*d_alat;
            unitcellX[35]=0*d_alat;
            unitcellY[35]=2.93986*d_alat;
            unitcellZ[35]=4.1661*d_alat;
            unitcellX[36]=0.288675*d_alat;
            unitcellY[36]=2.6726*d_alat;
            unitcellZ[36]=3.3946*d_alat;
            unitcellX[37]=0.57735*d_alat;
            unitcellY[37]=2.13808*d_alat;
            unitcellZ[37]=4.0118*d_alat;
            unitcellX[38]=0.866025*d_alat;
            unitcellY[38]=1.87082*d_alat;
            unitcellZ[38]=3.2403*d_alat;
            unitcellX[39]=1.1547*d_alat;
            unitcellY[39]=1.3363*d_alat;
            unitcellZ[39]=3.8575*d_alat;
            unitcellX[40]=1.44337*d_alat;
            unitcellY[40]=1.06904*d_alat;
            unitcellZ[40]=3.086*d_alat;
            unitcellX[41]=0.57735*d_alat;
            unitcellY[41]=3.47438*d_alat;
            unitcellZ[41]=3.5489*d_alat;
            unitcellX[42]=0.866025*d_alat;
            unitcellY[42]=3.20712*d_alat;
            unitcellZ[42]=2.7774*d_alat;
            unitcellX[43]=1.1547*d_alat;
            unitcellY[43]=2.6726*d_alat;
            unitcellZ[43]=3.3946*d_alat;
            unitcellX[44]=1.44338*d_alat;
            unitcellY[44]=2.40534*d_alat;
            unitcellZ[44]=2.6231*d_alat;
            unitcellX[45]=0*d_alat;
            unitcellY[45]=0.53452*d_alat;
            unitcellZ[45]=3.7032*d_alat;
            unitcellX[46]=0.288675*d_alat;
            unitcellY[46]=0.26726*d_alat;
            unitcellZ[46]=2.9317*d_alat;
            unitcellX[47]=0*d_alat;
            unitcellY[47]=1.87082*d_alat;
            unitcellZ[47]=3.2403*d_alat;
            unitcellX[48]=0.288675*d_alat;
            unitcellY[48]=1.60356*d_alat;
            unitcellZ[48]=2.4688*d_alat;
            unitcellX[49]=0.57735*d_alat;
            unitcellY[49]=1.06904*d_alat;
            unitcellZ[49]=3.086*d_alat;
            unitcellX[50]=0.866025*d_alat;
            unitcellY[50]=0.80178*d_alat;
            unitcellZ[50]=2.3145*d_alat;
            unitcellX[51]=1.1547*d_alat;
            unitcellY[51]=0.26726*d_alat;
            unitcellZ[51]=2.9317*d_alat;
            unitcellX[52]=1.44338*d_alat;
            unitcellY[52]=-2.98023e-08*d_alat;
            unitcellZ[52]=2.1602*d_alat;
            unitcellX[53]=0*d_alat;
            unitcellY[53]=3.20712*d_alat;
            unitcellZ[53]=2.7774*d_alat;
            unitcellX[54]=0.288675*d_alat;
            unitcellY[54]=2.93986*d_alat;
            unitcellZ[54]=2.0059*d_alat;
            unitcellX[55]=0.57735*d_alat;
            unitcellY[55]=2.40534*d_alat;
            unitcellZ[55]=2.6231*d_alat;
            unitcellX[56]=0.866025*d_alat;
            unitcellY[56]=2.13808*d_alat;
            unitcellZ[56]=1.8516*d_alat;
            unitcellX[57]=1.1547*d_alat;
            unitcellY[57]=1.60356*d_alat;
            unitcellZ[57]=2.4688*d_alat;
            unitcellX[58]=1.44338*d_alat;
            unitcellY[58]=1.3363*d_alat;
            unitcellZ[58]=1.6973*d_alat;
            unitcellX[59]=0.866025*d_alat;
            unitcellY[59]=3.47438*d_alat;
            unitcellZ[59]=1.3887*d_alat;
            unitcellX[60]=1.1547*d_alat;
            unitcellY[60]=2.93986*d_alat;
            unitcellZ[60]=2.0059*d_alat;
            unitcellX[61]=1.44338*d_alat;
            unitcellY[61]=2.6726*d_alat;
            unitcellZ[61]=1.2344*d_alat;
            unitcellX[62]=0*d_alat;
            unitcellY[62]=0.80178*d_alat;
            unitcellZ[62]=2.3145*d_alat;
            unitcellX[63]=0.288675*d_alat;
            unitcellY[63]=0.53452*d_alat;
            unitcellZ[63]=1.543*d_alat;
            unitcellX[64]=0.57735*d_alat;
            unitcellY[64]=0*d_alat;
            unitcellZ[64]=2.1602*d_alat;
            unitcellX[65]=0*d_alat;
            unitcellY[65]=2.13808*d_alat;
            unitcellZ[65]=1.8516*d_alat;
            unitcellX[66]=0.288675*d_alat;
            unitcellY[66]=1.87082*d_alat;
            unitcellZ[66]=1.0801*d_alat;
            unitcellX[67]=0.57735*d_alat;
            unitcellY[67]=1.3363*d_alat;
            unitcellZ[67]=1.6973*d_alat;
            unitcellX[68]=0.866025*d_alat;
            unitcellY[68]=1.06904*d_alat;
            unitcellZ[68]=0.9258*d_alat;
            unitcellX[69]=1.1547*d_alat;
            unitcellY[69]=0.53452*d_alat;
            unitcellZ[69]=1.543*d_alat;
            unitcellX[70]=1.44338*d_alat;
            unitcellY[70]=0.26726*d_alat;
            unitcellZ[70]=0.7715*d_alat;
            unitcellX[71]=0*d_alat;
            unitcellY[71]=3.47438*d_alat;
            unitcellZ[71]=1.3887*d_alat;
            unitcellX[72]=0.288675*d_alat;
            unitcellY[72]=3.20712*d_alat;
            unitcellZ[72]=0.6172*d_alat;
            unitcellX[73]=0.57735*d_alat;
            unitcellY[73]=2.6726*d_alat;
            unitcellZ[73]=1.2344*d_alat;
            unitcellX[74]=0.866025*d_alat;
            unitcellY[74]=2.40534*d_alat;
            unitcellZ[74]=0.4629*d_alat;
            unitcellX[75]=1.1547*d_alat;
            unitcellY[75]=1.87082*d_alat;
            unitcellZ[75]=1.0801*d_alat;
            unitcellX[76]=1.44338*d_alat;
            unitcellY[76]=1.60356*d_alat;
            unitcellZ[76]=0.3086*d_alat;
            unitcellX[77]=1.1547*d_alat;
            unitcellY[77]=3.20712*d_alat;
            unitcellZ[77]=0.6172*d_alat;
            unitcellX[78]=0*d_alat;
            unitcellY[78]=1.06904*d_alat;
            unitcellZ[78]=0.9258*d_alat;
            unitcellX[79]=0.288675*d_alat;
            unitcellY[79]=0.80178*d_alat;
            unitcellZ[79]=0.1543*d_alat;
            unitcellX[80]=0.57735*d_alat;
            unitcellY[80]=0.26726*d_alat;
            unitcellZ[80]=0.7715*d_alat;
            unitcellX[81]=0.866025*d_alat;
            unitcellY[81]=0*d_alat;
            unitcellZ[81]=-1.49012e-08*d_alat;
            unitcellX[82]=0*d_alat;
            unitcellY[82]=2.40534*d_alat;
            unitcellZ[82]=0.4629*d_alat;
            unitcellX[83]=0.57735*d_alat;
            unitcellY[83]=1.60356*d_alat;
            unitcellZ[83]=0.3086*d_alat;
            unitcellX[84]=1.1547*d_alat;
            unitcellY[84]=0.80178*d_alat;
            unitcellZ[84]=0.1543*d_alat;

            pre_make_xtal(0);
        }

        if(tokens[1]=="fcc")
        {   xtalf_rname = "gen.fccgsfe3.atoms";
            xtalf_wname = "generic.fccgsfe3.atoms";
            mockxtal = "FCC";

            enx = 2; eny = 4; enz = 2;
            no_atoms = 168*enx*eny*enz;
            cord_no = 168;

            xfact = 1.731999993;
            yfact = 3.741650105;
            zfact = 6.48074007;

            unitcellX[1]=0*d_alat;
            unitcellY[1]=0*d_alat;
            unitcellZ[1]=0*d_alat;
            unitcellX[2]=-1.19209e-07*d_alat;
            unitcellY[2]=2.53897*d_alat;
            unitcellZ[2]=6.24951*d_alat;
            unitcellX[3]=0.57735*d_alat;
            unitcellY[3]=1.73719*d_alat;
            unitcellZ[3]=6.0952*d_alat;
            unitcellX[4]=0.57735*d_alat;
            unitcellY[4]=1.06904*d_alat;
            unitcellZ[4]=6.32666*d_alat;
            unitcellX[5]=1.1547*d_alat;
            unitcellY[5]=0.93541*d_alat;
            unitcellZ[5]=5.94091*d_alat;
            unitcellX[6]=1.1547*d_alat;
            unitcellY[6]=0.26726*d_alat;
            unitcellZ[6]=6.17236*d_alat;
            unitcellX[7]=0.57735*d_alat;
            unitcellY[7]=3.60801*d_alat;
            unitcellZ[7]=6.09523*d_alat;
            unitcellX[8]=0*d_alat;
            unitcellY[8]=3.20712*d_alat;
            unitcellZ[8]=6.01805*d_alat;
            unitcellX[9]=0.57735*d_alat;
            unitcellY[9]=2.93986*d_alat;
            unitcellZ[9]=6.32668*d_alat;
            unitcellX[10]=1.1547*d_alat;
            unitcellY[10]=2.80623*d_alat;
            unitcellZ[10]=5.94093*d_alat;
            unitcellX[11]=0.57735*d_alat;
            unitcellY[11]=3.07349*d_alat;
            unitcellZ[11]=5.6323*d_alat;
            unitcellX[12]=0.57735*d_alat;
            unitcellY[12]=2.40534*d_alat;
            unitcellZ[12]=5.86375*d_alat;
            unitcellX[13]=1.1547*d_alat;
            unitcellY[13]=2.13808*d_alat;
            unitcellZ[13]=6.17238*d_alat;
            unitcellX[14]=1.1547*d_alat;
            unitcellY[14]=2.27171*d_alat;
            unitcellZ[14]=5.47799*d_alat;
            unitcellX[15]=1.1547*d_alat;
            unitcellY[15]=1.60356*d_alat;
            unitcellZ[15]=5.70945*d_alat;
            unitcellX[16]=1.1547*d_alat;
            unitcellY[16]=3.47438*d_alat;
            unitcellZ[16]=5.70947*d_alat;
            unitcellX[17]=1.1547*d_alat;
            unitcellY[17]=3.60801*d_alat;
            unitcellZ[17]=5.01509*d_alat;
            unitcellX[18]=1.1547*d_alat;
            unitcellY[18]=2.93986*d_alat;
            unitcellZ[18]=5.24654*d_alat;
            unitcellX[19]=1.78814e-07*d_alat;
            unitcellY[19]=0.66815*d_alat;
            unitcellZ[19]=6.24949*d_alat;
            unitcellX[20]=1.19209e-07*d_alat;
            unitcellY[20]=0.13363*d_alat;
            unitcellZ[20]=5.78656*d_alat;
            unitcellX[21]=0*d_alat;
            unitcellY[21]=2.00445*d_alat;
            unitcellZ[21]=5.78657*d_alat;
            unitcellX[22]=1.19209e-07*d_alat;
            unitcellY[22]=1.3363*d_alat;
            unitcellZ[22]=6.01803*d_alat;
            unitcellX[23]=0.57735*d_alat;
            unitcellY[23]=1.20267*d_alat;
            unitcellZ[23]=5.63228*d_alat;
            unitcellX[24]=0*d_alat;
            unitcellY[24]=1.46993*d_alat;
            unitcellZ[24]=5.32365*d_alat;
            unitcellX[25]=5.96046e-08*d_alat;
            unitcellY[25]=0.80178*d_alat;
            unitcellZ[25]=5.5551*d_alat;
            unitcellX[26]=0.57735*d_alat;
            unitcellY[26]=0.53452*d_alat;
            unitcellZ[26]=5.86373*d_alat;
            unitcellX[27]=1.1547*d_alat;
            unitcellY[27]=0.40089*d_alat;
            unitcellZ[27]=5.47797*d_alat;
            unitcellX[28]=0.57735*d_alat;
            unitcellY[28]=0.66815*d_alat;
            unitcellZ[28]=5.16935*d_alat;
            unitcellX[29]=0.57735*d_alat;
            unitcellY[29]=0*d_alat;
            unitcellZ[29]=5.4008*d_alat;
            unitcellX[30]=2.38419e-07*d_alat;
            unitcellY[30]=3.34075*d_alat;
            unitcellZ[30]=5.32367*d_alat;
            unitcellX[31]=1.19209e-07*d_alat;
            unitcellY[31]=2.6726*d_alat;
            unitcellZ[31]=5.55512*d_alat;
            unitcellX[32]=0.57735*d_alat;
            unitcellY[32]=2.53897*d_alat;
            unitcellZ[32]=5.16936*d_alat;
            unitcellX[33]=1.19209e-07*d_alat;
            unitcellY[33]=2.80623*d_alat;
            unitcellZ[33]=4.86074*d_alat;
            unitcellX[34]=0*d_alat;
            unitcellY[34]=2.13808*d_alat;
            unitcellZ[34]=5.09219*d_alat;
            unitcellX[35]=0.57735*d_alat;
            unitcellY[35]=1.87082*d_alat;
            unitcellZ[35]=5.40082*d_alat;
            unitcellX[36]=1.1547*d_alat;
            unitcellY[36]=1.73719*d_alat;
            unitcellZ[36]=5.01507*d_alat;
            unitcellX[37]=0.57735*d_alat;
            unitcellY[37]=2.00445*d_alat;
            unitcellZ[37]=4.70644*d_alat;
            unitcellX[38]=0.57735*d_alat;
            unitcellY[38]=1.3363*d_alat;
            unitcellZ[38]=4.93789*d_alat;
            unitcellX[39]=1.1547*d_alat;
            unitcellY[39]=1.06904*d_alat;
            unitcellZ[39]=5.24652*d_alat;
            unitcellX[40]=1.1547*d_alat;
            unitcellY[40]=1.20267*d_alat;
            unitcellZ[40]=4.55214*d_alat;
            unitcellX[41]=1.1547*d_alat;
            unitcellY[41]=0.53452*d_alat;
            unitcellZ[41]=4.78359*d_alat;
            unitcellX[42]=2.38419e-07*d_alat;
            unitcellY[42]=3.47438*d_alat;
            unitcellZ[42]=4.62928*d_alat;
            unitcellX[43]=0.57735*d_alat;
            unitcellY[43]=3.20712*d_alat;
            unitcellZ[43]=4.93791*d_alat;
            unitcellX[44]=1.1547*d_alat;
            unitcellY[44]=3.07349*d_alat;
            unitcellZ[44]=4.55216*d_alat;
            unitcellX[45]=0.57735*d_alat;
            unitcellY[45]=3.34075*d_alat;
            unitcellZ[45]=4.24353*d_alat;
            unitcellX[46]=0.57735*d_alat;
            unitcellY[46]=2.6726*d_alat;
            unitcellZ[46]=4.47498*d_alat;
            unitcellX[47]=1.1547*d_alat;
            unitcellY[47]=2.40534*d_alat;
            unitcellZ[47]=4.78361*d_alat;
            unitcellX[48]=1.1547*d_alat;
            unitcellY[48]=2.53897*d_alat;
            unitcellZ[48]=4.08923*d_alat;
            unitcellX[49]=1.1547*d_alat;
            unitcellY[49]=1.87082*d_alat;
            unitcellZ[49]=4.32068*d_alat;
            unitcellX[50]=1.1547*d_alat;
            unitcellY[50]=3.20712*d_alat;
            unitcellZ[50]=3.85777*d_alat;
            unitcellX[51]=-5.96046e-08*d_alat;
            unitcellY[51]=0.93541*d_alat;
            unitcellZ[51]=4.86071*d_alat;
            unitcellX[52]=0*d_alat;
            unitcellY[52]=0.26726*d_alat;
            unitcellZ[52]=5.09217*d_alat;
            unitcellX[53]=0.57735*d_alat;
            unitcellY[53]=0.13363*d_alat;
            unitcellZ[53]=4.70641*d_alat;
            unitcellX[54]=-1.19209e-07*d_alat;
            unitcellY[54]=0.40089*d_alat;
            unitcellZ[54]=4.39779*d_alat;
            unitcellX[55]=0*d_alat;
            unitcellY[55]=2.27171*d_alat;
            unitcellZ[55]=4.3978*d_alat;
            unitcellX[56]=0*d_alat;
            unitcellY[56]=1.60356*d_alat;
            unitcellZ[56]=4.62926*d_alat;
            unitcellX[57]=0.57735*d_alat;
            unitcellY[57]=1.46993*d_alat;
            unitcellZ[57]=4.24351*d_alat;
            unitcellX[58]=0*d_alat;
            unitcellY[58]=1.73719*d_alat;
            unitcellZ[58]=3.93488*d_alat;
            unitcellX[59]=-5.96046e-08*d_alat;
            unitcellY[59]=1.06904*d_alat;
            unitcellZ[59]=4.16633*d_alat;
            unitcellX[60]=0.57735*d_alat;
            unitcellY[60]=0.80178*d_alat;
            unitcellZ[60]=4.47496*d_alat;
            unitcellX[61]=1.1547*d_alat;
            unitcellY[61]=0.66815*d_alat;
            unitcellZ[61]=4.0892*d_alat;
            unitcellX[62]=0.57735*d_alat;
            unitcellY[62]=0.93541*d_alat;
            unitcellZ[62]=3.78058*d_alat;
            unitcellX[63]=0.57735*d_alat;
            unitcellY[63]=0.26726*d_alat;
            unitcellZ[63]=4.01203*d_alat;
            unitcellX[64]=1.1547*d_alat;
            unitcellY[64]=0*d_alat;
            unitcellZ[64]=4.32066*d_alat;
            unitcellX[65]=1.1547*d_alat;
            unitcellY[65]=0.13363*d_alat;
            unitcellZ[65]=3.62628*d_alat;
            unitcellX[66]=0*d_alat;
            unitcellY[66]=3.60801*d_alat;
            unitcellZ[66]=3.9349*d_alat;
            unitcellX[67]=0*d_alat;
            unitcellY[67]=2.93986*d_alat;
            unitcellZ[67]=4.16635*d_alat;
            unitcellX[68]=0.57735*d_alat;
            unitcellY[68]=2.80623*d_alat;
            unitcellZ[68]=3.7806*d_alat;
            unitcellX[69]=-1.19209e-07*d_alat;
            unitcellY[69]=3.07349*d_alat;
            unitcellZ[69]=3.47196*d_alat;
            unitcellX[70]=-1.19209e-07*d_alat;
            unitcellY[70]=2.40534*d_alat;
            unitcellZ[70]=3.70342*d_alat;
            unitcellX[71]=0.57735*d_alat;
            unitcellY[71]=2.13808*d_alat;
            unitcellZ[71]=4.01205*d_alat;
            unitcellX[72]=1.1547*d_alat;
            unitcellY[72]=2.00445*d_alat;
            unitcellZ[72]=3.6263*d_alat;
            unitcellX[73]=0.57735*d_alat;
            unitcellY[73]=2.27171*d_alat;
            unitcellZ[73]=3.31766*d_alat;
            unitcellX[74]=0.57735*d_alat;
            unitcellY[74]=1.60356*d_alat;
            unitcellZ[74]=3.54912*d_alat;
            unitcellX[75]=1.1547*d_alat;
            unitcellY[75]=1.3363*d_alat;
            unitcellZ[75]=3.85775*d_alat;
            unitcellX[76]=1.1547*d_alat;
            unitcellY[76]=1.46993*d_alat;
            unitcellZ[76]=3.16336*d_alat;
            unitcellX[77]=1.1547*d_alat;
            unitcellY[77]=0.80178*d_alat;
            unitcellZ[77]=3.39482*d_alat;
            unitcellX[78]=0.57735*d_alat;
            unitcellY[78]=3.47438*d_alat;
            unitcellZ[78]=3.54914*d_alat;
            unitcellX[79]=1.1547*d_alat;
            unitcellY[79]=3.34075*d_alat;
            unitcellZ[79]=3.16339*d_alat;
            unitcellX[80]=0.57735*d_alat;
            unitcellY[80]=3.60801*d_alat;
            unitcellZ[80]=2.85475*d_alat;
            unitcellX[81]=0.57735*d_alat;
            unitcellY[81]=2.93986*d_alat;
            unitcellZ[81]=3.08621*d_alat;
            unitcellX[82]=1.1547*d_alat;
            unitcellY[82]=2.6726*d_alat;
            unitcellZ[82]=3.39484*d_alat;
            unitcellX[83]=1.1547*d_alat;
            unitcellY[83]=2.80623*d_alat;
            unitcellZ[83]=2.70045*d_alat;
            unitcellX[84]=1.1547*d_alat;
            unitcellY[84]=2.13808*d_alat;
            unitcellZ[84]=2.93191*d_alat;
            unitcellX[85]=1.1547*d_alat;
            unitcellY[85]=3.47438*d_alat;
            unitcellZ[85]=2.469*d_alat;
            unitcellX[86]=0*d_alat;
            unitcellY[86]=1.20267*d_alat;
            unitcellZ[86]=3.47195*d_alat;
            unitcellX[87]=0*d_alat;
            unitcellY[87]=0.53452*d_alat;
            unitcellZ[87]=3.7034*d_alat;
            unitcellX[88]=0.57735*d_alat;
            unitcellY[88]=0.40089*d_alat;
            unitcellZ[88]=3.31765*d_alat;
            unitcellX[89]=5.96046e-08*d_alat;
            unitcellY[89]=0.66815*d_alat;
            unitcellZ[89]=3.00901*d_alat;
            unitcellX[90]=5.96046e-08*d_alat;
            unitcellY[90]=0*d_alat;
            unitcellZ[90]=3.24047*d_alat;
            unitcellX[91]=0*d_alat;
            unitcellY[91]=2.53897*d_alat;
            unitcellZ[91]=3.00904*d_alat;
            unitcellX[92]=0*d_alat;
            unitcellY[92]=1.87082*d_alat;
            unitcellZ[92]=3.24049*d_alat;
            unitcellX[93]=0.57735*d_alat;
            unitcellY[93]=1.73719*d_alat;
            unitcellZ[93]=2.85474*d_alat;
            unitcellX[94]=0*d_alat;
            unitcellY[94]=2.00445*d_alat;
            unitcellZ[94]=2.5461*d_alat;
            unitcellX[95]=5.96046e-08*d_alat;
            unitcellY[95]=1.3363*d_alat;
            unitcellZ[95]=2.77756*d_alat;
            unitcellX[96]=0.57735*d_alat;
            unitcellY[96]=1.06904*d_alat;
            unitcellZ[96]=3.08619*d_alat;
            unitcellX[97]=1.1547*d_alat;
            unitcellY[97]=0.93541*d_alat;
            unitcellZ[97]=2.70044*d_alat;
            unitcellX[98]=0.57735*d_alat;
            unitcellY[98]=1.20267*d_alat;
            unitcellZ[98]=2.3918*d_alat;
            unitcellX[99]=0.57735*d_alat;
            unitcellY[99]=0.53452*d_alat;
            unitcellZ[99]=2.62326*d_alat;
            unitcellX[100]=1.1547*d_alat;
            unitcellY[100]=0.26726*d_alat;
            unitcellZ[100]=2.93189*d_alat;
            unitcellX[101]=1.1547*d_alat;
            unitcellY[101]=0.40089*d_alat;
            unitcellZ[101]=2.2375*d_alat;
            unitcellX[102]=0*d_alat;
            unitcellY[102]=3.20712*d_alat;
            unitcellZ[102]=2.77758*d_alat;
            unitcellX[103]=0.57735*d_alat;
            unitcellY[103]=3.07349*d_alat;
            unitcellZ[103]=2.39183*d_alat;
            unitcellX[104]=1.19209e-07*d_alat;
            unitcellY[104]=3.34075*d_alat;
            unitcellZ[104]=2.08319*d_alat;
            unitcellX[105]=0*d_alat;
            unitcellY[105]=2.6726*d_alat;
            unitcellZ[105]=2.31465*d_alat;
            unitcellX[106]=0.57735*d_alat;
            unitcellY[106]=2.40534*d_alat;
            unitcellZ[106]=2.62328*d_alat;
            unitcellX[107]=1.1547*d_alat;
            unitcellY[107]=2.27171*d_alat;
            unitcellZ[107]=2.23753*d_alat;
            unitcellX[108]=0.57735*d_alat;
            unitcellY[108]=2.53897*d_alat;
            unitcellZ[108]=1.92889*d_alat;
            unitcellX[109]=0.57735*d_alat;
            unitcellY[109]=1.87082*d_alat;
            unitcellZ[109]=2.16035*d_alat;
            unitcellX[110]=1.1547*d_alat;
            unitcellY[110]=1.60356*d_alat;
            unitcellZ[110]=2.46898*d_alat;
            unitcellX[111]=1.1547*d_alat;
            unitcellY[111]=1.73719*d_alat;
            unitcellZ[111]=1.7746*d_alat;
            unitcellX[112]=1.1547*d_alat;
            unitcellY[112]=1.06904*d_alat;
            unitcellZ[112]=2.00605*d_alat;
            unitcellX[113]=1.1547*d_alat;
            unitcellY[113]=3.60801*d_alat;
            unitcellZ[113]=1.77462*d_alat;
            unitcellX[114]=0.57735*d_alat;
            unitcellY[114]=3.20712*d_alat;
            unitcellZ[114]=1.69744*d_alat;
            unitcellX[115]=1.1547*d_alat;
            unitcellY[115]=2.93986*d_alat;
            unitcellZ[115]=2.00607*d_alat;
            unitcellX[116]=1.1547*d_alat;
            unitcellY[116]=3.07349*d_alat;
            unitcellZ[116]=1.31169*d_alat;
            unitcellX[117]=1.1547*d_alat;
            unitcellY[117]=2.40534*d_alat;
            unitcellZ[117]=1.54314*d_alat;
            unitcellX[118]=0*d_alat;
            unitcellY[118]=0.13363*d_alat;
            unitcellZ[118]=2.54609*d_alat;
            unitcellX[119]=0*d_alat;
            unitcellY[119]=1.46993*d_alat;
            unitcellZ[119]=2.08317*d_alat;
            unitcellX[120]=0*d_alat;
            unitcellY[120]=0.80178*d_alat;
            unitcellZ[120]=2.31463*d_alat;
            unitcellX[121]=0.57735*d_alat;
            unitcellY[121]=0.66815*d_alat;
            unitcellZ[121]=1.92887*d_alat;
            unitcellX[122]=0*d_alat;
            unitcellY[122]=0.93541*d_alat;
            unitcellZ[122]=1.62024*d_alat;
            unitcellX[123]=0*d_alat;
            unitcellY[123]=0.26726*d_alat;
            unitcellZ[123]=1.8517*d_alat;
            unitcellX[124]=0.57735*d_alat;
            unitcellY[124]=0*d_alat;
            unitcellZ[124]=2.16033*d_alat;
            unitcellX[125]=0.57735*d_alat;
            unitcellY[125]=0.13363*d_alat;
            unitcellZ[125]=1.46595*d_alat;
            unitcellX[126]=0*d_alat;
            unitcellY[126]=2.80623*d_alat;
            unitcellZ[126]=1.62027*d_alat;
            unitcellX[127]=0*d_alat;
            unitcellY[127]=2.13808*d_alat;
            unitcellZ[127]=1.85172*d_alat;
            unitcellX[128]=0.57735*d_alat;
            unitcellY[128]=2.00445*d_alat;
            unitcellZ[128]=1.46597*d_alat;
            unitcellX[129]=0*d_alat;
            unitcellY[129]=2.27171*d_alat;
            unitcellZ[129]=1.15734*d_alat;
            unitcellX[130]=0*d_alat;
            unitcellY[130]=1.60356*d_alat;
            unitcellZ[130]=1.38879*d_alat;
            unitcellX[131]=0.57735*d_alat;
            unitcellY[131]=1.3363*d_alat;
            unitcellZ[131]=1.69742*d_alat;
            unitcellX[132]=1.1547*d_alat;
            unitcellY[132]=1.20267*d_alat;
            unitcellZ[132]=1.31167*d_alat;
            unitcellX[133]=0.57735*d_alat;
            unitcellY[133]=1.46993*d_alat;
            unitcellZ[133]=1.00304*d_alat;
            unitcellX[134]=0.57735*d_alat;
            unitcellY[134]=0.80178*d_alat;
            unitcellZ[134]=1.23449*d_alat;
            unitcellX[135]=1.1547*d_alat;
            unitcellY[135]=0.53452*d_alat;
            unitcellZ[135]=1.54312*d_alat;
            unitcellX[136]=1.1547*d_alat;
            unitcellY[136]=0.66815*d_alat;
            unitcellZ[136]=0.848735*d_alat;
            unitcellX[137]=1.1547*d_alat;
            unitcellY[137]=0*d_alat;
            unitcellZ[137]=1.08019*d_alat;
            unitcellX[138]=0*d_alat;
            unitcellY[138]=3.47438*d_alat;
            unitcellZ[138]=1.38881*d_alat;
            unitcellX[139]=0.57735*d_alat;
            unitcellY[139]=3.34075*d_alat;
            unitcellZ[139]=1.00306*d_alat;
            unitcellX[140]=0*d_alat;
            unitcellY[140]=3.60801*d_alat;
            unitcellZ[140]=0.694425*d_alat;
            unitcellX[141]=0*d_alat;
            unitcellY[141]=2.93986*d_alat;
            unitcellZ[141]=0.92588*d_alat;
            unitcellX[142]=0.57735*d_alat;
            unitcellY[142]=2.6726*d_alat;
            unitcellZ[142]=1.23451*d_alat;
            unitcellX[143]=1.1547*d_alat;
            unitcellY[143]=2.53897*d_alat;
            unitcellZ[143]=0.848755*d_alat;
            unitcellX[144]=0.57735*d_alat;
            unitcellY[144]=2.80623*d_alat;
            unitcellZ[144]=0.540125*d_alat;
            unitcellX[145]=0.57735*d_alat;
            unitcellY[145]=2.13808*d_alat;
            unitcellZ[145]=0.77158*d_alat;
            unitcellX[146]=1.1547*d_alat;
            unitcellY[146]=1.87082*d_alat;
            unitcellZ[146]=1.08021*d_alat;
            unitcellX[147]=1.1547*d_alat;
            unitcellY[147]=2.00445*d_alat;
            unitcellZ[147]=0.385825*d_alat;
            unitcellX[148]=1.1547*d_alat;
            unitcellY[148]=1.3363*d_alat;
            unitcellZ[148]=0.61728*d_alat;
            unitcellX[149]=0.57735*d_alat;
            unitcellY[149]=3.47438*d_alat;
            unitcellZ[149]=0.30867*d_alat;
            unitcellX[150]=1.1547*d_alat;
            unitcellY[150]=3.20712*d_alat;
            unitcellZ[150]=0.6173*d_alat;
            unitcellX[151]=1.1547*d_alat;
            unitcellY[151]=3.34075*d_alat;
            unitcellZ[151]=-0.0770848*d_alat;
            unitcellX[152]=1.1547*d_alat;
            unitcellY[152]=2.6726*d_alat;
            unitcellZ[152]=0.15437*d_alat;
            unitcellX[153]=0*d_alat;
            unitcellY[153]=0.40089*d_alat;
            unitcellZ[153]=1.15732*d_alat;
            unitcellX[154]=0*d_alat;
            unitcellY[154]=1.73719*d_alat;
            unitcellZ[154]=0.694405*d_alat;
            unitcellX[155]=0*d_alat;
            unitcellY[155]=1.06904*d_alat;
            unitcellZ[155]=0.92586*d_alat;
            unitcellX[156]=0.57735*d_alat;
            unitcellY[156]=0.93541*d_alat;
            unitcellZ[156]=0.540105*d_alat;
            unitcellX[157]=0*d_alat;
            unitcellY[157]=1.20267*d_alat;
            unitcellZ[157]=0.231475*d_alat;
            unitcellX[158]=0*d_alat;
            unitcellY[158]=0.53452*d_alat;
            unitcellZ[158]=0.46293*d_alat;
            unitcellX[159]=0.57735*d_alat;
            unitcellY[159]=0.26726*d_alat;
            unitcellZ[159]=0.77156*d_alat;
            unitcellX[160]=1.1547*d_alat;
            unitcellY[160]=0.13363*d_alat;
            unitcellZ[160]=0.385805*d_alat;
            unitcellX[161]=0.57735*d_alat;
            unitcellY[161]=0.40089*d_alat;
            unitcellZ[161]=0.077175*d_alat;
            unitcellX[162]=0*d_alat;
            unitcellY[162]=3.07349*d_alat;
            unitcellZ[162]=0.231495*d_alat;
            unitcellX[163]=0*d_alat;
            unitcellY[163]=2.40534*d_alat;
            unitcellZ[163]=0.46295*d_alat;
            unitcellX[164]=0.57735*d_alat;
            unitcellY[164]=2.27171*d_alat;
            unitcellZ[164]=0.077195*d_alat;
            unitcellX[165]=0*d_alat;
            unitcellY[165]=1.87082*d_alat;
            unitcellZ[165]=2.00421e-05*d_alat;
            unitcellX[166]=0.57735*d_alat;
            unitcellY[166]=1.60356*d_alat;
            unitcellZ[166]=0.30865*d_alat;
            unitcellX[167]=1.1547*d_alat;
            unitcellY[167]=1.46993*d_alat;
            unitcellZ[167]=-0.077105*d_alat;
            unitcellX[168]=1.1547*d_alat;
            unitcellY[168]=0.80178*d_alat;
            unitcellZ[168]=0.15435*d_alat;

            pre_make_xtal(0);
        }

        if(tokens[1]=="hcp")
        {   xtalf_rname = "gen.hcpgsfe3.atoms";
            xtalf_wname = "generic.hcpgsfe3.atoms";
            mockxtal = "HCP";

            enx = 4; eny = 2; enz = 2;
            no_atoms = 168*enx*eny*enz;
            cord_no = 168;

            xfact = 1.731999993;
            yfact = 3.741650105;
            zfact = 6.48074007;

            pre_make_xtal(0);
        }
    }
}


//---------------------------------------------------------------------------------------------------------------------------
void MD::pre_make_xtal(int a)
{   //make crystal structure based on alat
    //defining the necessary vectors
    vector <double> posx,posy,posz,num1,num2;
    int i,j,k,l,m=0,co=1;
    

    //write to a file
    ofstream xtal2(xtalf_wname.c_str());

    xtal2<<"# Atoms generated by XtalGen ©";
    xtal2<<endl<<endl<<no_atoms<<" "<<"atoms"<<endl<<endl;
    xtal2<<"1 atom types"<<endl<<endl;

    xtal2<<"0    "<<enx*d_alat*xfact<<"    xlo xhi"<<endl;
    xtal2<<"0    "<<eny*d_alat*yfact<<"    ylo yhi"<<endl;
    xtal2<<"0    "<<enz*d_alat*zfact<<"    zlo zhi"<<endl;

    if(a==1)
    {
        xtal2<<"0.0 0.0 0.0 xy xz yz"<<endl;
    }

    xtal2<<endl<<"Atoms"<<endl<<endl;

    for(i=1;i<=enx;i++)
    {
        for(j=1;j<=eny;j++)
        {
             for(k=1;k<=enz;k++)
             {
                  for(l=1;l<=cord_no;l++)
                  {
                        m++;
                        xtal2<<fixed<<setprecision(5)<<co<<" 1 "<<(unitcellX[l]+(d_alat*xfact*((double)i-1)))<<" ";
                        xtal2<<fixed<<setprecision(5)<<(unitcellY[l]+(d_alat*yfact*((double)j-1)))<<" ";
                        xtal2<<fixed<<setprecision(5)<<(unitcellZ[l]+(d_alat*zfact*((double)k-1)))<<endl;
                        co++;
                    }
                }
            }
        }

    xtal2.close();

    //copying the file to required folder
    ofstream xtal3("cpatoms.sh");

    xtal3<<"cd "<<mockxtal.c_str()<<"/"<<endl;
    xtal3<<"rm "<<xtalf_wname.c_str()<<endl;
    xtal3<<"cd ../"<<endl;
    xtal3<<"cp "<<xtalf_wname.c_str()<<" "<<mockxtal.c_str()<<"/"<<endl;

    xtal3.close();

    //make it runnable and running
    system("./cpatoms1.sh");

}
