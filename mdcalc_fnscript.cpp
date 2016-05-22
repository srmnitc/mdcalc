

#include "mdcalc_header.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>


//--------------------------------------------------------------------------------------------------------------------------------
void MD::make_script(){
	
	ofstream shell1("run.sh");

    shell1<<"rm log.lammps"<<endl;

    //running whichever calculation is selected
    //checking if cohesive energy calculation is selected

    if (tokens[0]=="cohesiveeng")
    {  
    	if(tokens[1]=="bcc")
                { 
                  shell1<<"cd BCC/"<<endl;
                  shell1<<"./1cohesive.sh"<<endl;
                  shell1<<"mv log.lammps ../log.lammps.cohesive.bcc"<<endl;
                  shell1<<"mv out_cohesive ../result.out"<<endl;
                  
                }
        else if(tokens[1]=="fcc")
                { shell1<<"cd FCC/"<<endl;
                  shell1<<"./1cohesive.sh"<<endl;
                  shell1<<"mv log.lammps ../log.lammps.cohesive.fcc"<<endl;
                  shell1<<"mv out_cohesive ../result.out"<<endl;
                  
                }
        else if(tokens[1]=="hcp")
                { shell1<<"cd HCP/"<<endl;
                  shell1<<"./1cohesive.sh"<<endl;
                  shell1<<"mv log.lammps ../log.lammps.cohesive.hcp"<<endl;
                  shell1<<"mv out_cohesive ../result.out"<<endl;
                  
                }
    }

    if (tokens[0]=="elasticconst")
    {  if(tokens[1]=="bcc")
                { shell1<<"cd BCC/"<<endl;
                  shell1<<"./3strain.sh"<<endl;
                  shell1<<"mv log.lammps ../log.lammps.strain.bcc"<<endl;
                  shell1<<"mv out_elastic ../result.out"<<endl;
               

                }
        else if(tokens[1]=="fcc")
                { shell1<<"cd FCC/"<<endl;
                  shell1<<"./3strain.sh"<<endl;
                  shell1<<"mv log.lammps ../log.lammps.strain.fcc"<<endl;
                  shell1<<"mv out_elastic ../result.out"<<endl;
              

                }
        else if(tokens[1]=="hcp")
                { shell1<<"cd HCP/"<<endl;
                  shell1<<"./3strain.sh"<<endl;
                  shell1<<"mv log.lammps ../log.lammps.strain.hcp"<<endl;
                  shell1<<"mv out_elastic ../result.out"<<endl;
                

                }
    }

    if (tokens[0]=="bulkmod")
    {  if(tokens[1]=="bcc")
                { shell1<<"cd BCC/"<<endl;
                  shell1<<"./2bulk.sh"<<endl;
                  shell1<<"mv log.lammps ../log.lammps.bulk.bcc"<<endl;
                  shell1<<"mv out_bulk ../result.out"<<endl;
                  

                }
        else if(tokens[1]=="fcc")
                { shell1<<"cd FCC/"<<endl;
                  shell1<<"./2bulk.sh"<<endl;
                  shell1<<"mv log.lammps ../log.lammps.bulk.fcc"<<endl;
                  shell1<<"mv out_bulk ../result.out"<<endl;
                  

                }
        else if(tokens[1]=="hcp")
                { shell1<<"cd HCP/"<<endl;
                  shell1<<"./2bulk.sh"<<endl;
                  shell1<<"mv log.lammps ../log.lammps.bulk.hcp"<<endl;
                  shell1<<"mv out_bulk ../result.out"<<endl;
                  

                }
    }

}
       