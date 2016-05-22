
#include "mdcalc_header.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>

void MD::read_result(){

	ifstream fin("result.out");
	double dummy;
	vector <float> result;

	while(fin>>dummy){
		result.push_back(dummy);
	}

	if (tokens[0]=="cohesiveeng")
    {  
    	cout<<endl<<endl;
    	cout<<"------------------------------------"<<endl;
    	cout<<"Cohesive Energy: "<<result[0]<<" eV"<<endl;
    	cout<<"------------------------------------"<<endl;
    }

    else if ((tokens[0]=="elasticconst")&&((tokens[1]=="bcc")||(tokens[1]=="fcc")))
    {  
    	cout<<endl<<endl;
    	cout<<"------------------------------------"<<endl;
    	cout<<"C11: "<<result[0]<<" GPa"<<endl;
    	cout<<"C12: "<<result[1]<<" GPa"<<endl;
    	cout<<"C14: "<<result[2]<<" GPa"<<endl;
    	cout<<"------------------------------------"<<endl;
    }

    else if ((tokens[0]=="elasticconst")&&(tokens[1]=="hcp"))
    {  
    	cout<<endl<<endl;
    	cout<<"------------------------------------"<<endl;
    	cout<<"C11: "<<result[0]<<" GPa"<<endl;
    	cout<<"C12: "<<result[1]<<" GPa"<<endl;
    	cout<<"C14: "<<result[2]<<" GPa"<<endl;
    	cout<<"C13: "<<result[3]<<" GPa"<<endl;
    	cout<<"C33: "<<result[4]<<" GPa"<<endl;
    	cout<<"------------------------------------"<<endl;
    }

    else if (tokens[0]=="bulkmod")
    {  
    	cout<<endl<<endl;
    	cout<<"------------------------------------"<<endl;
    	cout<<"Bulk Modulus: "<<result[0]<<" GPa"<<endl;
    	cout<<"------------------------------------"<<endl;
    }


}