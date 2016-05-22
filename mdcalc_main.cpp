

#include "mdcalc_header.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

int main(){


	string input_string;
	vector <string> tokenss;
	string buf;
	int indcr = 0;
	
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"                    Interface for LAMMPS                    "<<endl;
	cout<<"------------------------------------------------------------"<<endl;
	
	cout<<endl<<"md > ";
	getline(cin,input_string);
	stringstream ss(input_string);

	//split the string
	while (!ss.eof()){

		ss >> buf;
		tokenss.push_back(buf);
	}


	MD md1;
	md1.make_vector();

	if(tokenss[0]=="exit"){
		indcr = 1;
	}
	//keep the program running
	cout<<indcr<<endl;

	while (indcr!=1){

		//copy the vector
		md1.copy_vector(tokenss);

		//input style is
		//calc_keyword xtal_structure xtal_latticeconst potential_name

		//make crystal structure
		md1.calc_lat();
		md1.make_xtal();
		md1.make_script();
		md1.make_potential();
		system("./runfin.sh");
        system("./runfin2.sh");
        md1.read_result();

		
        cout<<endl<<endl;
		cout<<endl<<"md > ";
		getline(cin,input_string);
		stringstream ss(input_string);
		tokenss.clear();
		//split the string
		while (!ss.eof()){

			ss >> buf;
			tokenss.push_back(buf);
		}
	}

}