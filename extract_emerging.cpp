///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                       //
//                                                                                                                       //
// This small c++ program opens the pe-trajectories                                                                      //
// Then it searches (atm only) for the backscattered particles, reads the last two coordinates and the energy            //
// The last two coordinates are stored in a file with the energy                                                         //
// With this information the vector (direction) of the emerging(backscattered) particles can be calculated               //
//                                                                                                                       //
// This can be used to get the angular distribution of the emerging particles with the corresponding energies            //
//                                                                                                                       //
// compiled with 'g++ extract_emerging.cpp' and executed using './a.out'                                                 //
//                                                                                                                       //
//                                                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <math.h>
#include <ctime>

#define _USE_MATH_DEFINES

#ifndef PI
#define PI	3.14159265358979323846f
#endif

using namespace std;	// don't do this in the future!

// signum function
int sgn(double x)
{
	if (x > 0.0L)
		return 1.0L;
	else if (x < 0.0L)  
		return -1.0L;
	else  
		return 0.0L;
}

// Converts a string into a float, integer, etc.
template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

int main()
{
	
	//////////////////////////////////////////////////////////////////////
	// dont define all variables here.									//
	// Define them, where they are first used for a better oversight.	//
	// Some may be kept here at the start.								//
	//////////////////////////////////////////////////////////////////////

								///////////////////////////////////VARIABLES//////////////////////////////////////////////////
	string zeile;				// Bufferstring - the File is read line by line into the string "zeile"                     //
	string cutstr;				//                                                                                          //
	string cdsandenrgy[2][5];	// 2d array -> storing the energ and x,y,z-coordinates                                      //
	string temparray[6];		// temparray[k] k=0 ... TRAJ: Index (starting at 1)                                         //
								//              k=1 ... KPAR: Type of particle (electron=1, photon=2, positron=3)           //
								//              k=2 ... PARENT: Index of the parent trajectory (primary=0)                  //
								//              k=3 ... ICOL: Type of collision that created this trajectory (primary=0)    //
								//              k=4 ... EXIT: Exit type (backscattered=1, transmitted=2, absorbed=3)        //
	string file_in;				// inputfile and path (e.g. folder1/folder2/pe-trajectories.dat                             //
	string file_out1;			// outputfile 1 and path (e.g. folder1/folder2/pe-trajectories.dat                          //
	string file_out2;			// outputfile 2 and path (e.g. folder1/folder2/pe-trajectories.dat                          //
	string file_out3;			// outputfile 3 and path (e.g. folder1/folder2/pe-trajectories.dat                          //
	string file_out4;			// outputfile 4 and path (e.g. folder1/folder2/pe-trajectories.dat                          //
	string file_out5;			// outputfile 5 and path (e.g. folder1/folder2/pe-trajectories.dat                          //
	string file_out6;			// outputfile 6 and path (e.g. folder1/folder2/pe-trajectories.dat                          //
	string file_out7;			// outputfile 7                                                                             //
	string file_out8;			// outputfile 8                                                                             //
								//                                                                                          //
	string file_out10;			// outputfile 10                                                                            //
	string file_out11;			// outputfile 11                                                                            //
	string file_out12;			// outputfile 12                                                                            //
	string file_out13;			// outputfile 13                                                                            //
	string file_out14;			// outputfile 14                                                                            //
								//                                                                                          //
	int int_1 = 0;				//                                                                                          //
	int control1 = 0;			// Control_int for counting backscattered particles                                         //
	int control2 = 0;			// Control_int for counting transmitted particles                                           //
	int control3 = 0;			// Control_int for counting absorbed particles                                              //
	int flag1 = 1;				// Controlflag                                                                              //
	int ct1 = 0;				// Countvariable                                                                            //
								//                                                                                          //   
	bool dflag;					// bool-flag(for check of file-opening)                                                     //
								//                                                                                          //   
	long long energy = 0;		// maxrange = 9223372036854775807LL !!!!!!!!!!!!!!!!                                        //
	long long energy_t = 0;		//                                                                                          //
								//                                                                                          //   
	float e1, e2, e3 = 0;		// Exponents (Stores 15 from like 1E+15)                                                    //
	float nrg = 0;				// Energy                                                                                   //
	float norm = 0;				// Length of the vectors                                                                    //
	float string_t;				// temporary String                                                                         //
	float energy_p = 0.0;		// highest Energy                                                                           //
	float control_particle;		// controlparticle storage value                                                            //
								//                                                                                          //   
	double tilt_angle = -75.0;	// angle to tilt the vectors along the y-axis                                               //
	double x1, x2, x3 = 0.0;	// Coordinates Point 1                                                                      //
	double y1, y2, y3 = 0.0;	// Coordinates Point 2                                                                      //
	double d_v1, d_v2 = 0.0;	// dummy-variables (temp. storage)                                                          //
								//                                                                                          //
	time_t rawtime;				// actual time and date                                                                     //   
								///////////////////////////////////VARIABLES//////////////////////////////////////////////////

	time ( &rawtime );			// actual date and time //

	std::vector<double> kx1, kx2, kx3, ky1, ky2, ky3, e, e_cut, e_cut_negative;  	// array to store the last two x/y/z-coordinates and the energy //
																					// exiting vector(direction) can be calculated with that info   //

	file_in = "pe-trajectories.dat";										// inputfilename ('../../pe-trajectories.dat')                                   //

	file_out1 = "extract_emerging/backscatang_cartesian.dat";				// Outputfilename for Backscattered Energies over exiting angle (cartesian)      //
	file_out2 = "extract_emerging/binned_energy_spherical.dat";				// Outputfilename for Backscattered Energies over exiting angle (spherical)      //
	file_out3 = "extract_emerging/binned_probability_spherical.dat";		// Outputfilename for Backscattered probabilities over exiting angle (spherical) //
	file_out4 = "extract_emerging/energy_probability.dat";					// Outputfilename for energies (bins) over probability                           //
// 	file_out5 = "extract_emerging/angle_probability.dat";					// Outputfilename for probability (bins) over angle                              //
	file_out6 = "extract_emerging/info.dat";								// Outputfilename for the general information (number of particles, energies...) //
	file_out7 = "extract_emerging/energy_raw_data.dat";						// Outputfilename for backscattered energies (not binned -> ("raw data") )       //
// 	file_out8 = "extract_emerging/energy_angle.dat";						// Outputfilename - angle theta over energies (angle phi is summarized)          //

	file_out10 = "extract_emerging/probability_summarized.dat";				// binned probability summarized over angle phi  -> theta over prob. data        //
	file_out11 = "extract_emerging/energy_summarized.dat";					// binned probability summarized over angle phi  -> theta over prob. data        //

	file_out12 = "extract_emerging/energy_probability_pure.dat";			// hits over bins (just hits nothing done more with that)        //

	file_out13 = "extract_emerging/energy_probability_cut.dat";				// binned probability summarized over angle phi  -> theta over prob. data        //
	file_out14 = "extract_emerging/energy_probability_cut_negative.dat";	// binned probability summarized over angle phi  -> theta over prob. data        //


	fstream datei(file_in.c_str(), ios::in);								// open an fstream (and the file) to read                                        //
	fstream bckang(file_out1.c_str(), ios::out);
	fstream bckpol(file_out2.c_str(), ios::out);
	fstream bckprob(file_out3.c_str(), ios::out);
	fstream nrgprob(file_out4.c_str(), ios::out);
// 	fstream bckprob2(file_out5.c_str(), ios::out);
	fstream infofile(file_out6.c_str(), ios::out);
	fstream rawenerg(file_out7.c_str(), ios::out);
// 	fstream energyangle(file_out8.c_str(), ios::out);
	fstream prob_theta_binned(file_out10.c_str(), ios::out);
	fstream energy_theta_binned(file_out11.c_str(), ios::out);
	fstream energy_prob_pure(file_out12.c_str(), ios::out);
	fstream energy_prob_cut(file_out13.c_str(), ios::out);
	fstream energy_prob_cut_neg(file_out14.c_str(), ios::out);

	dflag = (bool)getline(datei, zeile);

	bckang << "#  x-coordinate1a y-coodinate1a z-coordinate1a (all in cm)" << endl;
	bckang << "#  x-coordinate1b y-coodinate1b z-coordinate1b energy" << endl;   

	if(dflag == 1)	// dflag == 1 ... file has been opened without problems //
	{
		cout <<  "\033[40;0;34mstarting... \033[0m" << endl;
		cout << "file: '\033[40;0;33m" << file_in << "\033[0m' opened" << endl;
		cout << endl;

		cout <<  "\033[40;0;31mReading file & extracting information (exiting vectors, energies)... \033[0m" << endl;

		while ( !datei.eof() )
		{
			if (zeile[0] == '0' and (bool)getline(datei, zeile) != 0) // First char of string = '0' -> beginning of a new particle //
			{
				for (int h1 = 0; h1 < 5; h1++)
				{
					cutstr = zeile.substr(8, 7);  // cuts a part of the string out to get the header of a file //
					temparray[h1] = cutstr;
					getline(datei, zeile);
				}

				// "      1" -> trajectory of a backscattered particle
				if (temparray[4] == "      1")
				{
					// GET BACKSCATTERED ENERGY OVER ANGULAR INFORMATION //
					flag1 = 0;
					ct1 = 0;

					while(flag1 == 0)
					{
						ct1++;

						cdsandenrgy[0][0] = zeile.substr (1, 13);	// cut x-coordinate1 out of the string 'zeile' //
						cdsandenrgy[0][1] = zeile.substr (15, 13);	// cut y-coordinate1 out of the string 'zeile' //
						cdsandenrgy[0][2] = zeile.substr (29, 13);	// cut z-coordinate1 out of the string 'zeile' //
						cdsandenrgy[0][3] = zeile.substr (44, 12);	// cut energy1  out of the string 'zeile'      //

						getline(datei, zeile);

						if (zeile[0] == '0')
						{
							if (ct1 > 3)  // prevents errors (backscattering with only one collission in the trajectoryfile) //
							{
								from_string<double>(x1, std::string(cdsandenrgy[1][0]), std::dec);	// Point 1_x //
								from_string<double>(x2, std::string(cdsandenrgy[1][1]), std::dec);	// Point 1_y //
								from_string<double>(x3, std::string(cdsandenrgy[1][2]), std::dec);	// Point 1_z //

								from_string<float>(nrg, std::string(cdsandenrgy[1][3]), std::dec);	// Energy //

								// energy_p ... peakenergy (the highest energy found from a backscattered particle //
								if (nrg > energy_p)
								{
									energy_p = nrg;
								}

								from_string<double>(y1, std::string(cdsandenrgy[0][0].substr (0, 9)), std::dec);	// Point 2_x //
								from_string<double>(y2, std::string(cdsandenrgy[0][1].substr (0, 9)), std::dec);	// Point 2_y //
								from_string<double>(y3, std::string(cdsandenrgy[0][2].substr (0, 9)), std::dec);	// Point 2_z //

								from_string<float>(e1, std::string(cdsandenrgy[0][0].substr (11, 2)), std::dec);	// Exponents 2_x //
								from_string<float>(e2, std::string(cdsandenrgy[0][1].substr (11, 2)), std::dec);	// Exponents 2_y //
								from_string<float>(e3, std::string(cdsandenrgy[0][2].substr (11, 2)), std::dec);	// Exponents 2_z //

								if (e1 < 28 or e2 < 28 or e3 < 28)
								{
									cout << endl;
									cout << "error - please check the .dat file" << endl;
									cout << "the last entry of the coordinates of the backscattered particle is below 1E+28" << endl;
									cout << "if ignored this may result into a miscalculation" << endl;
									cout << "e1: " << e1 << " ; e2: " << e2 << " ; e3: " << e3 << endl;
								}

								y1 *= pow (10.0, e1-28);	// reduce exponent to prevent errors with large exponents (too big for float)
								y2 *= pow (10.0, e2-28);	// reduce exponent to prevent errors with large exponents (too big for float)
								y3 *= pow (10.0, e3-28);	// reduce exponent to prevent errors with large exponents (too big for float)

								// Normalize the Vectors //
								norm = sqrt( (e1 - y1) * (e1 - y1) + (e2 - y2) * (e2 - y2) + (e3 - y3) * (e3 - y3) ); // calculate length of the vector //

								y1 = (x1 + (y1 - x1) / norm); // y1 ... the "end-coordinate" of the normalized x-coordinate //
								y2 = (x2 + (y2 - x2) / norm); // y2 ... the "end-coordinate" of the normalized y-coordinate //
								y3 = (x3 + (y3 - x3) / norm); // y3 ... the "end-coordinate" of the normalized z-coordinate //
								// Normalize the Vectors //

								// rotate the vectors
								/*
								// x-Axis-rotation
								x2 = x2*cos(tilt_angle * (M_PI / 180)) - x3 * sin(tilt_angle * (M_PI / 180));
								x3 = x2*sin(tilt_angle * (M_PI / 180)) + x3 * cos(tilt_angle * (M_PI / 180));

								y2 = y2*cos(tilt_angle * (M_PI / 180)) - y3*sin(tilt_angle * (M_PI / 180));
								y3 = y2*sin(tilt_angle * (M_PI / 180)) + y3*cos(tilt_angle * (M_PI / 180));
								*/
								
								// y-Axis-rotation
								d_v1 =  x1 * cos(tilt_angle * (M_PI / 180)) + x3 * sin(tilt_angle * (M_PI / 180));
								d_v2 = -x1 * sin(tilt_angle * (M_PI / 180)) + x3 * cos(tilt_angle * (M_PI / 180));

								x1 = d_v1;
								x3 = d_v2;

								d_v1 =  y1 * cos(tilt_angle * (M_PI / 180)) + y3 * sin(tilt_angle * (M_PI / 180));
								d_v2 = -y1 * sin(tilt_angle * (M_PI / 180)) + y3 * cos(tilt_angle * (M_PI / 180));
								
								y1 = d_v1;
								y3 = d_v2;

								// z-Axis-rotation
								/*
								d_v1 = x1 * cos(tilt_angle * M_PI / 180.0) - x2 * sin(tilt_angle * M_PI / 180.0);
								d_v2 = x1 * sin(tilt_angle * M_PI / 180.0) + x2 * cos(tilt_angle * M_PI / 180.0);

								x1 = d_v1;
								x2 = d_v2;

								d_v1 = y1 * cos(tilt_angle*M_PI / 180.0) - y2 * sin(tilt_angle * M_PI / 180.0);
								d_v2 = y1 * sin(tilt_angle*M_PI / 180.0) + y2 * cos(tilt_angle * M_PI / 180.0);

								y1 = d_v1;
								y2 = d_v2;
								*/
								// rotate the vectors

								// Store Values into the arrays //
								kx1.push_back(x1);	// add (push) the values into the array 
								kx2.push_back(x2);	// add (push) the values into the array
								kx3.push_back(x3);	// add (push) the values into the array

								e.push_back(nrg);	// add (push) the energy into the array

								ky1.push_back(y1);	// add (push) the values into the array
								ky2.push_back(y2);	// add (push) the values into the array
								ky3.push_back(y3);	// add (push) the values into the array
								// Store Values into the arrays //
							}
							else
							{
								cout << "too less points of the trajectory (<3) for the particle #: "<< control1 <<" in the 'pe-trajectories.dat'" << endl;
								cout << "check the entry of the particle of the backscattered particle number: " << control1 << endl;
							}

							flag1 = 1;
							break;
						}

						cdsandenrgy[1][0] = cdsandenrgy[0][0];	// x-coordinate2  //
						cdsandenrgy[1][1] = cdsandenrgy[0][1];	// y-coordinate2  //
						cdsandenrgy[1][2] = cdsandenrgy[0][2];	// z-coordinate2  //
						cdsandenrgy[1][3] = cdsandenrgy[0][3];	// energy2        //
					}
					control1++;

				// GET BACKSCATTERED ENERGY OVER ANGULAR INFORMATION //
				}

				// "      2" -> trajectory of a transmitted particle
				if (temparray[4] == "      2")
				{
					// GET TRANSMITTED ENERGY OVER ANGULAR INFORMATION //
					//cout << "Trajektorie: " << temparray[0] << "; Exit: " << temparray[4]<< endl;
					control2++;
				}

				// "      3" -> trajectory of an absorbed particle found
				if (temparray[4] == "      3")
				{
					// GET ABSORBED INFORMATION //
					//cout << "absorbTrajektorie: " << temparray[0] << "; Exit: " << temparray[4]<< endl;
					control3++;
				}
			}
			else
			{
				getline(datei, zeile);
			}

			int_1++;

			if (int_1 == 100)  // don't cout every particles processed to speed it up (cout every particle would take much more time) //
			{
				cout << "\r"  << "Particle " <<  temparray[0];
				int_1 = 0;
			}
		}

		cout << "\r"  << "Particle " <<  temparray[0];

		cout << endl;
		cout << endl;

		////////////////////////////write values into the file (energy normalized over the highest energy)//////////////////////
		cout << "writing x_{1/2}, y_{1/2}, z_{1/2} and energy/highest into " << "'\033[40;0;33m" << file_out1 << "\033[0m' ...";

		for (int temp1 = 0; temp1 < control1; temp1++)
		{
			bckang << scientific << setw(13) << kx1[temp1] << " " << setw(13) << kx2[temp1] << " " << setw(13) << kx3[temp1] << " " << e[temp1] / energy_p << endl;
			bckang << scientific << setw(13) << ky1[temp1] << " " << setw(13) << ky2[temp1] << " " << setw(13) << ky3[temp1] << endl;
		}

		//////////////////////////// write values into the file (energy normalized over the highest energy) //////////////////////

		cout << endl;
		cout << "##########################################################" << endl;
		cout << temparray[0] << " particles processed" << endl;
		cout << endl;
		cout << "highest energy found : " << energy_p << " [eV]" << endl;
		cout << endl;
		cout << "Sum(Backscattered): " << control1  << " particles" << endl;
		cout << "Sum(Transmitted): "   << control2  << " particles" << endl;
		cout << "Sum(Absorbed): "      << control3  << " particles" << endl;
		cout << endl;
		cout << "Sum of all (Backscattered, Transmitted and Absorbed): " << control1 + control2 + control3 << " particles"<< endl;
		cout << "##########################################################"<< endl;

																					////////////////////////////////////////////////////////////////////
		from_string<float>(control_particle, std::string(temparray[0]), std::dec);	// convert the entry of the array temparray[0] into a float       //
																					// this entry gives the amount of the processed particles         //
																					// if this value differs from the sum of all patricles(backscat., //
																					// transmitted and absorbed particles something has gone wrong    //
																					////////////////////////////////////////////////////////////////////

		if ( (control1 + control2 + control3) - control_particle != 0)	// error -> we lost some particles somewhere //
		{
			cout << "error - some particles might have been lost!" << endl;
		}

		////////////////BINNING - Energy and Probability - Calculation of the spherical coordinates (sum into bins)/////////////////////
		cout <<  "\033[40;0;31mSpherical binning...  \033[0m" << endl;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// kx1[i] ... x-Coordinate1 of the vector of the exiting(backscattered) Particle                        //
		// kx2[i] ... y-Coordinate1 of the vector of the exiting(backscattered) Particle                        //
		// kx3[i] ... z-Coordinate1 of the vector of the exiting(backscattered) Particle                        //
		//                                                                                                      //
		//                                                                                                      //
		// ky1[i] ... x-Coordinate2 of the vector of the exiting(backscattered) Particle                        //
		// ky2[i] ... y-Coordinate2 of the vector of the exiting(backscattered) Particle                        //
		// ky3[i] ... z-Coordinate2 of the vector of the exiting(backscattered) Particle                        //
		//                                                                                                      //
		//                                                                                                      //
		// The e.g. x-Component of Vector of the exiting(backscattered) i particle is given with: ky1[i]-kx1[i] //
		//                                                                                                      //
		//   e[i] ... Energy of the backscattered Particle i in units of [eV]                                   //
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// control_particle = 500000;

		cout << "Sum of all incident particles set to: " << control_particle << endl;
		cout << "if SE are on - please change the variable control_particles in the source code!!" << endl;

									//////////////////////////////////////////////////////////////////////////////////////////////
		int intv_phi = 60.0;		// set interval amount of (360 / int_phi) bins ;     int_phi should be a natural number     //
		int intv_theta = 45.0;		// set interval amount of (180 / int_theta) bins                                            //
		double tan_safe = 0.0;		//                                                                                          //
		double cos_safe = 0.0;		//                                                                                          //
		double d_var1 = 0;			// checkvariable double 1                                                                   //
		int maxbins = 0;			// amount of bins overall                                                                   //
									//                                                                                          //
		int ccount = 0;				// global count variable                                                                    //
		int count_2 = 0;			// checkvariable for all processes particles gathered in bins and so on                     //
		int i_var1 = 0;				// checkvariable int 1                                                                      //
		int hits_count = 0;			// Countvariable for the hits (for hits over angular)                                       //
									//                                                                                          //
		double phi, theta;			// variables for angle phi and theta                                                        //
									//                                                                                          //
		float energyy = 0;			// stores the gathered energy in one bin                                                    //
									//////////////////////////////////////////////////////////////////////////////////////////////

		cout << "Number of bins:   theta [0°;180°]: " << intv_theta << " ; " << "phi [-180°;180°]: " << intv_phi << endl;
		cout << "Invervals:  theta: " << 180 / intv_theta << "° ; " << "phi: " << 360 / intv_phi << "°" << endl;
		cout << endl;

		maxbins = intv_theta * intv_phi;

		double hits_array[maxbins];		// stores the hits summarized in one bin temporarily //
		double energy_array[maxbins];	// stores energy summarized in one bin temporarily   //
		double theta_array[maxbins];	// stores theta angles                               //
		double phi_array[maxbins];		// stores phi angles                                 //

		bckprob << setw(15) << "phi" << " " << setw(15) << "theta" << " " << setw(15) <<  "hits" << endl;

		for (phi = -180.0; phi < 180.0; phi += (360.0 / intv_phi))
		{
			for (theta = 0; theta < 180.0; theta += (180 / intv_theta))
			{
				for (int temp2 = 0; temp2 <= control1; temp2++)
				{
					tan_safe = 0;
					cos_safe = 0;

					if ( (ky1[temp2] - kx1[temp2]) > 0 )	// distinction of cases for angel phi : 'x > 0' //
					{
						tan_safe = atan( ( ky2[temp2] - kx2[temp2] ) / ( ky1[temp2] - kx1[temp2] ) ) * (180.0 / M_PI);	// arcustan //
					}

					if ( (ky1[temp2] - kx1[temp2]) < 0 )	// distinction of cases for angel phi : 'x < 0' //
					{
						if( (ky2[temp2] - kx2[temp2]) < 0 )	// distinction of cases for angel phi : 'x < 0 and y < 0 '//
						{
							tan_safe = (atan( ( ky2[temp2] - kx2[temp2] ) / ( ky1[temp2] - kx1[temp2] ) ) - M_PI) * (180.0 / M_PI);	// arcustan //
						}


						if( (ky2[temp2] - kx2[temp2]) >= 0 )	// distinction of cases for angel phi : 'x < 0 and y >= 0' //
						{
							tan_safe = (atan( ( ky2[temp2] - kx2[temp2] ) / ( ky1[temp2] - kx1[temp2] ) ) + M_PI) * (180.0 / M_PI);	// arcustan //
						}
					}

					if ( (ky1[temp2] - kx1[temp2]) == 0 )	// distinction of cases for angel phi : 'x = 0' //
					{
						tan_safe = sgn(ky2[temp2] - kx2[temp2]) * (M_PI / 2) * (180.0 / M_PI);
					}

					if ( tan_safe > phi and tan_safe <= phi + (360.0 / intv_phi) )	// check if angle tan is inside the bin //
					{

						cos_safe = acos( ky3[temp2] - kx3[temp2] ) * (180.0 / M_PI);	// compute angle theta(='cos_safe') arcuscosinus //

						if ( cos_safe > theta and cos_safe <= theta + (180.0 / intv_theta) )	// check if spherical angle theta is inside the bin //
						{
							energyy += e[temp2];	// add energy to the bin because the vector corresponding to theta and phi is inside the bin       //
							ccount ++;				// increment the ccount to save the information of the processed (backscattered) particles         //
							hits_count++;			// increment hits_count for the hits added to THIS bin (gather probability over angle information) //
						}
					}
				}


				energy_array[count_2] = energyy;							// Store Energy of the bin (sum of all energy going 'through' this bin [eV]     //
				hits_array[count_2] = hits_count;							// Store Hits of the bin                                                        //
				theta_array[count_2] = theta + (180 / (2 * intv_theta));	// Store angle theta of the bin                                                 //
				phi_array[count_2] = phi + (360 / (2 * intv_phi));			// Store angle phi of the bin                                                   //

				count_2++;													// overall processed particles                                                  //

				cout << "\r"  << "Bins: " << count_2 << " / " << maxbins;

				energyy = 0;    // clear the gathered energy for the next bin //
				hits_count = 0; // clear the gathered hits for the next bin   //
			} 
		} 

		cout << endl;
		////////////////BINNING - Energy and Probability - Calculation of the spherical coordinates (sum into bins)/////////////////////



		////////////////////////////////////////////////////////////////////CUTTING/////////////////////////////////////////////////////////////////////////////////////
		// generating backscattered energy over probability with particles between phi:[-pi:pi] and theta:[theta_low:theta_high] -> stored in 'energy_probability_cut.dat'
		// all particles who don't cross this theta area are stored in 'energy_probability_cut_negative.dat'
		// all entries added from 'energy_probability_cut.dat' and 'energy_probability_cut_negative.dat' should result in the overall spectrum ('energy_probability.dat'

		e_cut = e;			// copy the vector containing all found energies of backscattered energies
		e_cut_negative = e;	// copy the vector containing all found energies of backscattered energies
		count_2 = 0;

		int theta_low = 0;		// lower angle we need the information of backscattered energy over probability
		int theta_high = 75;	// upper angle we need the information of backscattered energy over probability

		cout << "##########################################################"<< endl;

		cout << "collecting info(backscat. energy) about particles within the range [" << theta_low << "°:" << theta_high << "°] ...." << endl;

		for (phi = -180.0; phi < 180.0; phi += (360.0 / intv_phi) )
		{
			for (theta = 0; theta < 180; theta += (180 / intv_theta) )
			{
				for (int temp2 = 0; temp2 <= control1; temp2++)
				{
					tan_safe = 0;
					cos_safe = 0;

					if ( (ky1[temp2] - kx1[temp2]) > 0 ) // distinction of cases for angel phi : 'x > 0' //
					{
						tan_safe = atan( ( ky2[temp2] - kx2[temp2] ) / ( ky1[temp2] - kx1[temp2] ) ) * (180.0 / M_PI);   // arcustan //
					}

					if ( (ky1[temp2] - kx1[temp2]) < 0 )    // distinction of cases for angel phi : 'x < 0' //
					{
						if( (ky2[temp2] - kx2[temp2]) < 0 )  // distinction of cases for angel phi : 'x < 0 and y < 0 '//
						{
							tan_safe = (atan( ( ky2[temp2] - kx2[temp2] ) / ( ky1[temp2] - kx1[temp2] ) ) - M_PI) * (180.0 / M_PI);   // arcustan //
						}


						if( (ky2[temp2] - kx2[temp2]) >= 0 )  // distinction of cases for angel phi : 'x < 0 and y >= 0' //
						{
							tan_safe = (atan( ( ky2[temp2] - kx2[temp2] ) / ( ky1[temp2] - kx1[temp2] ) ) + M_PI) * (180.0 / M_PI);   // arcustan //
						}
					}

					if ( (ky1[temp2] - kx1[temp2]) == 0 ) // distinction of cases for angel phi : 'x = 0' //
					{
						tan_safe = sgn(ky2[temp2] - kx2[temp2]) * (M_PI / 2) * (180.0 / M_PI);
					}

					if ( tan_safe > phi and tan_safe <= phi + (360.0 / intv_phi) )  // check if angle tan is inside the bin //
					{
						cos_safe = acos( ky3[temp2] - kx3[temp2] ) * (180.0 / M_PI); // compute angle theta(='cos_safe') arcuscosinus //

						if ( cos_safe > theta and cos_safe <= theta + (180.0 / intv_theta) ) // check if spherical angle theta is inside the bin //
						{
							if (cos_safe > theta_low and cos_safe < theta_high)  // check if particle crosses the range we don't/want
							{
								e_cut_negative[temp2] = 0; //delete particle energy in negative list
							}
							else
							{
								e_cut[temp2] = 0; //delete particle energy in "normal" list
							}
						}
					}
				}

				count_2++;

				cout << "\r"  << "Bins: " << count_2 << " / " << maxbins;
			}
		} 

		cout << endl;
		cout << "done." << endl;

		int count_eliminated = 0;
		int count_eliminated_negative = 0;

		for (int temp2 = 0; temp2 <= control1; temp2++)
		{
			if (e_cut[temp2] != 0)
			{
				count_eliminated++;
			}
			if (e_cut_negative[temp2] != 0)
			{
				count_eliminated_negative++;
			}
		}

		if (control1 - (count_eliminated + count_eliminated_negative) > 0)
		{
			cout << "error; particles missing(cutting): " << control1 - (count_eliminated + count_eliminated_negative) << endl;
		}

		int energybins_cut = 1000;							// amount of energybins the cutting data(backscattered energy over probability) is collected//
		double probability_cut[energybins_cut];				// arrays the binned data is stored
		double probability_cut_negative[energybins_cut];	// array the binned(negative) data is stored

		cout <<  "\033[40;0;31mcollecting cutted data(backscat. energy over prob. into "<<  energybins_cut << " bins)...  \033[0m" << endl;

		for (int temp2 = 0; temp2 <= control1; temp2++)	// process all (backscattered) particles //
		{
			for (int temp3 = 0; temp3 < energybins_cut; temp3 ++)	// for loop as much bins there are //
			{
				if (e_cut[temp2] > temp3 * (energy_p / energybins_cut) and e_cut[temp2] <= (temp3 + 1) * (energy_p / energybins_cut))
				{
					probability_cut[temp3] ++;
				}

				if (e_cut_negative[temp2] > temp3 * (energy_p / energybins_cut) and e_cut_negative[temp2] <= (temp3 + 1) * (energy_p / energybins_cut))
				{
				probability_cut_negative[temp3] ++;
				}
			}

			cout << "\r"  << "particle: " << temp2;
		}

		cout << endl;

		energy_prob_cut << "# backscatter over prob. for particles inside phi:[-pi:pi], theta:[" << theta_low << ":" << theta_high << "]" << endl;
		energy_prob_cut_neg << "# backscatter over prob. for particles inside phi:[-pi:pi], theta:[" << theta_high << ":90]" << endl;

		cout <<  "\033[40;0;31mwriting data...  \033[0m" << endl;

		for (int temp3 = 0; temp3 < energybins_cut; temp3 ++)
		{
			energy_prob_cut << setw(15) << ((temp3 + 0.5) * (energy_p / energybins_cut)) << setw(15) << probability_cut[temp3] / (control_particle * (energy_p / energybins_cut)) << endl;
			energy_prob_cut_neg << setw(15) << ((temp3 + 0.5) * (energy_p / energybins_cut)) << setw(15) << probability_cut_negative[temp3] / (control_particle * (energy_p / energybins_cut)) << endl;
		}
		////////////////////////////////////////////////////////////////////CUTTING/////////////////////////////////////////////////////////////////////////////////////



		//////////////////////writing the (binned) angles and energies into the file - energies can be manipulated here/////////////////////////////////
		float solid_angle_bin = 0.0;		// solid angle of the actual bin (varies with the angle theta) multiplied with the amount of all particles(incident amount of particles)//
		float control_solid_angle = 0.0;	// variable in which the solid angle of every bin is summarized. This should be 4 * PI //

		bckpol << "#binned energy for each solid angle. Contains the energy going through the solid angle of the bin" << endl;
		bckpol << setw(15) << "phi" << " " << setw(15) << "theta" << " " << setw(15) << "energy" << endl;

		for (int t_var1 = 0; t_var1 < maxbins; t_var1++)
		{

			solid_angle_bin = control_particle * ((cos((theta_array[t_var1] - (180.0 / (2.0 * intv_theta))) * (M_PI / 180)) - cos((theta_array[t_var1] + (180.0 / (2.0 * intv_theta))) * (M_PI / 180))) * (((phi_array[t_var1] + (360.0 / (2.0 * intv_phi))) * (M_PI / 180)) - (phi_array[t_var1] - (360.0 / (2.0 * intv_phi))) * (M_PI / 180))); // solid angle of the bin multiplied with the total amount of particles //

			bckpol  <<  setw(15) << phi_array[t_var1] << " " << setw(15) << theta_array[t_var1] << " " << setw(15) << energy_array[t_var1] / (solid_angle_bin * energy_p) << endl;  // energy      //
			bckprob <<  setw(15) << phi_array[t_var1] << " " << setw(15) << theta_array[t_var1] << " " << setw(15) << hits_array[t_var1] / solid_angle_bin << endl;    // probability //

			control_solid_angle += solid_angle_bin;
		}

		if (fabs((control_solid_angle / control_particle) - (4 * M_PI)) > 0.01)					// absolute value of the difference of the sum of all bins(solid angle) and the 4*PI     //
		{																						// 0.01 is quite high. Usually (control_solid_angle/control_particle) is around 12.5664  //
			cout << "#1 sum of the solid angle of all bins differs to much from 4*PI" << endl;	// 0.01 is due to rounding. Sum of all bins will never be 4*PI :)                        //
			cout << "It should be " << 4 * M_PI << " but is " << control_solid_angle / control_particle << endl;
			cout << "size of the bins may be wrong!" << endl;
		}

		cout << ccount << " particles binned into " << maxbins << " bins. Written in the file: '\033[40;0;33m" << file_out2 << "\033[0m'" << endl;
		//////////////////////writing the (binned) angles and energies into the file - energies can be manipulated here/////////////////////////////////



		/////////////////////////////////////////sum energies into bins/////////////////////////////////////////////////////
		int energybins = 1000;				// amount of bins to collect the data: backscattered energy over probability //
		int particles_done = 0;				// particles processed                                                       //
		double energy_bins[energybins];		// stores phi angles                                                         //
		double probability[energybins];		// stores phi angles                                                         //

		double backscattered_energy = 0.0;

		cout << "##########################################################"<< endl;

		cout <<  "\033[40;0;31mCreating energy/probability data...  \033[0m" << endl;

		cout << "collecting energies(from 0 to " << energy_p << " eV) into " << energybins << " bins ..." << endl;

		for (int temp2 = 0; temp2 <= control1; temp2++)			// process all (backscattered) particles //
		{
			for (int temp3 = 0; temp3 < energybins; temp3 ++)	// for loop as much bins there are //
			{
				if (e[temp2] > temp3 * (energy_p / energybins) and e[temp2] <= (temp3 + 1) * (energy_p / energybins))
				{
					probability[temp3] ++;
					particles_done++;
				}
			}

			cout << "\r"  << "particle: " << temp2;

			backscattered_energy += e[temp2];
		}

		cout << endl;
		cout << "particles_done , control1  "<< particles_done << " , " << control1 << endl;

		if (particles_done != control1)
		{
			cout << "some particles are missing and were not processed ... :(" << endl;
		}
		/////////////////////////////////////////sum energies into bins/////////////////////////////////////////////////////



		//////////writing the energies and probabilites into the file - everything can be manipulated here (eg. normalized)////////////////
		cout << endl;
		cout << "writing file into: '\033[40;0;33m" << file_out4 << "\033[0m' ... " << endl;	/*  energy over probability data written here! */

		float control_fl_1 = 0.0;
		float control_fl_2 = 0.0;

		for (int temp3 = 0; temp3 < energybins; temp3 ++)
		{
			nrgprob << setw(15) << ((temp3 + 0.5)*(energy_p / energybins)) << setw(15) << probability[temp3] / (control_particle * (energy_p / energybins)) << endl;
			energy_prob_pure << ((temp3 + 0.5) * (energy_p / energybins)) << setw(15) << probability[temp3] << endl;
			control_fl_2 += probability[temp3] / (control_particle * energy_p);
			control_fl_1 += ((temp3 + 0.5) * (energy_p / energybins)) * probability[temp3];
		}

		if ( fabs(control_fl_2 - control1 / (control_particle * energy_p)) > 0.001)  // fabs due to rounding errors :/ //
		{
			cout << "error - fabs(control_fl_2 - control1 / (control1 * energy_p)) > 0.001 " << endl;
		}
		//////////writing the energies and probabilites into the file - everything can be manipulated here (eg. normalized)////////////////


/*
		//////////creating the 'raw-energy-data' -> all backscattered energies sorted from lowest to highest////////////////
		cout << "##########################################################"<< endl;

		cout <<  "\033[40;0;31mCreating raw-energy data...  \033[0m" << endl;

		float temp_energy;
		int flag = 0;
		int count_cycles = 1;

		cout << "arrange energies" << endl;

		while(flag == 0)  // do until all energies are arranged (from lowest to highest energy) //
		{
			flag = 1;

			for (int temp2 = 0; temp2 < control1; temp2++)	// process all (backscattered) particles //
			{
				if(e[temp2] > e[temp2+1])	// switch energies (sorting) //
				{
					temp_energy = e[temp2];
					e[temp2] = e[temp2+1];
					e[temp2 + 1] = temp_energy;
					flag = 0;
				}

			}

			cout << "\r" << count_cycles;
			count_cycles++;
		}

		for (int temp3 = 0; temp3 < control1; temp3 ++)	// write the data into the file //
		{
			rawenerg << e[temp3] << endl;
		}

		cout << endl;
		cout << "writing file into: '\033[40;0;33m" << file_out7 << "\033[0m' ... " << endl;
		//////////creating the 'raw-energy-data' -> all backscattered energies sorted from lowest to highest////////////////
*/


		//////////calculating/writing Backscattered probability/energy over angle for angle theta (angle phi is summarized)////////////////
		int t_var2;

		float control_solid_angle2 = 0.0;	// variable in which the solid angle of every bin is summarized. This should be 4*PI //
		float temp_energy_sum = 0.0;
		float temp_probability_sum = 0.0;

		std::vector<double> binned_energy, binned_probability;  // array to store the summarized probabilities/energies over the angle phi           //

		cout << "##########################################################"<< endl;
		cout <<  "\033[40;0;31mCreating energy/angle prob/angle summarized over phi...  \033[0m" << endl;

		for (int t_var1 = 0; t_var1 < intv_theta; t_var1++)
		{
			for (t_var2 = 0; t_var2 < intv_phi; t_var2++)
			{
				solid_angle_bin = intv_phi * control_particle * ((cos((theta_array[t_var1 + t_var2*intv_theta] - (180.0 / (2.0 * intv_theta))) * (M_PI / 180)) - cos((theta_array[t_var1 + t_var2 * intv_theta] + (180.0 / (2.0 * intv_theta))) * (M_PI / 180))) * (((phi_array[t_var1 + t_var2 * intv_theta] + (360.0 / (2.0 * intv_phi))) * (M_PI / 180)) - (phi_array[t_var1 + t_var2 * intv_theta] - (360.0 / (2.0 * intv_phi))) * (M_PI / 180)));	// solid angle of the bin (varies with angle theta) //

				control_solid_angle2 += solid_angle_bin;

				temp_probability_sum += hits_array[t_var1 + t_var2 * intv_theta] / solid_angle_bin;		// add probabilities for an fix angle theta for all phi-angles //
				temp_energy_sum += energy_array[t_var1 + t_var2 * intv_theta] / solid_angle_bin;		// add energies for an fix angle theta for all phi-angles //
			}

			binned_probability.push_back(temp_probability_sum);	// add (push) the binned probabilities summarized over angle phi into the array 'binned_probabilities' //
			binned_energy.push_back(temp_energy_sum);			// add (push) the binned energies summarized over angle phi into the array 'binned_energy' //

			temp_probability_sum = 0.0;
			temp_energy_sum = 0.0;
		}

		for (int ii = 0; ii < intv_theta; ii++)
		{
			prob_theta_binned << setw(15) << (ii + 0.5) * (180 / intv_theta) << setw(15) << binned_probability[ii] << endl;	// write summarized probabilities into file //
			energy_theta_binned << setw(15) << (ii + 0.5) * (180 / intv_theta) << setw(15) << binned_energy[ii] << endl;	// write summarized energies into file      //
		}

		if (fabs((control_solid_angle2 / (control_particle * intv_phi)) - (4 * M_PI)) > 0.01)		// absolute value of the difference of the sum of all bins(solid angle) and the 4 * PI   //
		{																							// 0.01 is quite high. Usually (control_solid_angle/control_particle) is around 12.5664  //
			cout << "#3 sum of the solid angle of all bins differs to much from 4*PI" << endl;		// 0.01 is due to rounding. Sum of all bins will never be 4 * PI :)                      //
			cout << "It should be " << 4 * M_PI << " but is " << control_solid_angle/control_particle << endl;
			cout << "size of the bins may be wrong!" << endl;
		}
		//////////calculating/writing Backscattered probability/energy over angle for angle theta (angle phi is summarized)////////////////


		//////////writing the infofile (containing the number of (all, backscat., transmit.)particles and probabilities////////////////
		cout <<  "\033[40;0;31mCreating information file...  \033[0m" << endl;

		cout << "writing general informations into: '\033[40;0;33m" << file_out6 << "\033[0m' ... " << endl;

		infofile << left << setw(50) << "File processed: " << file_in << endl;
		infofile << left << setw(50) << "Date of processing: " << ctime (&rawtime) << endl;
		infofile << endl;
		infofile << endl;
		infofile << left << setw(50) << "Sum of all particles: " << control1 + control2 + control3 << endl;
		infofile << left << setw(50) << "Sum of backscattered particles: " << control1 << endl;
		infofile << left << setw(50) << "Sum of transmitted particles: " << control2 << endl;
		infofile << left << setw(50) << "Sum of absorbed particles: " << control3 << endl;
		infofile << endl;
		infofile << left << setw(50) << "highest energy of a backscattered electron: " << energy_p << " eV" << endl;
		infofile << endl;
		infofile << left << setw(50) << "Amount of bins for the angle phi [-PI;PI]: " << intv_phi << endl;
		infofile << left << setw(50) << "Amount of bins for the angle theta [0;PI]: " << intv_theta << endl;
		infofile << endl;
		infofile << left << setw(50) << "Tilt-angle: " << tilt_angle << " °" << endl;
		infofile << endl;
		infofile << left << setw(50) << "(backscattered particles)/(Incident particles): " << ((float)control1) / (control1 + control2 + control3) << endl;
		infofile << left << setw(50) << "  (transmitted particles)/(Incident particles): " << ((float)control2) / (control1 + control2 + control3) << endl;
		infofile << left << setw(50) << "     (absorbed particles)/(Incident particles): " << ((float)control3) / (control1 + control2 + control3) << endl;
		infofile << endl;
		infofile << endl;
		infofile << left << setw(50) << "(backscattered energy)/(incident energy * sum of all particles): " << backscattered_energy / (control_particle * energy_p) << endl;

		cout << "##########################################################"<< endl;
		cout <<  "\033[40;0;34mFinished! \033[0m" << endl;
		}
		else // flag == 0 ... error reading the file //
		{
			cout << "error reading file '" <<  file_in << "'" << endl;
		}
		//////////writing the infofile (containing the number of (all, backscat., transmit.)particles and probabilities////////////////

		// closings
		datei.close();
		bckang.close();
		bckpol.close();
		bckprob.close();
		nrgprob.close();
	// 	bckprob2.close();
		infofile.close();
		rawenerg.close();
	// 	energyangle.close();
		prob_theta_binned.close();
		energy_theta_binned.close();
		energy_prob_pure.close();
		energy_prob_cut.close();
		energy_prob_cut_neg.close();
}
